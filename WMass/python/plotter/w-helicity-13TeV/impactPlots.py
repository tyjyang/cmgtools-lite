# USAGE: python impactPlots.py fitresults_poim1_exp0_bbb0.root -o plots/fit/absY/results/2018-12-10 --nuis 'CMS.*' --pois 'Wplus.*left.*'
import ROOT, os, datetime, re, operator, math
from array import array
from subMatrix import niceName
ROOT.gROOT.SetBatch(True)

from make_diff_xsec_cards import get_ieta_ipt_from_process_name


if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)

    from optparse import OptionParser
    parser = OptionParser(usage='%prog fitresults.root [options] ')
    parser.add_option('-o','--outdir',     dest='outdir',     default='',   type='string', help='outdput directory to save the matrix')
    parser.add_option(     '--pois',       dest='pois',       default='',   type='string', help='pois for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option(     '--nuis',       dest='nuis',       default='',   type='string', help='nuis for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option(     '--nuisgroups', dest='nuisgroups', default='',   type='string', help='nuis groups for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option(     '--suffix',     dest='suffix',     default='',   type='string', help='suffix for the correlation matrix')
    parser.add_option(     '--target',     dest='target',     default='mu', type='string', help='target POI (can be mu,xsec,xsecnorm)')
    parser.add_option(     '--zrange',     dest='zrange',     default='', type='string', help='Pass range for z axis as "min,max". If not given, it is set automatically based on the extremes of the histogram')
    parser.add_option(     '--canvasSize', dest='canvasSize', default='', type='string', help='Pass canvas dimensions as "height,width" ')
    parser.add_option(     '--margin',     dest='margin',     default='', type='string', help='Pass canvas margin as "left,right,top,bottom" ')
    parser.add_option(     '--nContours', dest='nContours',    default=0, type=int, help='Number of contours in palette. Default is 20')
    parser.add_option(     '--palette'  , dest='palette',      default=0, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_option(     '--parNameCanvas', dest='parNameCanvas',    default='', type='string', help='The canvas name is built using the parameters selected with --nuis or nuisgroups. If they are many, better to pass a name, like QCDscales or PDF for example')
    (options, args) = parser.parse_args()

    ROOT.TColor.CreateGradientColorTable(3,
                                      array ("d", [0.00, 0.50, 1.00]),
                                      ##array ("d", [1.00, 1.00, 0.00]),
                                      ##array ("d", [0.70, 1.00, 0.34]),
                                      ##array ("d", [0.00, 1.00, 0.82]),
                                      array ("d", [0.00, 1.00, 1.00]),
                                      array ("d", [0.34, 1.00, 0.65]),
                                      array ("d", [0.82, 1.00, 0.00]),
                                      255,  0.95)

    if len(options.nuis) and len(options.nuisgroups):
        print 'You can specify either single nuis (--nuis) or poi groups (--nuisgroups), not both!'
        sys.exit()

    if options.outdir:
        options.outdir = options.outdir + "/" + options.target + "/"
        ROOT.gROOT.SetBatch()
        if not os.path.isdir(options.outdir):
            os.system('mkdir -p {od}'.format(od=options.outdir))
        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=options.outdir))

    if len(options.nuis)==0 and len(options.nuisgroups)==0:
        print "Will plot the impact for all the single nuisances. It may be a big matrix"
        nuis_regexps = ['.*']
    else:
        nuis_regexps = list(options.nuis.split(',')) if len(options.nuis) else list(options.nuisgroups.split(','))
        print "Filtering nuisances or nuisance groups with the following regex: ",nuis_regexps

    if len(options.pois)==0:
        pois_regexps = ['.*']
    else:
        pois_regexps = list(options.pois.split(','))
    
    hessfile = ROOT.TFile(args[0],'read')
    
    group = 'group_' if len(options.nuisgroups) else ''
    if   options.target=='xsec':     target = 'pmaskedexp'
    elif options.target=='xsecnorm': target = 'pmaskedexpnorm'
    else:                            target = 'mu'

    th2name = 'nuisance_{group}impact_{sfx}'.format(group=group,sfx=target)
    impMat = hessfile.Get(th2name)
    if impMat==None:
        print "ERROR: Cannot find the impact TH2 named ",th2name," in the input file. Maybe you didn't run --doImpacts?\nSkipping."
        sys.exit()
    
    pois = []; pois_indices = []
    for ib in xrange(1,impMat.GetNbinsX()+1):
        for poi in pois_regexps:
            if re.match(poi, impMat.GetXaxis().GetBinLabel(ib)):
                pois.append(impMat.GetXaxis().GetBinLabel(ib))
                pois_indices.append(ib)

    nuisances = []; nuisances_indices = []
    for ib in xrange(1,impMat.GetNbinsY()+1):
        for n in nuis_regexps:
            if re.match(n, impMat.GetYaxis().GetBinLabel(ib)):
                nuisances.append(impMat.GetYaxis().GetBinLabel(ib))
                nuisances_indices.append(ib)

    mat = {}
    for ipoi, poi in enumerate(pois):
        for inuis, nuis in enumerate(nuisances):
            mat[(poi,nuis)] = impMat.GetBinContent(pois_indices[ipoi], nuisances_indices[inuis])

    ## sort the pois and nuisances alphabetically, except for pdfs, which are sorted by number
    pois = sorted(pois, key= lambda x: int(x.split('_')[-1]) if '_Ybin_' in x else 0)
    pois = sorted(pois, key= lambda x: get_ieta_ipt_from_process_name(x) if ('_ieta_' in x and '_ipt_' in x) else 0)
    if len(options.nuisgroups)==0:
        ## for mu* QCD scales, distinguish among muR and muRXX with XX in 1-10
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('pdf','')) if 'pdf' in x else 0)
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('muRmuF','')) if ('muRmuF' in x and x != "muRmuF")  else 0)
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('muR','')) if (''.join([j for j in x if not j.isdigit()]) == 'muR' and x != "muR") else 0)
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('muF','')) if (''.join([j for j in x if not j.isdigit()]) == 'muF' and x != "muF") else 0)
        nuisances = sorted(nuisances, key= lambda x: int(x.split('EffStat')[1]) if 'EffStat' in x else 0)

    #print "sorted pois = ", pois
    #print "\n\nsorted nuisances = ", nuisances

    ch = 1200
    cw = 1200
    if options.canvasSize:
        ch = int(options.canvasSize.split(',')[0])
        cw = int(options.canvasSize.split(',')[1])
    c = ROOT.TCanvas("c","",ch,cw)
    c.SetGridx()
    c.SetGridy()
    if options.nContours: ROOT.gStyle.SetNumberContours(options.nContours)
    if options.palette:   ROOT.gStyle.SetPalette(options.palette)

    clm = 0.15
    crm = 0.15
    cbm = 0.15
    ctm = 0.1
    if options.margin:
        clm,crm,ctm,cbm = (float(x) for x in options.margin.split(','))
    c.SetLeftMargin(clm)
    c.SetRightMargin(crm)
    c.SetBottomMargin(cbm)
    c.SetTopMargin(ctm)

    ## make the new, smaller TH2F correlation matrix
    nbinsx = len(pois); nbinsy = len(nuisances)
    th2_sub = ROOT.TH2F('sub_imp_matrix', '', nbinsx, 0., nbinsx, nbinsy, 0., nbinsy)
    th2_sub.GetXaxis().SetTickLength(0.)
    th2_sub.GetYaxis().SetTickLength(0.)

    ## pretty nested loop. enumerate the tuples
    for i,x in enumerate(pois):
        for j,y in enumerate(nuisances):
            ## set it into the new sub-matrix
            th2_sub.SetBinContent(i+1, j+1, mat[(x,y)])
            ## set the labels correctly
            new_x = niceName(x)
            new_y = niceName(y)
            th2_sub.GetXaxis().SetBinLabel(i+1, new_x)
            th2_sub.GetYaxis().SetBinLabel(j+1, new_y)

    rmax = max(abs(th2_sub.GetMaximum()),abs(th2_sub.GetMinimum()))    
    th2_sub.GetZaxis().SetRangeUser(-rmax,rmax)
    if options.zrange:
        rmin,rmax = (float(x) for x in options.zrange.split(','))
        th2_sub.GetZaxis().SetRangeUser(rmin,rmax)
    th2_sub.GetXaxis().LabelsOption("v")
    ROOT.gStyle.SetPaintTextFormat('1.3f')
    if options.target != "mu":
        ROOT.gStyle.SetPaintTextFormat('.3g')
    if len(pois)<30 and len(nuisances)<30: th2_sub.Draw('colz0 text45')
    else: th2_sub.Draw('colz0')
    if len(nuisances)>30: th2_sub.GetYaxis().SetLabelSize(0.02)

    if options.parNameCanvas: 
        cname = options.parNameCanvas
    else:
        nuisName = 'NuisGroup'+options.nuisgroups if len(options.nuisgroups) else options.nuis
        nuisName = nuisName.replace(',','AND').replace('.','').replace('*','').replace('$','').replace('^','').replace('|','').replace('[','').replace(']','')
        poisName = options.pois.replace(',','AND').replace('.','').replace('*','').replace('$','').replace('^','').replace('|','').replace('[','').replace(']','')
        cname = "{nn}_On_{pn}".format(nn=nuisName,pn=poisName)

    if options.outdir:
        for i in ['pdf', 'png']:
            suff = '' if not options.suffix else '_'+options.suffix
            c.SaveAs(options.outdir+'/smallImpacts{suff}_{target}_{cn}.{i}'.format(suff=suff,target=target,i=i,cn=cname))
        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=options.outdir))
