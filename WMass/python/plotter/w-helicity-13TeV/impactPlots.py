# USAGE: python impactPlots.py fitresults_poim1_exp0_bbb0.root -o plots/fit/absY/results/2018-12-10 --nuis 'CMS.*' --pois 'Wplus.*left.*'
import ROOT, os, datetime, re, operator, math
from array import array
from subMatrix import niceName
ROOT.gROOT.SetBatch(True)

from make_diff_xsec_cards import get_ieta_ipt_from_process_name
from make_diff_xsec_cards import getDiffXsecBinning
from make_diff_xsec_cards import templateBinning

from postFitPlots import prepareLegend

import utilities
utilities = utilities.util()


def niceSystName(label):
    if 'lepScale' in label: niceName = 'lepton scale'
    elif 'OtherBkg' in label: niceName = 'other bkg'
    elif 'pdfs' in label: niceName = 'PDFs #oplus #alpha_{S}'
    elif 'binByBinStat' in label: niceName = 'MC statistics'
    elif 'EffStat' in label: niceName = 'efficiency stat.'
    elif 'Fakes' in label: niceName = 'QCD bkg.'
    elif 'OtherExp' in label: niceName = 'other experimental'
    elif 'lumi' in label: niceName = 'luminosity'
    elif 'QCDTheo' in label: niceName = '#mu_{F}, #mu_{R}, #mu_{F}#mu_{R}'
    elif 'QEDTheo' in label: niceName = 'FSR'
    elif 'stat' in label: niceName = 'statistical'
    elif 'Total' in label: niceName = 'Total'
    elif 'EffSyst' in label: niceName = 'efficiency syst.'
    elif 'L1Prefire' in label: niceName = 'L1 prefire'
    else: niceName = label
    return niceName

def latexLabel(label):
    bin = int(label.split(' ')[-1]) if any(pol in label for pol in ['left','right','long']) else -1
    deltaYW = 0.2 ### should use binningYW.txt, but let's  not make it too complicated
    minY,maxY=bin*deltaYW,(bin+1)*deltaYW
    binLbl = ' {minY}<|Y_\mathrm{{ W }}|<{maxY}'.format(minY=minY,maxY=maxY)
    lblNoBin = ' '.join(label.split(' ')[:-1])
    lblNoBin = lblNoBin.replace('+','^+').replace(' left','_L').replace(' right','_R')
    lblNoBin += binLbl
    return lblNoBin

def prepareLegendV2(textSize=0.035,xmin=0.35,xmax=0.9 ):
    (x1,y1,x2,y2) = (xmin, 0.62, xmax, 0.87)
    leg = ROOT.TLegend(x1,y1,x2,y2)
    leg.SetNColumns(3)
    leg.SetFillColor(0)
    leg.SetFillColorAlpha(0,0.6)
    leg.SetShadowColor(0)
    leg.SetLineColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(textSize)
    return leg

if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)

    from optparse import OptionParser
    parser = OptionParser(usage='%prog fitresults.root [options] ')
    parser.add_option('-o','--outdir',     dest='outdir',     default='',   type='string', help='outdput directory to save the matrix')
    parser.add_option(     '--pois',       dest='pois',       default='',   type='string', help='pois for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option(     '--nuis',       dest='nuis',       default='',   type='string', help='nuis for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option(     '--nuisgroups', dest='nuisgroups', default='',   type='string', help='nuis groups for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option(     '--suffix',     dest='suffix',     default='',   type='string', help='suffix for the correlation matrix')
    parser.add_option(     '--target',     dest='target',     default='mu', type='string', help='target POI (can be mu,xsec,xsecnorm, or other things)')
    parser.add_option(     '--zrange',     dest='zrange',     default='', type='string', help='Pass range for z axis as "min,max". If not given, it is set automatically based on the extremes of the histogram')
    parser.add_option(     '--canvasSize', dest='canvasSize', default='', type='string', help='Pass canvas dimensions as "width,height" ')
    parser.add_option(     '--margin',     dest='margin',     default='', type='string', help='Pass canvas margin as "left,right,top,bottom" ')
    parser.add_option(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (let it be odd, so the central strip is white if not using --abs-value and the range is symmetric)')
    parser.add_option(     '--palette'  , dest='palette',      default=0, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_option(     '--invertPalette', dest='invertePalette' , default=False , action='store_true',   help='Inverte color ordering in palette')
    parser.add_option(     '--parNameCanvas', dest='parNameCanvas',    default='', type='string', help='The canvas name is built using the parameters selected with --nuis or nuisgroups. If they are many, better to pass a name, like QCDscales or PDF for example')
    parser.add_option(     '--abs-value', dest='absValue' , default=False , action='store_true',   help='Use absolute values for impacts (groups are already positive)')
    parser.add_option('-a','--absolute',   dest='absolute',   default=False, action='store_true', help='absolute uncertainty (default is relative)')
    parser.add_option(     '--latex',      dest='latex',      default=False, action='store_true', help='target POI (can be mu,xsec,xsecnorm)')
    parser.add_option('-y','--ybinfile',   dest='ybinfile',   default='',  type='string', help='do 1D summary plot as a function of YW using this file with the yw binning')
    parser.add_option(     '--etaptbinfile',   dest='etaptbinfile',   default='',  type='string', help='do 1D summary plot as a function of |eta| using this file with the |eta| binning')
    parser.add_option(     '--splitOutByTarget', dest='splitOutByTarget' , default=False , action='store_true',   help='Create a subfolder appending target to output directory, where plots will be saved')
    parser.add_option(     '--pt-min-signal'  , dest='ptMinSignal',  default=-1, type=float, help='Only for 2D xsec: specify the minimum pt for the bins considered as signal. If not given, take the very first pt value from the gen pt binning')
    parser.add_option('-c', '--channel',     dest='channel',     default='',   type='string', help='Specify channel (el|mu) to make legend. It is no longer guessed from the name of POIs. If empty, generic "l" will be used in legends instead of (#mu|e)')
    parser.add_option(     '--ybinsBkg', dest='ybinsBkg', type='string', default="10,11", help='Define which Y bins are to be considered as background. With format 14,15 ')
    parser.add_option(     '--longBkg'     , dest='longBkg'  , default=False         , action='store_true',   help='if True, longitudinal component was treated as background, so the POIs are missing. Manage inputs accordingly')
    parser.add_option(      '--skipPreliminary', dest='skipPreliminary', default=False, action='store_true', help='Do not add "Preliminary" to text on top left of canvas')
    (options, args) = parser.parse_args()

    # palettes:
    # 69 + inverted, using TColor::InvertPalette(), kBeach
    # 70 + inverted, using TColor::InvertPalette(), kBlackBody
    # 107, kCool
    # 57 kBird
    # 55 kRainBow

    ROOT.TColor.CreateGradientColorTable(3,
                                         array ("d", [0.00, 0.50, 1.00]),
                                         ##array ("d", [1.00, 1.00, 0.00]),
                                         ##array ("d", [0.70, 1.00, 0.34]),
                                         ##array ("d", [0.00, 1.00, 0.82]),
                                         array ("d", [0.00, 1.00, 1.00]),
                                         array ("d", [0.34, 1.00, 0.65]),
                                         array ("d", [0.82, 1.00, 0.00]),
                                         255,  0.95)

    if options.absValue:
        ROOT.TColor.CreateGradientColorTable(2,
                                             array ("d", [0.00, 1.00]),
                                             ##array ("d", [1.00, 1.00, 0.00]),
                                             ##array ("d", [0.70, 1.00, 0.34]),
                                             ##array ("d", [0.00, 1.00, 0.82]),
                                             array ("d", [1.00, 1.00]),
                                             array ("d", [1.00, 0.65]),
                                             array ("d", [1.00, 0.00]),
                                             255,  0.95)

    if len(options.nuis) and len(options.nuisgroups):
        print 'You can specify either single nuis (--nuis) or poi groups (--nuisgroups), not both!'
        sys.exit()

    if len(options.channel) and options.channel not in ["el","mu"]:
        print "Warning: if you specify a channel with -c, it must be 'el' or 'mu'"
        quit()

    # flavour is used in legend, so it is a TLatex string
    # for generic lepton, the Latex \ell is not supported (once should use TMathText which allows to write in Latex style and has \ell, but
    # apparently it is not supported for PDFs (one solution is that one should save in .eps and convert the file to pdf format later)
    # see: https://root-forum.cern.ch/t/tlatex-vs-ell/7767/11
    flavour = "l" 
    if len(options.channel):
        if options.channel == "el": flavour = "e"
        else                      : flavour = "#mu"

    if options.outdir:
        # emanuele: why this?
        #options.outdir = options.outdir + "/" + options.target + "/"
        # mciprian: for better bookkeeping, because it was easier to look at many plots
        if options.splitOutByTarget:
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
    valuesAndErrors = utilities.getFromHessian(args[0], params = pois_regexps)

    group = 'group_' if len(options.nuisgroups) else ''
    if   options.target=='xsec':           target = 'pmaskedexp'
    elif options.target=='xsecnorm':       target = 'pmaskedexpnorm'
    elif options.target=='unpolxsec':      target = 'sumpois'
    elif options.target=='unpolxsecnorm':  target = 'sumpoisnorm'
    elif options.target=='asym':           target = 'chargepois'
    elif options.target=='unpolasym':      target = 'chargemetapois'
    elif options.target=='A0' or options.target=='A4': target = 'polpois'
    elif options.target=='etaptasym':  target = 'chargepois'
    elif any(options.target==x for x in ['etaxsec','ptxsec']):          target = 'sumpois'
    elif any(options.target==x for x in ['etaxsecnorm', 'ptxsecnorm']): target = 'sumpoisnorm'
    elif any(options.target==x for x in ['etaasym', 'ptasym']):         target = 'chargemetapois'
    else:                              target = 'mu'

    # patch for diff.xsec, because I usually pass pois as regular expressions matching the charge
    if any(options.target == x for x in ['etaasym','etaptasym']):
        pois_regexps = list(x.replace('plus','').replace('minus','') for x in pois_regexps)    

    #if options.ybinfile:
    #    pois_regexps = ['W{ch}.*{pol}.*'.format(ch=charge,pol=pol) for charge in ['plus','minus'] for pol in ['left','right','long']]

    th2name = 'nuisance_{group}impact_{sfx}'.format(group=group,sfx=target)
    impMat = hessfile.Get(th2name)
    if impMat==None:
        print "ERROR: Cannot find the impact TH2 named ",th2name," in the input file. Maybe you didn't run --doImpacts?\nSkipping."
        sys.exit()
    
    # apparently, target chargemetapois for diff.xsec contains labels ending with chargemetaasym and chargemetatotalxsec, 
    # so one should filter one of them according to options.target
    # same for chargepois, that has both chargeasym and chargetotalxsec

    # use dictionary for pois_indices, so to keep track of poi's name and allow for some sorting (for example to allow plus and minus in same plot)
    pois = []; pois_indices = {}
    for ib in xrange(1,impMat.GetNbinsX()+1):
        for poi in pois_regexps:
            if re.match(poi, impMat.GetXaxis().GetBinLabel(ib)):
                if options.target=='etaasym':
                    if re.match(".*chargemetaasym", impMat.GetXaxis().GetBinLabel(ib)):
                        pois.append(impMat.GetXaxis().GetBinLabel(ib))
                        pois_indices[impMat.GetXaxis().GetBinLabel(ib)] = ib
                elif options.target=='etaptasym':
                    if re.match(".*chargeasym", impMat.GetXaxis().GetBinLabel(ib)):
                        pois.append(impMat.GetXaxis().GetBinLabel(ib))
                        pois_indices[impMat.GetXaxis().GetBinLabel(ib)] = ib
                else:
                    pois.append(impMat.GetXaxis().GetBinLabel(ib))
                    pois_indices[impMat.GetXaxis().GetBinLabel(ib)] = ib


    nuisances = []; nuisances_indices = []
    for ib in xrange(1,impMat.GetNbinsY()+1):
        for n in nuis_regexps:
            if re.match(n, impMat.GetYaxis().GetBinLabel(ib)):
                nuisances.append(impMat.GetYaxis().GetBinLabel(ib))
                nuisances_indices.append(ib)

    mat = {}

    pois = sorted(pois, key= lambda x: int(x.split('_')[-2]) if '_Ybin_' in x else 0)
    pois = sorted(pois, key= lambda x: get_ieta_ipt_from_process_name(x) if ('_ieta_' in x and '_ipt_' in x) else 0)
    # sort by charge if needed
    pois = sorted(pois, key = lambda x: 0 if "plus" in x else 1 if "minus" in x else 2)

    for ipoi, poi in enumerate(pois):
        for inuis, nuis in enumerate(nuisances):
            mat[(poi,nuis)] = impMat.GetBinContent(pois_indices[poi], nuisances_indices[inuis])

    ## sort the pois and nuisances alphabetically, except for pdfs, which are sorted by number
    if len(options.nuisgroups)==0:
        ## for mu* QCD scales, distinguish among muR and muRXX with XX in 1-10
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('pdf','')) if 'pdf' in x else 0)
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('muRmuF','')) if ('muRmuF' in x and x != "muRmuF")  else 0)
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('muR','')) if (''.join([j for j in x if not j.isdigit()]) == 'muR' and x != "muR") else 0)
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('muF','')) if (''.join([j for j in x if not j.isdigit()]) == 'muF' and x != "muF") else 0)
        nuisances = sorted(nuisances, key= lambda x: utilities.getNEffStat(x) if 'EffStat' in x else 0)

    #print "sorted pois = ", pois
    #print "\n\nsorted nuisances = ", nuisances

    ch = 1200
    cw = 1200
    if options.canvasSize:
        cw = int(options.canvasSize.split(',')[0])
        ch = int(options.canvasSize.split(',')[1])
    c = ROOT.TCanvas("c","",cw,ch)
    c.SetGridx()
    c.SetGridy()
    if options.nContours: ROOT.gStyle.SetNumberContours(options.nContours)
    if options.palette:   ROOT.gStyle.SetPalette(options.palette)
    if options.invertePalette:    ROOT.TColor.InvertPalette()


    clm = 0.2
    crm = 0.2
    cbm = 0.2
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
    
    poiName_target = {"mu":        "signal strength",
                      "xsec":      "cross section",
                      "xsecnorm":  "normalized cross section",
                      "unpolxsec": "unpolarized cross section",
                      "unpolxsecnorm": "unpolarized normalized cross section",
                      "asym":      "charge asymmetry",
                      "unpolasym": "unpolarized charge asymmetry",
                      "A0":        "A0",
                      "A4":        "A4",
                      "etaptasym":   "charge asymmetry",
                      "etaasym":   "charge asymmetry",
                      "etaxsec":   "cross section",
                      "etaxsecnorm": "normalized cross section",
                      "ptasym":   "charge asymmetry",
                      "ptxsec":   "cross section",
                      "ptxsecnorm": "normalized cross section",
                      }
    th2_sub.GetZaxis().SetTitle("impact on POI for {p} {units}".format(units='' if options.absolute else '(%)', p=poiName_target[options.target]))
    th2_sub.GetZaxis().SetTitleOffset(1.2)

    ## pretty nested loop. enumerate the tuples
    for i,x in enumerate(pois):
        for j,y in enumerate(nuisances):
            ## set it into the new sub-matrix
            if options.absolute: 
                val = mat[(x,y)]
            else: 
                val = 100*mat[(x,y)]/valuesAndErrors[x][0] if valuesAndErrors[x][0] !=0 else 0.0
            th2_sub.SetBinContent(i+1, j+1, abs(val) if options.absValue else val)
            ## set the labels correctly
            new_x = niceName(x)
            new_y = niceName(y)
            th2_sub.GetXaxis().SetBinLabel(i+1, new_x)
            th2_sub.GetYaxis().SetBinLabel(j+1, new_y)
            

    rmax = max(abs(th2_sub.GetMaximum()),abs(th2_sub.GetMinimum()))
    if not options.absolute: rmax = min(20.,rmax)
    if options.absValue:
        th2_sub.GetZaxis().SetRangeUser(0,rmax)
    else :
        th2_sub.GetZaxis().SetRangeUser(-rmax,rmax)
    if options.zrange:
        rmin,rmax = (float(x) for x in options.zrange.split(','))
        th2_sub.GetZaxis().SetRangeUser(rmin,rmax)
    th2_sub.GetXaxis().LabelsOption("v")
    if options.absolute: 
        ROOT.gStyle.SetPaintTextFormat('1.3f')
        if options.target != "mu":
            ROOT.gStyle.SetPaintTextFormat('.3g')
    else: 
        ROOT.gStyle.SetPaintTextFormat('1.1f')

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

    suff = '' if not options.suffix else '_'+options.suffix
    poisNameNice = options.parNameCanvas if len(options.parNameCanvas) else  poisName.replace('(','').replace(')','').replace('\|','or')
    if options.latex:
        txtfilename = 'smallImpacts{rel}{suff}_{target}_{nn}_On_{pn}.tex'.format(rel='Abs' if options.absolute else 'Rel',suff=suff,target=options.target,i=i,nn=nuisName,pn=poisNameNice)
        if options.outdir: txtfilename = options.outdir + '/' + txtfilename
        txtfile = open(txtfilename,'w')
        txtfile.write("\\begin{tabular}{l "+"   ".join(['r' for i in xrange(th2_sub.GetNbinsX())])+"} \\hline \n")
        txtfile.write('                    ' + " & ".join(['{poi}'.format(poi=latexLabel(th2_sub.GetXaxis().GetBinLabel(i+1))) for i in xrange(th2_sub.GetNbinsX())]) + "\\\\ \\hline \n")
        for j in xrange(th2_sub.GetNbinsY()):
            txtfile.write('{label:<20}& '.format(label=th2_sub.GetYaxis().GetBinLabel(j+1)) + " & ".join(['{syst:^.2f}'.format(syst=th2_sub.GetBinContent(i+1,j+1)) for i in xrange(th2_sub.GetNbinsX())]) + "\\\\ \n")
        txtfile.write("\\end{tabular}\n")
        txtfile.close()
        print "Saved the latex table in file ",txtfilename

    if options.outdir and not options.latex and not (options.ybinfile or options.etaptbinfile):
        for i in ['pdf', 'png']:
            suff = '' if not options.suffix else '_'+options.suffix
            c.SaveAs(options.outdir+'/smallImpacts{rel}{suff}_{target}_{cn}.{i}'.format(rel='Abs' if options.absolute else 'Rel',suff=suff,target=options.target,i=i,cn=cname))

        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=options.outdir))

    if options.ybinfile:
        ybinfile = options.ybinfile
        ybinfile = open(ybinfile, 'r')
        ybins = eval(ybinfile.read())
        ybinfile.close()

        bkgYBins = []
        if options.ybinsBkg:
            bkgYBins = list(int(i) for i in options.ybinsBkg.split(','))        

        summaries = {}
        groups = [th2_sub.GetYaxis().GetBinLabel(j+1) for j in xrange(th2_sub.GetNbinsY())]
        charges = ['allcharges'] if 'asym' in options.target else ['plus','minus']
        polarizations = ['unpolarized'] if re.match('^unpol|^A\d',options.target) else ['left','right','long']
        cp = 'plus_left' # groups assume a common Y binning
        for charge in charges:
            for pol in polarizations:
                ing  = 0
                ing2 = 0
                for ing_tmp,nuisgroup in enumerate(groups):
                    # in case we add other lines, let's avoid messing up all the colors of the old ones
                    # I define a new integer counter that is not updated on certain conditions
                    h = ROOT.TH1D(charge+'_'+pol+'_'+nuisgroup,'',len(ybins[cp])-1,array('d',ybins[cp]))
                    summaries[(charge,pol,nuisgroup)] = h
                    summaries[(charge,pol,nuisgroup)].SetMarkerSize(2)
                    if nuisgroup in ["QCDTheo","L1Prefire"]:
                        summaries[(charge,pol,nuisgroup)].SetMarkerColor(utilities.safecolor(len(groups)+ing2+1))
                        summaries[(charge,pol,nuisgroup)].SetLineColor(utilities.safecolor(len(groups)+ing2+1))
                        ing2 += 1
                    else:
                        summaries[(charge,pol,nuisgroup)].SetMarkerColor(utilities.safecolor(ing+1))
                        summaries[(charge,pol,nuisgroup)].SetLineColor(utilities.safecolor(ing+1))
                        ing += 1
                    summaries[(charge,pol,nuisgroup)].SetLineWidth(2)
                    summaries[(charge,pol,nuisgroup)].GetXaxis().SetRangeUser(0.,2.41)
                    summaries[(charge,pol,nuisgroup)].GetXaxis().SetTitle('|Y_{W}|')
                    if options.absolute:
                        summaries[(charge,pol,nuisgroup)].GetYaxis().SetRangeUser(5.e-4,1.)
                        summaries[(charge,pol,nuisgroup)].GetYaxis().SetTitle('Uncertainty')
                    else:
                        summaries[(charge,pol,nuisgroup)].GetYaxis().SetRangeUser(5.e-3,500.)
                        summaries[(charge,pol,nuisgroup)].GetYaxis().SetTitle('Relative uncertainty (%)')
                    summaries[(charge,pol,nuisgroup)].GetXaxis().SetTitleSize(0.06)
                    summaries[(charge,pol,nuisgroup)].GetXaxis().SetLabelSize(0.04)
                    summaries[(charge,pol,nuisgroup)].GetYaxis().SetTitleSize(0.06)
                    summaries[(charge,pol,nuisgroup)].GetYaxis().SetLabelSize(0.04)
                    summaries[(charge,pol,nuisgroup)].GetYaxis().SetTitleOffset(0.7)
        for j,nuisgroup in enumerate(groups):
            for i in xrange(th2_sub.GetNbinsX()):
                lbl = th2_sub.GetXaxis().GetBinLabel(i+1)
                if '+' in lbl: charge='plus'
                elif '-' in lbl: charge='minus'
                else: charge='allcharges' 
                pol = 'unpolarized' if len(polarizations)==1 else lbl.split()[-2]
                ybin = int(lbl.split()[-1])
                summaries[(charge,pol,nuisgroup)].SetBinContent(ybin+1,th2_sub.GetBinContent(i+1,j+1))

        fitErrors = {} # this is needed to have the error associated to the TH2 x axis label
        for i,x in enumerate(pois):
            new_x = niceName(x).replace('el: ','').replace('#mu: ','').replace('l: ','')
            if options.absolute:
                fitErrors[new_x] = abs(valuesAndErrors[x][1]-valuesAndErrors[x][0])
            else:
                fitErrors[new_x] = 100*abs((valuesAndErrors[x][1]-valuesAndErrors[x][0])/valuesAndErrors[x][0])

        etaptbinfile = options.ybinfile.replace("binningYW.txt","binningPtEta.txt")
        etaPtBinningVec = getDiffXsecBinning(etaptbinfile, "reco")
        recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])
        ptRangeText = "p_{T}^{%s} #in [%.3g, %.3g] GeV" % (flavour, recoBins.ptBins[0], recoBins.ptBins[-1])

        cs = ROOT.TCanvas("cs","",1800,900)
        cs.SetLeftMargin(0.1)
        cs.SetRightMargin(0.05)
        cs.SetBottomMargin(0.15)
        cs.SetTopMargin(0.1)
        cs.SetTicky(1)
        cs.SetGridy()
        for charge in charges:
            for ipol,pol in enumerate(polarizations):
                #leg = prepareLegend(xmin=0.5,legWidth=0.40,textSize=0.04)
                leg = prepareLegendV2(textSize=0.04,xmin=0.12,xmax=0.94)  # similar to function above, but should use the space is in a slightly better way
                if charge=='allcharges': sign=''
                else: sign='+' if charge is 'plus' else '-' 
                # leave a space before - sign, otherwise - is too close to W (in the png it gets too far instead, ROOT magic!)
                thischannel = "W_{{{pol}}}^{{{chsign}}} #rightarrow {fl}#nu".format(pol=pol,chsign=sign if sign != "-" else (" "+sign),fl=flavour)
                header = "#bf{{Uncertainties on {p} for {ch}     {ptt}}}".format(p=poiName_target[options.target], ch=thischannel, ptt=ptRangeText)
                leg.SetHeader(header)            
                leg.SetNColumns(4)
                lat = ROOT.TLatex()
                lat.SetNDC(); lat.SetTextFont(42)
                quadrsum = summaries[(charge,pol,groups[0])].Clone('quadrsum_{ch}_{pol}'.format(ch=charge,pol=pol))
                totalerr = quadrsum.Clone('totalerr_{ch}_{pol}'.format(ch=charge,pol=pol))
                for ing,ng in enumerate(groups):
                    drawopt = 'pl' if ing==0 else 'pl same'
                    summaries[(charge,pol,ng)].Draw(drawopt)
                    if   ng=='binByBinStat': summaries[(charge,pol,ng)].SetMarkerStyle(ROOT.kFullCircle)
                    elif ng=='stat'        : summaries[(charge,pol,ng)].SetMarkerStyle(ROOT.kFullCircle); summaries[(charge,pol,ng)].SetMarkerColor(ROOT.kBlack); summaries[(charge,pol,ng)].SetLineColor(ROOT.kBlack);
                    elif ng=='luminosity'  : summaries[(charge,pol,ng)].SetMarkerStyle(ROOT.kFullSquare)
                    elif ng=='QEDTheo'     : summaries[(charge,pol,ng)].SetMarkerStyle(ROOT.kOpenSquareDiagonal)
                    elif ng=='L1Prefire'     : summaries[(charge,pol,ng)].SetMarkerStyle(ROOT.kFourSquaresPlus)
                    else: summaries[(charge,pol,ng)].SetMarkerStyle(ROOT.kFullTriangleUp+ing)
                    leg.AddEntry(summaries[(charge,pol,ng)], niceSystName(ng), 'pl')
                    # now compute the quadrature sum of all the uncertainties (neglecting correlations among nuis groups)
                    for y in xrange(quadrsum.GetNbinsX()):
                        quadrsum.SetBinContent(y+1,math.hypot(quadrsum.GetBinContent(y+1),summaries[(charge,pol,ng)].GetBinContent(y+1)))
                quadrsum.SetMarkerStyle(ROOT.kFullCrossX); quadrsum.SetMarkerSize(3); quadrsum.SetMarkerColor(ROOT.kBlack); quadrsum.SetLineColor(ROOT.kBlack); quadrsum.SetLineStyle(ROOT.kDashed)
                ## this is only for debuggin purposes, do not show for the final plots
                # quadrsum.Draw('pl same')
                # now fill the real total error from the fit (with correct correlations)
                for y in xrange(totalerr.GetNbinsX()):
                    if y in bkgYBins: continue
                    if pol=='long' and options.longBkg: continue
                    totalerr.SetBinContent(y+1,fitErrors['W{sign} {pol} {bin}'.format(sign=sign,pol=pol,bin=y)])
                totalerr.SetMarkerStyle(ROOT.kFullDoubleDiamond); totalerr.SetMarkerSize(3); totalerr.SetMarkerColor(ROOT.kRed+1); totalerr.SetLineColor(ROOT.kRed+1);
                totalerr.Draw('pl same')
                #leg.AddEntry(quadrsum, 'Quadr. sum. impacts', 'pl')
                leg.AddEntry(totalerr, 'Total uncertainty', 'pl')
                leg.Draw('same')
                lat.DrawLatex(0.1, 0.92, '#bf{CMS}' + ('' if options.skipPreliminary else ' #it{Preliminary}'))
                lat.DrawLatex(0.78, 0.92, '35.9 fb^{-1} (13 TeV)')
                for i in ['pdf', 'png', '.C']:
                    suff = '' if not options.suffix else '_'+options.suffix
                    cs.SetLogy()
                    cs.SaveAs(options.outdir+'/ywImpacts{rel}{suff}_{target}_{ch}{pol}.{i}'.format(rel='Abs' if options.absolute else 'Rel',suff=suff,target=options.target,i=i,ch=charge,pol=pol))



    if options.etaptbinfile:
        
        etaPtBinningVec = getDiffXsecBinning(options.etaptbinfile, "gen")
        genBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])

        summaries = {}
        groups = [th2_sub.GetYaxis().GetBinLabel(j+1) for j in xrange(th2_sub.GetNbinsY())]
        charges = ['allcharges'] if 'asym' in options.target else ['plus','minus']
        for charge in charges:
            ing = 0
            ing2 = 0
            for ing_tmp,nuisgroup in enumerate(groups):
                # in case we add other lines, let's avoid messing up all the colors of the old ones
                # I define a new integer counter that is not updated on certain conditions
                h = None
                if options.target in ["ptxsec", "ptxsecnorm", "ptasym"]:
                    # ptBinsInRange = []
                    # for x in genBins.ptBins:
                    #     if x >= options.ptMinSignal: ptBinsInRange.append(float(x))
                    # print ptBinsInRange
                    # h = ROOT.TH1D(charge+'_'+nuisgroup,'',len(ptBinsInRange),array('d',ptBinsInRange))
                    h = ROOT.TH1D(charge+'_'+nuisgroup,'',genBins.Npt,array('d',genBins.ptBins))
                else:
                    h = ROOT.TH1D(charge+'_'+nuisgroup,'',genBins.Neta,array('d',genBins.etaBins))
                #print "Last bin right edge: " + str(h.GetXaxis().GetBinUpEdge(h.GetNbinsX()))
                summaries[(charge,nuisgroup)] = h
                summaries[(charge,nuisgroup)].SetMarkerSize(2)
                if nuisgroup in ["QEDTheo","L1Prefire"]:
                    summaries[(charge,nuisgroup)].SetMarkerColor(utilities.safecolor(len(groups)+ing2+1))
                    summaries[(charge,nuisgroup)].SetLineColor(utilities.safecolor(len(groups)+ing2+1))
                    ing2 += 1
                else:
                    summaries[(charge,nuisgroup)].SetMarkerColor(utilities.safecolor(ing+1))
                    summaries[(charge,nuisgroup)].SetLineColor(utilities.safecolor(ing+1))
                    ing += 1
                summaries[(charge,nuisgroup)].SetLineWidth(2)
                summaries[(charge,nuisgroup)].GetXaxis().SetRangeUser(0.,2.4)
                summaries[(charge,nuisgroup)].GetXaxis().SetTitle('lepton |#eta|')
                if options.target in ["ptxsec", "ptxsecnorm", "ptasym"]:
                    summaries[(charge,nuisgroup)].GetXaxis().SetTitle('lepton p_{T} [GeV]')
                if options.absolute:
                    summaries[(charge,nuisgroup)].GetYaxis().SetRangeUser(5.e-4,1.)
                    summaries[(charge,nuisgroup)].GetYaxis().SetTitle('Uncertainty')
                else:
                    summaries[(charge,nuisgroup)].GetYaxis().SetRangeUser(1.e-3,500.)
                    summaries[(charge,nuisgroup)].GetYaxis().SetTitle('Relative uncertainty (%)')
                summaries[(charge,nuisgroup)].GetXaxis().SetTitleSize(0.06)
                summaries[(charge,nuisgroup)].GetXaxis().SetLabelSize(0.04)
                summaries[(charge,nuisgroup)].GetYaxis().SetTitleSize(0.06)
                summaries[(charge,nuisgroup)].GetYaxis().SetLabelSize(0.04)
                summaries[(charge,nuisgroup)].GetYaxis().SetTitleOffset(0.7)
        for j,nuisgroup in enumerate(groups):
            for i in xrange(th2_sub.GetNbinsX()):
                lbl = th2_sub.GetXaxis().GetBinLabel(i+1)
                # we expect something like "el: W+ i#eta, ip_T = 3, 5" we need to get 3 
                etabin = int(((lbl.split('=')[-1]).split(',')[0]).rstrip().lstrip())
                #print "lbl = %s   bineta = %s" % (lbl, str(etabin))
                if 'W+' in lbl: charge='plus'
                elif 'W-' in lbl: charge='minus'
                else: charge='allcharges' 
                if options.target in ["ptxsec", "ptxsecnorm", "ptasym"]:
                    if options.ptMinSignal > 0 and genBins.ptBins[int(etabin)] < options.ptMinSignal:
                        summaries[(charge,nuisgroup)].SetBinContent(etabin+1,0.0)
                    else:
                        summaries[(charge,nuisgroup)].SetBinContent(etabin+1,th2_sub.GetBinContent(i+1,j+1))
                else:
                    summaries[(charge,nuisgroup)].SetBinContent(etabin+1,th2_sub.GetBinContent(i+1,j+1))

        fitErrors = {} # this is needed to have the error associated to the TH2 x axis label
        for i,x in enumerate(pois):
            #print "Debug: poi = " + x
            ivarMatch = "_ieta_" if options.target not in ["ptxsec", "ptxsecnorm", "ptasym"] else "_ipt_"
            bineta = (x.split(ivarMatch)[1]).split("_")[0]
            new_x = "W{ch} {bin}".format(ch="+" if "plus" in x else "-" if "minus" in x else "", bin=bineta)
            # print "new_x = " + new_x
            if options.target in ["ptxsec", "ptxsecnorm", "ptasym"]:
                if options.ptMinSignal > 0 and genBins.ptBins[int(bineta)] < options.ptMinSignal:
                    continue
            if options.absolute:
                fitErrors[new_x] = abs(valuesAndErrors[x][1]-valuesAndErrors[x][0])
            else:
                fitErrors[new_x] = 100*abs((valuesAndErrors[x][1]-valuesAndErrors[x][0])/valuesAndErrors[x][0])

        ptmin = genBins.ptBins[0] if options.ptMinSignal < 0 else options.ptMinSignal
        ptmax = genBins.ptBins[-1]
        etamin = genBins.etaBins[0]
        etamax = genBins.etaBins[-1]
        # this will be for strips at constant pt (if the proper options were passed)
        isSinglePtStrip = False
        if options.target not in ["etaxsec", "etaxsecnorm", "etaasym", "ptxsec", "ptxsecnorm", "ptasym"]:
            print "Going to draw a strip at constant pt"
            isSinglePtStrip = True
            ptbin_hasChanged = False
            theptbin = 0
            for i in xrange(th2_sub.GetNbinsX()):
                lbl = th2_sub.GetXaxis().GetBinLabel(i+1)
                # we expect something like "el: W+ i#eta, ip_T = 3, 5" we need to get 5 
                ptbin = int(((lbl.split('=')[-1]).split(',')[1]).rstrip().lstrip())
                if i == 0: 
                    theptbin = ptbin  # for the first item, assign the value
                else: 
                    if theptbin != ptbin:
                        ptbin_hasChanged = True                        
                        theptbin = ptbin
                #print "lbl = %s   binpt = %s" % (lbl, str(ptbin))
            if ptbin_hasChanged:
                print "Warning: based on the options, you wanted to select a strip along eta at constant pt"
                print "However, I see that you are selecting more pt bins, check the option for the regular expression to select POIs"
                quit()
            #print "theptbin = " + str(theptbin)
            ptmin = genBins.ptBins[theptbin]
            ptmax = genBins.ptBins[theptbin+1]
            print "Plot for pt-bin %d: ptmin, ptmax = %.3g, %.3g " % (theptbin, ptmin, ptmax)
        ptRangeText = "p_{T}^{%s} #in [%.3g, %.3g] GeV" % (flavour, ptmin, ptmax)
        etaRangeText = "| #eta^{%s} | #in [%.3g, %.3g]" % (flavour, etamin, etamax)
        if etamin == 0.0:
            etaRangeText = "| #eta^{%s} | < %.3g" % (flavour, etamax)            
        #print ptRangeText
        #for key in fitErrors:
        #    print "fitErrors: key = " + key

        cs = ROOT.TCanvas("cs","",1800,900)
        cs.SetLeftMargin(0.1)
        cs.SetRightMargin(0.05)
        cs.SetBottomMargin(0.15)
        cs.SetTopMargin(0.1)
        cs.SetGridy()
        cs.SetTicky(1)
        for charge in charges:
            leg = prepareLegendV2(textSize=0.04,xmin=0.12,xmax=0.94)
            if charge=='allcharges': 
                sign=''
            else: 
                sign='+' if charge is 'plus' else '-'
            # leave a space before - sign, otherwise - is too close to W in the pdf (in the png it gets too far instead, ROOT magic!)
            thischannel = "W^{{{chsign}}} #rightarrow {fl}#nu".format(chsign=sign if sign != "-" else (" "+sign),fl=flavour)
            varRangeText = etaRangeText if any(options.target == x for x in ["ptxsec", "ptxsecnorm", "ptasym"]) else ptRangeText  
            header = "#bf{{Uncertainties on {p} for {ch}     {ptt}}}".format(p=poiName_target[options.target], ch=thischannel, ptt=varRangeText)
            leg.SetHeader(header)            
            leg.SetNColumns(4)
            lat = ROOT.TLatex()
            lat.SetNDC(); 
            lat.SetTextFont(42)
            quadrsum = summaries[(charge,groups[0])].Clone('quadrsum_{ch}'.format(ch=charge))
            totalerr = quadrsum.Clone('totalerr_{ch}'.format(ch=charge))
            for ing,ng in enumerate(groups):
                drawopt = 'pl' if ing==0 else 'pl same'
                summaries[(charge,ng)].Draw(drawopt)                
                if ing == 0: 
                    if options.target in ["ptxsec", "ptxsecnorm", "ptasym"]:
                        summaries[(charge,ng)].GetXaxis().SetRangeUser(ptmin,ptmax)
                    else:
                        summaries[(charge,ng)].GetXaxis().SetRangeUser(0.,2.4)
                if   ng=='binByBinStat': 
                    summaries[(charge,ng)].SetMarkerStyle(ROOT.kFullCircle)
                elif ng=='stat':
                    #print "bin  stat.uncertainty"
                    #for ibin in range(1,summaries[(charge,ng)].GetNbinsX()+1):
                    #    print "%s   %s" % (str(ibin-1),str(summaries[(charge,ng)].GetBinContent(ibin)))
                    summaries[(charge,ng)].SetMarkerStyle(ROOT.kFullCircle); 
                    summaries[(charge,ng)].SetMarkerColor(ROOT.kBlack); 
                    summaries[(charge,ng)].SetLineColor(ROOT.kBlack);
                elif ng=='luminosity'  : 
                    summaries[(charge,ng)].SetMarkerStyle(ROOT.kFullSquare)
                elif ng=='QEDTheo'  : 
                    summaries[(charge,ng)].SetMarkerStyle(ROOT.kOpenSquareDiagonal)
                elif ng=='L1Prefire'  : 
                    summaries[(charge,ng)].SetMarkerStyle(ROOT.kFourSquaresPlus)
                else: 
                    summaries[(charge,ng)].SetMarkerStyle(ROOT.kFullTriangleUp+ing)
                leg.AddEntry(summaries[(charge,ng)], niceSystName(ng), 'pl')
                # now compute the quadrature sum of all the uncertainties (neglecting correlations among nuis groups)
                for y in xrange(quadrsum.GetNbinsX()):
                    if options.target in ["ptxsec", "ptxsecnorm", "ptasym"]:
                        if options.ptMinSignal > 0 and genBins.ptBins[int(y)] < options.ptMinSignal:
                            quadrsum.SetBinContent(y+1,0.0)
                        else:
                            quadrsum.SetBinContent(y+1,math.hypot(quadrsum.GetBinContent(y+1),summaries[(charge,ng)].GetBinContent(y+1)))
                    else:
                        quadrsum.SetBinContent(y+1,math.hypot(quadrsum.GetBinContent(y+1),summaries[(charge,ng)].GetBinContent(y+1)))
            quadrsum.SetMarkerStyle(ROOT.kFullCrossX); 
            quadrsum.SetMarkerSize(3); 
            quadrsum.SetMarkerColor(ROOT.kBlack); 
            quadrsum.SetLineColor(ROOT.kBlack); 
            quadrsum.SetLineStyle(ROOT.kDashed)
            # quadrsum.Draw('pl same') 
            # now fill the real total error from the fit (with correct correlations)
            for y in xrange(totalerr.GetNbinsX()):
                if options.target in ["ptxsec", "ptxsecnorm", "ptasym"]:
                    if options.ptMinSignal > 0 and genBins.ptBins[int(y)] < options.ptMinSignal:
                        totalerr.SetBinContent(y+1,0.0)
                    else:
                        totalerr.SetBinContent(y+1,fitErrors['W{sign} {bin}'.format(sign=sign,bin=y)])
                else:
                    totalerr.SetBinContent(y+1,fitErrors['W{sign} {bin}'.format(sign=sign,bin=y)])
            totalerr.SetMarkerStyle(ROOT.kFullDoubleDiamond); 
            totalerr.SetMarkerSize(3); 
            totalerr.SetMarkerColor(ROOT.kRed+1); 
            totalerr.SetLineColor(ROOT.kRed+1);
            totalerr.Draw('pl same')
            #leg.AddEntry(quadrsum, 'Quadr. sum. impacts', 'pl')
            leg.AddEntry(totalerr, 'Total uncertainty', 'pl')
            leg.Draw('same')
            lat.DrawLatex(0.1, 0.92, '#bf{CMS}' + ('' if options.skipPreliminary else ' #it{Preliminary}'))
            lat.DrawLatex(0.78, 0.92, '35.9 fb^{-1} (13 TeV)')
            for i in ['pdf', 'png']:
                suff = '' if not options.suffix else '_'+options.suffix
                cs.SetLogy()
                cs.SaveAs(options.outdir+'/etaImpacts{rel}{suff}_{target}_{ptbin}{ch}.{i}'.format(rel='Abs' if options.absolute else 'Rel',
                                                                                                  suff=suff,target=options.target,i=i,
                                                                                                  ptbin=("ipt"+str(theptbin)+"_") if isSinglePtStrip else "",
                                                                                                  ch=charge))

