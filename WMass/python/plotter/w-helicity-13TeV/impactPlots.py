# USAGE: python impactPlots.py fitresults_poim1_exp0_bbb0.root -o plots/fit/absY/results/2018-12-10 --nuis 'CMS.*' --pois 'Wplus.*left.*'
import ROOT, os, datetime, re, operator, math
from array import array
from subMatrix import niceName
ROOT.gROOT.SetBatch(True)
import utilities
utilities = utilities.util()

def latexLabel(label):
    bin = int(label.split(' ')[-1]) if any(pol in label for pol in ['left','right','long']) else -1
    deltaYW = 0.2 ### should use binningYW.txt, but let's  not make it too complicated
    minY,maxY=bin*deltaYW,(bin+1)*deltaYW
    binLbl = ' {minY}<|Y_\mathrm{{ W }}|<{maxY}'.format(minY=minY,maxY=maxY)
    lblNoBin = ' '.join(label.split(' ')[:-1])
    lblNoBin = lblNoBin.replace('+','^+').replace(' left','_L').replace(' right','_R')
    lblNoBin += binLbl
    return lblNoBin

if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat('.3f')

    from optparse import OptionParser
    parser = OptionParser(usage='%prog fitresults.root [options] ')
    parser.add_option('-o','--outdir',     dest='outdir',     default='',   type='string', help='outdput directory to save the matrix')
    parser.add_option(     '--pois',       dest='pois',       default='',   type='string', help='pois for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option(     '--nuis',       dest='nuis',       default='',   type='string', help='nuis for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option(     '--nuisgroups', dest='nuisgroups', default='',   type='string', help='nuis groups for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option(     '--suffix',     dest='suffix',     default='',   type='string', help='suffix for the correlation matrix')
    parser.add_option(     '--target',     dest='target',     default='mu', type='string', help='target POI (can be mu,xsec,xsecnorm)')
    parser.add_option('-a','--absolute',   dest='absolute',   default=False, action='store_true', help='absolute uncertainty (default is relative)')
    parser.add_option(     '--latex',      dest='latex',      default=False, action='store_true', help='target POI (can be mu,xsec,xsecnorm)')
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
    valuesAndErrors = utilities.getFromHessian(args[0])

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
    if len(options.nuisgroups)==0:
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('pdf','')) if 'pdf' in x else 0)
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('muRmuF','')) if 'muRmuF' in x else 0)
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('muR','')) if ''.join([j for j in x if not j.isdigit()]) == 'muR'  else 0)
        nuisances = sorted(nuisances, key= lambda x: int(x.replace('muF','')) if ''.join([j for j in x if not j.isdigit()]) == 'muF'  else 0)
        nuisances = sorted(nuisances, key= lambda x: int(x.split('EffStat')[1]) if 'EffStat' in x else 0)

    #print "sorted pois = ", pois
    #print "\n\nsorted nuisances = ", nuisances

    c = ROOT.TCanvas("c","",1200,1200)
    c.SetGridx()
    c.SetGridy()

    c.SetLeftMargin(0.20)
    c.SetRightMargin(0.20)
    c.SetBottomMargin(0.20)

    ## make the new, smaller TH2F correlation matrix
    nbinsx = len(pois); nbinsy = len(nuisances)
    th2_sub = ROOT.TH2F('sub_imp_matrix', '', nbinsx, 0., nbinsx, nbinsy, 0., nbinsy)
    th2_sub.GetXaxis().SetTickLength(0.)
    th2_sub.GetYaxis().SetTickLength(0.)
    th2_sub.GetZaxis().SetTitle("impact on POI {units}".format(units='' if options.absolute else '(%)'))

    ## pretty nested loop. enumerate the tuples
    for i,x in enumerate(pois):
        for j,y in enumerate(nuisances):
            ## set it into the new sub-matrix
            if options.absolute: 
                val = mat[(x,y)]
            else: 
                val = 100*mat[(x,y)]/valuesAndErrors[x+'_'+target][0] if valuesAndErrors[x+'_'+target][0] !=0 else 0.0
            th2_sub.SetBinContent(i+1, j+1, val)
            ## set the labels correctly
            new_x = niceName(x)
            new_y = niceName(y)
            th2_sub.GetXaxis().SetBinLabel(i+1, new_x)
            th2_sub.GetYaxis().SetBinLabel(j+1, new_y)
            

    rmax = max(abs(th2_sub.GetMaximum()),abs(th2_sub.GetMinimum()))
    th2_sub.GetZaxis().SetRangeUser(-rmax,rmax)
    th2_sub.GetXaxis().LabelsOption("v")
    if options.absolute: 
        ROOT.gStyle.SetPaintTextFormat('1.3f')
    else: 
        ROOT.gStyle.SetPaintTextFormat('1.1f')
    if len(pois)<30 and len(nuisances)<30: th2_sub.Draw('colz0 text45')
    else: th2_sub.Draw('colz0')
    if len(nuisances)>30: th2_sub.GetYaxis().SetLabelSize(0.02)


    nuisName = 'NuisGroup'+options.nuisgroups if len(options.nuisgroups) else options.nuis
    nuisName = nuisName.replace(',','AND').replace('.','').replace('*','').replace('$','').replace('^','').replace('|','').replace('[','').replace(']','')
    poisName = options.pois.replace(',','AND').replace('.','').replace('*','').replace('$','').replace('^','').replace('|','').replace('[','').replace(']','')

    suff = '' if not options.suffix else '_'+options.suffix
    poisNameNice = poisName.replace('(','').replace(')','').replace('\|','or')
    if options.latex:
        txtfilename = 'smallImpacts{rel}{suff}_{target}_{nn}_On_{pn}.tex'.format(rel='Abs' if options.absolute else 'Rel',suff=suff,target=target,i=i,nn=nuisName,pn=poisNameNice)
        if options.outdir: txtfilename = options.outdir + '/' + txtfilename
        txtfile = open(txtfilename,'w')
        txtfile.write("\\begin{tabular}{l "+"   ".join(['r' for i in xrange(th2_sub.GetNbinsX())])+"} \\hline \n")
        txtfile.write('                    ' + " & ".join(['{poi}'.format(poi=latexLabel(th2_sub.GetXaxis().GetBinLabel(i+1))) for i in xrange(th2_sub.GetNbinsX())]) + "\\\\ \\hline \n")
        for j in xrange(th2_sub.GetNbinsY()):
            txtfile.write('{label:<20}& '.format(label=th2_sub.GetYaxis().GetBinLabel(j+1)) + " & ".join(['{syst:^.2f}'.format(syst=th2_sub.GetBinContent(i+1,j+1)) for i in xrange(th2_sub.GetNbinsX())]) + "\\\\ \n")
        txtfile.write("\\end{tabular}\n")
        txtfile.close()
        print "Saved the latex table in file ",txtfilename

    if options.outdir and not options.latex:
        for i in ['pdf', 'png']:
            c.SaveAs(options.outdir+'/smallImpacts{rel}{suff}_{target}_{nn}_On_{pn}.{i}'.format(rel='Abs' if options.absolute else 'Rel',suff=suff,target=target,i=i,nn=nuisName,pn=poisNameNice))
        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=options.outdir))
