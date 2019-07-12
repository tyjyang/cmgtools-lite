import ROOT, itertools, datetime, math, copy, os, commands, sys
from array import array


## usage:  python templatesRoofit.py --infile /eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/ntuplesRecoil/TREES_SIGNAL_1l_recoil_fullTrees/friends/ --nlo --outdir ~/www/private/w-helicity-13TeV/cosThetaFits/2018-10-10-fineBin/

## usage add --batch to run on the batch

def fixFits():
    bf_txt = open(options.fixFits+'/badfits.txt', 'r')
    bf_ls  = bf_txt.readlines()
    #bf_ls = [i.replace('norm_','').strip() for i in bf_ls]
    bf_ls = [i.strip() for i in bf_ls]
    bf_ls = list(set(bf_ls))
    nbadfits = 0
    badfits = {}
    for bf in bf_ls:
        if '16_19' in bf or '16_0' in bf: continue

        print 'really found a bad fit', bf

        var = bf.split('_')[-1]
        varf = [i for i in os.listdir(options.fixFits) if var+'.root' in i][0]
        if not badfits.has_key(varf):
            badfits[varf] = []
            badfits[varf].append(bf) #.split('_')[:3])
        else:
            badfits[varf].append(bf) #.split('_')[:3])
            nbadfits += 1

    print 'found', nbadfits, 'bad fits'
    print badfits

    ## we also need the three analytic helicity fraction functions
    
    ana_r = ROOT.TF1('ana_r', '3./8.*(1+x)^2')
    ana_l = ROOT.TF1('ana_l', '3./8.*(1-x)^2')
    ana_0 = ROOT.TF1('ana_0', '3./4*(TMath::Sqrt(1-x*x))^2')

    fixedHists = []
    c = ROOT.TCanvas()
    lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextSize(0.03)

    for k,v in badfits.items():
        tmp_f = ROOT.TFile(options.fixFits+'/'+k,'update')
        if not len(tmp_f.GetListOfKeys()): 
            print 'this file is screwed up', tmp_f.GetName()
            continue

        for h in v:
            print 'from file', tmp_f.GetName(), 'reading histogram', h
            tmp_hist_orig = tmp_f.Get(h)
            tmp_hist = tmp_hist_orig.Clone(h+'_fixed')

            ## bit of a duplication here
            fitterR = ROOT.TF1('helR_'+h,'[2] * ( (1.-[1]-[0])*{wr} + [1]*{wl} + [0]*{w0} )'.format(wr=ana_r.GetExpFormula(), wl=ana_l.GetExpFormula(), w0=ana_0.GetExpFormula()), -1., 1.)
            fitterL = ROOT.TF1('helL_'+h,'[2] * ( (1.-[1]-[0])*{wr} + [1]*{wl} + [0]*{w0} )'.format(wr=ana_r.GetExpFormula(), wl=ana_l.GetExpFormula(), w0=ana_0.GetExpFormula()), -1., 1.)
            fitter0 = ROOT.TF1('hel0_'+h,'[2] * ( (1.-[1]-[0])*{wr} + [1]*{wl} + [0]*{w0} )'.format(wr=ana_r.GetExpFormula(), wl=ana_l.GetExpFormula(), w0=ana_0.GetExpFormula()), -1., 1.)
            fitterR .SetLineColor(ROOT.kBlue); fitterR .SetLineWidth(2)
            fitterL .SetLineColor(ROOT.kBlue); fitterL .SetLineWidth(2)
            fitter0 .SetLineColor(ROOT.kBlue); fitter0 .SetLineWidth(2)
    
            fitterR.SetParLimits(0, 0., 0.45); fitterR.SetParLimits(1, 0.25, 1.)
            fitterL.SetParLimits(0, 0., 0.45); fitterL.SetParLimits(1, 0.25, 1.)
            fitter0.SetParLimits(0, 0., 0.45); fitter0.SetParLimits(1, 0.25, 1.)

            tmp_hist.Fit(fitterR.GetName(), '', '', -1., 1.)
            tmp_hist.Fit(fitterL.GetName(), '', '', -1., 1.)
            tmp_hist.Fit(fitter0.GetName(), '', '', -1., 1.)

            f0 = fitterR.GetParameter(0); f0_err = fitterR.GetParError(0)
            fL = fitterR.GetParameter(1); fL_err = fitterR.GetParError(1)
            fR = 1.-fL-f0               ; fR_err = math.sqrt(f0_err**2 + fL_err**2)
    
            tmp_hist.Draw()

            lat.DrawLatex(0.45, 0.85, 'f0: {a1:.4f} +- {a2:.4f}'.format(a1=f0,a2=f0_err))
            lat.DrawLatex(0.45, 0.80, 'fL: {a1:.4f} +- {a2:.4f}'.format(a1=fL,a2=fL_err))
            lat.DrawLatex(0.45, 0.75, 'fR: {a1:.4f} +- {a2:.4f}'.format(a1=fR,a2=fR_err))
    
            c.SaveAs(options.fixFits+'/fixedFits/'+tmp_hist.GetName()+'.png')
            c.SaveAs(options.fixFits+'/fixedFits/'+tmp_hist.GetName()+'.pdf')

            ## now update the fractions in the 2d histograms
            ch = h.split('_')[3]
            vfull = h.split('_')[-1]
            vind = int(vfull[3:]) if 'pdf' in vfull else 0 if 'nominal' in vfull else 61 if 'alphaSUp' in vfull else 62
            tmp_fractionL = tmp_f.Get('fractionL_'+ch+'_'+str(vind))
            tmp_fractionR = tmp_f.Get('fractionR_'+ch+'_'+str(vind))
            tmp_fraction0 = tmp_f.Get('fraction0_'+ch+'_'+str(vind))

            binx, biny = int(h.split('_')[1])+1, int(h.split('_')[2])+1 ## +1 needed because i was stupid. who knows why. a great mystery of science

            f0_old = tmp_fraction0.GetBinContent(binx,biny)
            fR_old = tmp_fractionR.GetBinContent(binx,biny)
            fL_old = tmp_fractionL.GetBinContent(binx,biny)

            oldstr = 'f0: {f0:.4f}, fR: {fR:.4f}, fL: {fL:.4f}'.format(f0=f0_old, fR=fR_old, fL=fL_old)
            newstr = 'f0: {f0:.4f}, fR: {fR:.4f}, fL: {fL:.4f}'.format(f0=f0    , fR=fR    , fL=fL    )

            print '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
            print 'for binx {bx} and biny {by}, changing the fractions from '.format(bx=binx-1, by=biny-1)
            print 'OLD: {old} to'.format(old=oldstr)
            print 'NEW: {new}   '.format(new=newstr)
            print '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
            
            ## actually doing the change
            f0_old = tmp_fraction0.SetBinContent(binx,biny, f0); f0_old = tmp_fraction0.SetBinError(binx,biny, f0_err)
            fR_old = tmp_fractionR.SetBinContent(binx,biny, fR); fR_old = tmp_fractionR.SetBinError(binx,biny, fR_err)
            fL_old = tmp_fractionL.SetBinContent(binx,biny, fL); fL_old = tmp_fractionL.SetBinError(binx,biny, fL_err)
            
            tmp_fractionL.Write('', ROOT.TObject.kWriteDelete)
            tmp_fractionR.Write('', ROOT.TObject.kWriteDelete)
            tmp_fraction0.Write('', ROOT.TObject.kWriteDelete)

            fixedHists.append(copy.deepcopy(tmp_hist)) 

        tmp_f.Close()

    os.system('mkdir -p '+options.fixFits+'/fixedFits/')
    os.system('cp ~/public/index.php '+options.fixFits+'/fixedFits/')


    fixedFile = ROOT.TFile(options.fixFits+'/fixedFits/fixedFits.root', 'recreate')
    for f in fixedHists:
        f.Write()
    fixedFile.Close()
    

def makeCondorFile():
    condor_file = open(options.outdir+'/submitFile.condor','w')
    condor_file.write('''Universe = vanilla
Executable = {od}/dummy.sh
use_x509userproxy = $ENV(X509_USER_PROXY)
Log        = {od}/$(ProcId).log   
Output     = {od}/$(ProcId).out   
Error      = {od}/$(ProcId).error
getenv      = True
environment = "LS_SUBCWD={here}"
+MaxRuntime = 36000\n\n
'''.format(od=options.outdir,here=os.environ['PWD'] ) )
    if os.environ['USER'] in ['mdunser', 'psilva']:
        condor_file.write('+AccountingGroup = "group_u_CMST3.all"\n\n')

    for ipdf in range(63):
        condor_file.write('arguments = templatesRoofit.py --infile {inf} --outdir {od} --nlo --ipdf {p}\n'.format(inf=options.infile, od=options.outdir, p=ipdf))
        condor_file.write('queue 1\n\n')
    
    condor_file.close()

def formatHisto(hist):
    hist.GetXaxis().SetTitleOffset(1.02)
    hist.GetXaxis().SetTitleSize(0.06)
    hist.GetXaxis().SetLabelSize(0.06)

    hist.GetYaxis().SetTitleOffset(1.02)
    hist.GetYaxis().SetTitleSize(0.06)
    hist.GetYaxis().SetLabelSize(0.06)

    hist.GetZaxis().SetLabelSize(0.06)

def symmetrizeFractions(fractionR, fractionL, fraction0, nsmooth=0):
    ## symmetrize the fractions
    newcontents = {}
    nbinsY = fractionR.GetXaxis().GetNbins()
    nbinsZ = fractionR.GetYaxis().GetNbins()

    fractionR_sym = fractionR.Clone(fractionR.GetName()+'_sym')
    fractionL_sym = fractionL.Clone(fractionL.GetName()+'_sym')
    fraction0_sym = fraction0.Clone(fraction0.GetName()+'_sym')
    fractionR_sym.Reset()
    fractionL_sym.Reset()
    fraction0_sym.Reset()

    ## careful here, this symmetrization makes f0 flat in Y

    for iy in range(1,nbinsY+1):
        for iz in range(1,nbinsZ+1):
            newcontents[('R', iy, iz)] = 1./2.*sum( fractionR.GetBinContent(i, iz) for i in [iy, nbinsY+1-iy] )
            newcontents[('L', iy, iz)] = 1./2.*sum( fractionL.GetBinContent(i, iz) for i in [iy, nbinsY+1-iy] )
            newcontents[('0', iy, iz)] = 1./nbinsY*sum( fraction0.GetBinContent(i, iz) for i in range(1,nbinsZ+1) )
    for iy in range(1,nbinsY+1):
        for iz in range(1,nbinsZ+1):
            fractionR_sym.SetBinContent(iy, iz, newcontents[('R', iy, iz)])
            fractionL_sym.SetBinContent(iy, iz, newcontents[('L', iy, iz)])
            fraction0_sym.SetBinContent(iy, iz, newcontents[('0', iy, iz)])
    
    for i in range(nsmooth):
        fractionR.Smooth()
        fractionL.Smooth()
        fraction0.Smooth()

    return [fractionR_sym, fractionL_sym, fraction0_sym]

def saveCanv(c, name):
    c.SaveAs(name+'.pdf')
    c.SaveAs(name+'.png')

from optparse import OptionParser
parser = OptionParser(usage='python %prog --infile directoryWithRootFile --outdir outputPlotDirectory ( --preFSR / --lorenzo )')
parser.add_option('--roofit' , dest='doRooFit' , default=False, action='store_true', help='Use fancy RooFit generic pdf fitting instead of simple chi2.')
parser.add_option('--generic', dest='doGeneric', default=True , action='store_true', help='Use generic pdf instead of whatever else')
parser.add_option('--nlo'    , dest='doNLO'    , default=False, action='store_true', help='Use the amc@nlo sample.')
parser.add_option('--bins-atlas', dest='atlasBins'    , default=False, action='store_true', help='Use the binning of the ATLAS paper.')
parser.add_option('-i','--infile', dest='infile', default='', type='string', help='Specify the input file. should be a single (big) friendtree.')
parser.add_option('-o','--outdir', dest='outdir', default='', type='string', help='Specify the output directory for the plots. It makes a lot of plots.')
parser.add_option('-l','--lorenzo', dest='doLorenzo', default=False, action='store_true', help='use lorenzo\'s trees')
parser.add_option('--preFSR', default=True, action='store_true', help='use preFSR leptons/Ws')
parser.add_option('--batch', default=False, action='store_true', help='run on the batch')
parser.add_option('--ipdf', default=-1, type='int', help='run only this pdf')
parser.add_option('--fixFits', default='', type='string', help='fix cosine theta fits that failed from input directory')
(options, args) = parser.parse_args()


ROOT.gROOT.SetBatch()

if options.fixFits:
    if not os.path.isdir(options.fixFits):
        raise RuntimError, 'you should give the input directory in which you want to fix the fits as argument of the option!'
    else:
        fixFits()
        sys.exit(0)


if not (options.outdir or options.infile):
    raise RuntimeError, 'You have to give an input file and an output directory!'

plotsdir = '{od}/chi2_muEl/'.format(od=options.outdir)

if not os.path.isdir(plotsdir):
    os.system('mkdir -p {pd}'.format(pd=plotsdir))
    os.system('cp ~/index.php {pd}/'.format(pd=plotsdir))

if options.batch:
    makeCondorFile()
    os.system('cp dummy.sh {pd}/../'.format(pd=plotsdir))
    print 'condor_submit {od}/submitFile.condor'.format(od=options.outdir)
    sys.exit(0)

maxYW = 10.

## some tree specific things. we need the name of the 
## variables of the W in the tree

if options.doLorenzo:
    var_wy  = 'WpreFSR_y'
    var_wpt = 'WpreFSR_qt'
    var_wch = 'mu_charge'
    var_dec = 'WpreFSR_decayId'
    var_cos = 'WpreFSR_cosCS'
    var_phi = 'WpreFSR_phiCS'
    
    ## weightstring
    gen_weight = 'weights[0]/abs(weights[0])'

else:
    if options.preFSR:
        var_wy  = 'abs(prefsrw_y)'
        var_wpt = 'prefsrw_pt'
        var_wch = 'prefsrw_charge'
        var_dec = 'prefsrw_decayId'
        var_cos = 'prefsrw_costcs'

    else:
        var_wy  = 'abs(genw_y)'
        var_wpt = 'genw_pt'
        var_wch = 'genw_charge'
        var_dec = 'genw_decayId'
        var_cos = 'genw_costcs'
    
    ## weightstring
    gen_weight = 'weightGen/abs(weightGen)'


## get the tree from the file

if os.path.isfile(options.infile):
    infile = ROOT.TFile(options.infile, 'read')
    tree = infile.Get('Friends' if not options.doLorenzo else 'tree')
else:
    tree = ROOT.TChain('Friends' if not options.doLorenzo else 'tree')
    ## rootlist = list( i for i in os.listdir(options.infile) if '.root' in i)
    #rootlist = list( i for i in commands.getoutput('eos ls '+options.infile).split('\n') if '.root' in i)
    rootlist = [i for i in os.listdir(options.infile) if '.root' in i and 'WJetsToLNu' in i]
    print 'adding these files', rootlist
    for f in rootlist:
        tree.Add('root://eoscms.cern.ch//'+options.infile+'/'+f)


## first build the quantiles in both rapidity and pT


if not options.atlasBins:
    ## make the quantiles for YW
    ## ===================================
    ## fixed fine binning print('FILLING FOR THE W-rapidity QUANTILES')
    ## fixed fine binning 
    ## fixed fine binning nqy = 15
    ## fixed fine binning h_tmp_wy  = ROOT.TH1F('h_tmp_wy','quantile calculation Yw', 1000, 0., maxYW)
    ## fixed fine binning tree.Draw(var_wy+'>>h_tmp_wy', '('+var_wch + ' > 0 && abs(genw_decayId)==14 )*'+gen_weight)
    ## fixed fine binning xqy = array('d', [i/float(nqy)+1./float(nqy) for i in range(nqy)])
    ## fixed fine binning yqy = array('d', [1. for i in range(nqy)])
    ## fixed fine binning foobar = h_tmp_wy.GetQuantiles(nqy,yqy,xqy);
    ## fixed fine binning 
    ## fixed fine binning yqnewy = array('d', [0.] + [float('{a:.3f}'.format(a=i)) for i in yqy])
    
    yqnewy = [0.0+i*0.25 for i in range(17)] + [maxYW] ## this is absy
    
    wrap_nbins = len(yqnewy)-1
    wrap_bins = array('d', yqnewy)
    
    ## make the quantiles for W-pT
    ## ===================================
    
    print('FILLING FOR THE W-pT QUANTILES')
    ## fixed binning h_tmp_wpt = ROOT.TH1F('h_tmp_wpt','quantile calculation', 1000, 0., 100.)
    ## fixed binning tree.Draw(var_wpt+'>>h_tmp_wpt', '(abs(genw_decayId)==14 || abs(genw_decayId)== 12) *'+gen_weight)
    ## fixed binning xq = array('d', [i/float(nq)+1./float(nq) for i in range(int(nq))])
    ## fixed binning yq = array('d', [1 for i in range(nq)])
    ## fixed binning h_tmp_wpt.GetQuantiles(nq,yq,xq);
    ## fixed binning 
    ## fixed binning 
    ## fixed binning yqnew = array('d', [0.] + [float('{a:.3f}'.format(a=i)) for i in yq])
    ## fixed binning hnew = ROOT.TH1F('h_wpt_quantiles', 'wpt rebinned', len(yq), yqnew)
    nq = 20
    yqnew = array('d', [0.0, 1.971, 2.949, 3.838, 4.733, 5.674, 6.684, 7.781, 8.979, 10.303, 11.777, 13.435, 15.332, 17.525, 20.115, 23.245, 27.173, 32.414, 40.151, 53.858, 100.0])
    #nq =  10 ## number of quantiles
    #yqnew = array('d', [0.0, 2.949, 4.733, 6.684, 8.979, 11.777, 15.332, 20.115, 27.173, 40.151, 100.0])
    
    wpt_nbins = nq
    wpt_bins = yqnew

else:
    wrap_bins   = array('d', [0.20*i for i in range(15)] + [3., 3.25, 3.5, 3.75, 4., 4.25, 4.5, 5., maxYW])
    wrap_nbins  = len(wrap_bins)-1
    wpt_bins  = array('d', range(11) + [12., 14., 16., 18., 20., 25., 30., 35., 40., 50., 60., 70., 100.])
    wpt_nbins = len(wpt_bins)-1

print('this is the YW binning', wrap_bins)
print('these are the w-pt quantiles: ', wpt_bins)

## ===================================
## finished with the quantile production. we now have wpt_bins and wrap_bins



## we also need the three analytic helicity fraction functions

ana_r = ROOT.TF1('ana_r', '3./8.*(1+x)^2')
ana_l = ROOT.TF1('ana_l', '3./8.*(1-x)^2')
ana_0 = ROOT.TF1('ana_0', '3./4*(TMath::Sqrt(1-x*x))^2')

## done making the functions

lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextSize(0.03)


## build normalized templates of CM for each of the helicities
cost_nbins =  50
cost_bins = array('d', [2*i/float(cost_nbins)-1 for i in range(cost_nbins+1)])
    
nfills = int(1e7)
h_ana_wr = ROOT.TH1F('h_analytic_wr', 'h_analytic_wr', cost_nbins, -1., 1.); h_ana_wr.FillRandom(ana_r.GetName(), nfills); h_ana_wr.Scale(1./nfills);
h_ana_wl = ROOT.TH1F('h_analytic_wl', 'h_analytic_wl', cost_nbins, -1., 1.); h_ana_wl.FillRandom(ana_l.GetName(), nfills); h_ana_wl.Scale(1./nfills);
h_ana_w0 = ROOT.TH1F('h_analytic_w0', 'h_analytic_w0', cost_nbins, -1., 1.); h_ana_w0.FillRandom(ana_0.GetName(), nfills); h_ana_w0.Scale(1./nfills);
h_ana_wr .GetXaxis().SetTitle('cos #theta')
h_ana_wl .GetXaxis().SetTitle('cos #theta')
h_ana_w0 .GetXaxis().SetTitle('cos #theta')

## write the histograms into a file
date = datetime.date.today().isoformat()
ofn = 'fractions_roofit_histos_chi2_{date}_muEl_plusMinus{lor}{fsr}'.format(date=date,lor = '_lorenzo' if options.doLorenzo else '',fsr='_preFSR' if options.preFSR else 'dressed')
if options.ipdf > -1:
    ofn = options.outdir+'/'+ofn
    ofn+= '_pdf'+str(options.ipdf) if options.ipdf  and options.ipdf < 61 else '_nominal' if not options.ipdf else '_alphaSUp' if options.ipdf == 61 else '_alphaSDown'
ofn +='.root'
outfile = ROOT.TFile(ofn, 'recreate')
outfile.cd()
    

for pdfWgt in range(63):

    if options.ipdf > -1 and options.ipdf != pdfWgt:
        continue

    pdfWgtString = 'nominal' if not pdfWgt else 'pdf'+str(pdfWgt) if pdfWgt < 61 else 'alphaSUp' if pdfWgt == 61 else 'alphaSDown'

    print 'AT PDF WEIGHT {w}'.format(w=pdfWgtString)

    ## first we need a bunch of histograms in which we store the weights
    ## as a function of cos(theta)* and the helicity fractions :
    
    h3_rapVsCMvsWPtPlus  = ROOT.TH3D('h3_rapVsCMvsWPtPlus' +str(pdfWgt),'h3_rapVsCMvsWPtPlus' , cost_nbins, cost_bins, wrap_nbins, wrap_bins, wpt_nbins, wpt_bins); h3_rapVsCMvsWPtPlus .Sumw2()
    h3_rapVsCMvsWPtMinus = ROOT.TH3D('h3_rapVsCMvsWPtMinus'+str(pdfWgt),'h3_rapVsCMvsWPtMinus', cost_nbins, cost_bins, wrap_nbins, wrap_bins, wpt_nbins, wpt_bins); h3_rapVsCMvsWPtMinus.Sumw2()
    
    h3_rapVsCMvsWPtPlus .GetXaxis().SetTitle('cos #theta'); h3_rapVsCMvsWPtPlus .GetYaxis().SetTitle('Y_{W}'); h3_rapVsCMvsWPtPlus .GetZaxis().SetTitle('W-p_{T}')
    h3_rapVsCMvsWPtMinus.GetXaxis().SetTitle('cos #theta'); h3_rapVsCMvsWPtMinus.GetYaxis().SetTitle('Y_{W}'); h3_rapVsCMvsWPtMinus.GetZaxis().SetTitle('W-p_{T}')
    
    nbinsX = h3_rapVsCMvsWPtPlus.GetXaxis().GetNbins()
    nbinsY = h3_rapVsCMvsWPtPlus.GetYaxis().GetNbins()
    nbinsZ = h3_rapVsCMvsWPtPlus.GetZaxis().GetNbins()
    
    typeString = ' preFSR leptons' if options.preFSR else 'dressed leptons'
    typeString = typeString + pdfWgtString
    
    ## project pT versus rapidity for the fractions
    fractionR_plus = h3_rapVsCMvsWPtPlus.Project3D('zy').Clone('fractionR_plus_'+str(pdfWgt)); fractionR_plus.SetTitle('W^{+}: fractions R: '+typeString); fractionR_plus.Reset()
    fractionL_plus = h3_rapVsCMvsWPtPlus.Project3D('zy').Clone('fractionL_plus_'+str(pdfWgt)); fractionL_plus.SetTitle('W^{+}: fractions L: '+typeString); fractionL_plus.Reset()
    fraction0_plus = h3_rapVsCMvsWPtPlus.Project3D('zy').Clone('fraction0_plus_'+str(pdfWgt)); fraction0_plus.SetTitle('W^{+}: fractions 0: '+typeString); fraction0_plus.Reset()
    formatHisto(fractionR_plus)
    formatHisto(fractionL_plus)
    formatHisto(fraction0_plus)
    
    fractionR_minus = h3_rapVsCMvsWPtMinus.Project3D('zy').Clone('fractionR_minus_'+str(pdfWgt)); fractionR_minus.SetTitle('W^{-}:fractions R: '+typeString); fractionR_minus.Reset()
    fractionL_minus = h3_rapVsCMvsWPtMinus.Project3D('zy').Clone('fractionL_minus_'+str(pdfWgt)); fractionL_minus.SetTitle('W^{-}:fractions L: '+typeString); fractionL_minus.Reset()
    fraction0_minus = h3_rapVsCMvsWPtMinus.Project3D('zy').Clone('fraction0_minus_'+str(pdfWgt)); fraction0_minus.SetTitle('W^{-}:fractions 0: '+typeString); fraction0_minus.Reset()
    formatHisto(fractionR_minus)
    formatHisto(fractionL_minus)
    formatHisto(fraction0_minus)
    
    ## done making the histograms
    
    genWeightPDF = '({pdfw}*{gw})'.format(pdfw='1.' if not pdfWgt else 'hessWgt'+str(pdfWgt) if pdfWgt < 61 else 'qcd_alphaSUp' if pdfWgt == 61 else 'qcd_alphaSDn', gw=gen_weight)

    print('FILLING THE 3D HISTOGRAM!! with weight', genWeightPDF)
    genIdCut = '(genw_decayId == 14 || genw_decayId == 12)' if not options.doLorenzo else ' 1.'
    chargePos = var_wch + (' > ' if not options.doLorenzo else ' < ')+ '0' 
    chargeNeg = var_wch + (' < ' if not options.doLorenzo else ' > ')+ '0' 
    tree.Draw(var_wpt+':'+var_wy+':'+var_cos+'>>h3_rapVsCMvsWPtPlus' +str(pdfWgt), '('+chargePos + ' && '+genIdCut+' )* '+genWeightPDF, '')
    tree.Draw(var_wpt+':'+var_wy+':'+var_cos+'>>h3_rapVsCMvsWPtMinus'+str(pdfWgt), '('+chargeNeg + ' && '+genIdCut+' )* '+genWeightPDF, '')
    print('... done filling the 3D histogram. writing it for future use.')
    h3_rapVsCMvsWPtPlus  .Write()
    h3_rapVsCMvsWPtMinus .Write()
    
    ## END FILLING WPT YW CM HISTOGRAM
    ## ======================================
    
    ## f_wr = ROOT.TF1('f_analytic_wr', '3./8.*(1+x)^2'              , -1., 1.)
    ## f_wl = ROOT.TF1('f_analytic_wl', '3./8.*(1-x)^2'              , -1., 1.)
    ## f_w0 = ROOT.TF1('f_analytic_w0', '3./4*(TMath::Sqrt(1-x*x))^2', -1., 1.)
    
    c = ROOT.TCanvas('canv'+str(pdfWgt), '')
    nfits = 0
    badfits = {}
    allfits = {}
    
    ROOT.gStyle.SetOptStat(0)
    
    for iy in range(1,nbinsY+1):
        for iz in range(1,nbinsZ+1):
            for wch in ['plus', 'minus']:
                nfits += 1
                #if nfits > 2: continue
                c.Clear()
                pos = wch == 'plus'
                name = '{iy}_{iz}_{ch}'.format(iy=iy-1,iz=iz-1,ch=wch)
                print('at bin {n}'.format(n=name))
                c.SetName ('canv_'+name+'_'+pdfWgtString)
                c.SetTitle('canv_'+name+pdfWgtString)
                h3_rapVsCMvsWPt = h3_rapVsCMvsWPtPlus if wch == 'plus' else h3_rapVsCMvsWPtMinus
                h_cm = h3_rapVsCMvsWPt.ProjectionX(name, iy if iy else 1, iy if iy else nbinsY, iz if iz else 1, iz if iz else nbinsZ)
                #h_cm.SetTitle('cosTheta_'+name)
                ylo = h3_rapVsCMvsWPt.GetYaxis().GetBinLowEdge(iy)
                yhi = h3_rapVsCMvsWPt.GetYaxis().GetBinUpEdge(iy)
                plo = h3_rapVsCMvsWPt.GetZaxis().GetBinLowEdge(iz)
                phi = h3_rapVsCMvsWPt.GetZaxis().GetBinUpEdge(iz)
                h_cm.SetTitle('{ch}: y_{{W}} #in [{ylo:.2f},{yhi:.2f}] , p_{{T}}^{{W}} #in [{plo:.2f},{phi:.2f}]'.format(ch='W^{+}' if wch == 'plus' else 'W^{-}', ylo=ylo,yhi=yhi,plo=plo,phi=phi)+typeString )
    
                fitterR = ROOT.TF1('helR_'+name,'[2] * ( (1.-[1]-[0])*{wr} + [1]*{wl} + [0]*{w0} )'.format(wr=ana_r.GetExpFormula(), wl=ana_l.GetExpFormula(), w0=ana_0.GetExpFormula()), -1., 1.)
                fitterL = ROOT.TF1('helL_'+name,'[2] * ( (1.-[1]-[0])*{wr} + [1]*{wl} + [0]*{w0} )'.format(wr=ana_r.GetExpFormula(), wl=ana_l.GetExpFormula(), w0=ana_0.GetExpFormula()), -1., 1.)
                fitter0 = ROOT.TF1('hel0_'+name,'[2] * ( (1.-[1]-[0])*{wr} + [1]*{wl} + [0]*{w0} )'.format(wr=ana_r.GetExpFormula(), wl=ana_l.GetExpFormula(), w0=ana_0.GetExpFormula()), -1., 1.)
    
                fitterR.SetParLimits(0, 0., 0.55); fitterR.SetParLimits(1, 0.25, 1.)
                fitterL.SetParLimits(0, 0., 0.55); fitterL.SetParLimits(1, 0.25, 1.)
                fitter0.SetParLimits(0, 0., 0.55); fitter0.SetParLimits(1, 0.25, 1.)
    
                h_cm_norm = h_cm.Clone('norm_'+h_cm.GetName()+'_'+pdfWgtString)
                h_cm_norm.Scale(1./h_cm_norm.Integral())
                h_cm_norm.Draw()
    
                h_cm_norm.Fit(fitterR.GetName(), '', '', -1., 1.)
                h_cm_norm.Fit(fitterL.GetName(), '', '', -1., 1.)
                h_cm_norm.Fit(fitter0.GetName(), '', '', -1., 1.)

                h_cm_norm.Write()
    
                f0 = fitterR.GetParameter(0)
                fL = fitterR.GetParameter(1)
                fR = 1.-fL-f0 #fitterR.GetParameter(0)
    
                # if f0 < 0.01: f0 = 0.01
    
                f0_err = fitterR.GetParError(0)
                fL_err = fitterR.GetParError(1)
                fR_err = math.sqrt(f0_err**2 + fL_err**2)
    
                print('fractions: \t fL: {fL:.4f} +- {fLe:.4f} \t fR: {fR:.4f} +- {fRe:.4f} \t f0: {f0:.4f} +- {f0e:.4f}'.format(fL=fL, fLe=fL_err, fR=fR, fRe=fR_err, f0=f0, f0e=f0_err))
                (fractionR_plus if pos else fractionR_minus).SetBinContent(iy, iz, fR); (fractionR_plus if pos else fractionR_minus).SetBinError(iy, iz, fR_err)
                (fractionL_plus if pos else fractionL_minus).SetBinContent(iy, iz, fL); (fractionL_plus if pos else fractionL_minus).SetBinError(iy, iz, fL_err)
                (fraction0_plus if pos else fraction0_minus).SetBinContent(iy, iz, f0); (fraction0_plus if pos else fraction0_minus).SetBinError(iy, iz, f0_err)
    
                lat.DrawLatex(0.45, 0.85, 'f0: {a1:.4f} +- {a2:.4f}'.format(a1=f0,a2=f0_err))
                lat.DrawLatex(0.45, 0.80, 'fL: {a1:.4f} +- {a2:.4f}'.format(a1=fL,a2=fL_err))
                lat.DrawLatex(0.45, 0.75, 'fR: {a1:.4f} +- {a2:.4f}'.format(a1=fR,a2=fR_err))
    
                chi2 = fitterR.GetChisquare(); ndf = fitterR.GetNDF(); 
                print('nev in the fitting: {n}'.format(n=h_cm.Integral()))
                if h_cm.Integral() < 1500.:
                    print('VERY FEW EVENTS!!!!')
                    print(name)
                    print('nev = {v}'.format(v=h_cm.Integral()))
                if chi2/ndf > 5.0:
                    print('VERY BAD FIT!!!! in bin', name, 'for pdf weight', pdfWgtString, '. chi2/ndf = {f}'.format(f=chi2/ndf))
                    badfits[name] = chi2/ndf
                    os.system('echo '+h_cm_norm.GetName()+' >> '+options.outdir+'/badfits.txt')
                    ## save the plots only for the bad fits!
                    lat.DrawLatex(0.45, 0.68, '#chi^{{2}}/ndf: {a1:.2f}/{a2} = {a3:.2f}'.format(a1=chi2, a2=ndf, a3=chi2/ndf))
                    c.SaveAs('{pd}/cosTheta_chi2_{n}_{pdf}.pdf'.format(pd=plotsdir,n=name, pdf= pdfWgtString))
                    c.SaveAs('{pd}/cosTheta_chi2_{n}_{pdf}.png'.format(pd=plotsdir,n=name, pdf= pdfWgtString))
                allfits[name] = chi2/ndf
    
                ## lat.DrawLatex(0.45, 0.68, '#chi^{{2}}/ndf: {a1:.2f}/{a2} = {a3:.2f}'.format(a1=chi2, a2=ndf, a3=chi2/ndf))
                # don't make all the plots... c.SaveAs('{pd}/cosTheta_chi2_{n}.pdf'.format(pd=plotsdir,n=name))
                # don't make all the plots... c.SaveAs('{pd}/cosTheta_chi2_{n}.png'.format(pd=plotsdir,n=name))
    
    
    chi2dist = ROOT.TH1F('chi2', '#chi^{{2}} distribution - {n}'.format(n=('NLO' if options.doNLO else 'LO')), 50, 0., 5.); chi2dist.SetLineWidth(2); chi2dist.SetLineColor(ROOT.kAzure-2)
    chi2dist.GetXaxis().SetTitle('#chi^{2}')
    chi2dist.GetYaxis().SetTitle('# of fits')
    for n,chi2 in allfits.items():
        chi2dist.Fill( min(chi2, chi2dist.GetXaxis().GetBinCenter(chi2dist.GetXaxis().GetLast())) )
    chi2dist.Draw()
    lat.DrawLatex(0.7, 0.6, 'mean #chi^{{2}}: {m:.2f}'.format(m=chi2dist.GetMean()))
    
    c.SaveAs('{pd}/chi2_chi2_pdfWgt{w}.pdf'.format(pd=plotsdir,w=pdfWgtString))
    c.SaveAs('{pd}/chi2_chi2_pdfWgt{w}.png'.format(pd=plotsdir,w=pdfWgtString))
    
        
    ## don't symmetrize for absY (fractionR_plus_sym  , fractionL_plus_sym , fraction0_plus_sym ) = symmetrizeFractions(copy.deepcopy(fractionR_plus ), copy.deepcopy(fractionL_plus ), copy.deepcopy(fraction0_plus ))
    ## don't symmetrize for absY (fractionR_minus_sym , fractionL_minus_sym, fraction0_minus_sym) = symmetrizeFractions(copy.deepcopy(fractionR_minus), copy.deepcopy(fractionL_minus), copy.deepcopy(fraction0_minus))
    
    ## ## this doesn't do anything
    ## for i in range(0):
    ##     fractionR_plus_sym .Smooth(1, 'k5b')
    ##     fractionL_plus_sym .Smooth(1, 'k5b')
    ##     fraction0_plus_sym .Smooth(1, 'k5b')
    ##     fractionR_minus_sym.Smooth(1, 'k5b')
    ##     fractionL_minus_sym.Smooth(1, 'k5b')
    ##     fraction0_minus_sym.Smooth(1, 'k5b')
    
        
    c.SetLeftMargin  (0.15)
    c.SetRightMargin (0.15)
    c.SetBottomMargin(0.15)
    
    fractionR_plus .GetZaxis().SetRangeUser(0., 0.5)
    fractionR_minus.GetZaxis().SetRangeUser(0., 0.5)
    fractionL_plus .GetZaxis().SetRangeUser(0., 0.8)
    fractionL_minus.GetZaxis().SetRangeUser(0., 0.8)
    fraction0_plus .GetZaxis().SetRangeUser(0., 0.4)
    fraction0_minus.GetZaxis().SetRangeUser(0., 0.4)
    
    fractionR_plus .Draw('colz')
    c.SaveAs('{pd}/fractionR_plus_pdf{w}.pdf'.format(pd=plotsdir,w=pdfWgtString))
    c.SaveAs('{pd}/fractionR_plus_pdf{w}.png'.format(pd=plotsdir,w=pdfWgtString))
    fractionL_plus .Draw('colz')
    c.SaveAs('{pd}/fractionL_plus_pdf{w}.pdf'.format(pd=plotsdir,w=pdfWgtString))
    c.SaveAs('{pd}/fractionL_plus_pdf{w}.png'.format(pd=plotsdir,w=pdfWgtString))
    fraction0_plus .Draw('colz')
    c.SaveAs('{pd}/fraction0_plus_pdf{w}.pdf'.format(pd=plotsdir,w=pdfWgtString))
    c.SaveAs('{pd}/fraction0_plus_pdf{w}.png'.format(pd=plotsdir,w=pdfWgtString))
    fractionR_minus.Draw('colz')
    c.SaveAs('{pd}/fractionR_minus_pdf{w}.pdf'.format(pd=plotsdir,w=pdfWgtString))
    c.SaveAs('{pd}/fractionR_minus_pdf{w}.png'.format(pd=plotsdir,w=pdfWgtString))
    fractionL_minus.Draw('colz')
    c.SaveAs('{pd}/fractionL_minus_pdf{w}.pdf'.format(pd=plotsdir,w=pdfWgtString))
    c.SaveAs('{pd}/fractionL_minus_pdf{w}.png'.format(pd=plotsdir,w=pdfWgtString))
    fraction0_minus.Draw('colz')
    c.SaveAs('{pd}/fraction0_minus_pdf{w}.pdf'.format(pd=plotsdir,w=pdfWgtString))
    c.SaveAs('{pd}/fraction0_minus_pdf{w}.png'.format(pd=plotsdir,w=pdfWgtString))
    
    fractionLmR_plus = fractionL_plus .Clone('fractionLmR_plus_pdf{w}' .format(w=pdfWgtString))
    fractionLmR_minus= fractionL_minus.Clone('fractionLmR_minus_pdf{w}'.format(w=pdfWgtString))
    fractionLmR_plus .SetTitle('W^{+}: fractions L-R'+typeString)
    fractionLmR_minus.SetTitle('W^{-}: fractions L-R'+typeString)
    fractionLmR_plus .Add(fractionR_plus ,-1.)
    fractionLmR_minus.Add(fractionR_minus,-1.)
    fractionLmR_plus .GetZaxis().SetRangeUser(0., 0.8)
    fractionLmR_minus.GetZaxis().SetRangeUser(0., 0.8)
    
    fractionLmR_plus.Draw('colz')
    c.SaveAs('{pd}/fractionLmR_plus_pdf{w}.pdf'.format(pd=plotsdir, w=pdfWgtString))
    c.SaveAs('{pd}/fractionLmR_plus_pdf{w}.png'.format(pd=plotsdir, w=pdfWgtString))
    fractionLmR_minus.Draw('colz')
    c.SaveAs('{pd}/fractionLmR_minus_pdf{w}.pdf'.format(pd=plotsdir, w=pdfWgtString))
    c.SaveAs('{pd}/fractionLmR_minus_pdf{w}.png'.format(pd=plotsdir, w=pdfWgtString))
    
    fractionR_plus.Write()
    fractionL_plus.Write()
    fraction0_plus.Write()
    fractionR_minus.Write()
    fractionL_minus.Write()
    fraction0_minus.Write()
    
    ## fractionR_plus_sym .Write()
    ## fractionL_plus_sym .Write()
    ## fraction0_plus_sym .Write()
    ## fractionR_minus_sym.Write()
    ## fractionL_minus_sym.Write()
    ## fraction0_minus_sym.Write()
    
    #h_tmp_wy.Write()
    
outfile.Close()


