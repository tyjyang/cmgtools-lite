import ROOT,os,re
import array
ROOT.gROOT.SetBatch(True)

BINSX = [];
BINSY = [];

def makeParametersHisto(name):
    nFuncPars = 5
    nBinsZ = nFuncPars + 2 # 1 underflow + 5 function parameters + 1 overflow
    BINSZ = range(nBinsZ+1)
    parsh = ROOT.TH3F(name,'',
                      len(BINSX)-1,array.array('f',BINSX),
                      len(BINSY)-1,array.array('f',BINSY),
                      len(BINSZ)-1,array.array('f',BINSZ))
    parsh.GetZaxis().SetBinLabel(1,'norad')
    parsh.GetZaxis().SetBinLabel(nBinsZ,'rad outside the cone')
    for ipar in xrange(1,nFuncPars):
        parsh.GetZaxis().SetBinLabel(ipar+1,'par_{ipar}'.format(ipar=ipar))
    return parsh

def makeCanvas():
    ROOT.gStyle.SetOptStat(0)
    c2 = ROOT.TCanvas('foo','', 600, 600)
    c2.Range(0,0,1,1);
    c2.SetFillColor(0);
    c2.SetBorderMode(0);
    c2.SetBorderSize(2);
    c2.SetTickx(1);
    c2.SetTicky(1);
    c2.SetLeftMargin(0.03);
    c2.SetRightMargin(0.05);
    c2.SetTopMargin(0.05);
    c2.SetBottomMargin(0.0);
    c2.SetFrameFillStyle(0);
    c2.SetFrameBorderMode(0);
    return c2

def makePads(subpadsYLow,subpadsYUp):
    pads = []
    for ip in xrange(len(subpadsYLow)):
        pad = ROOT.TPad("pad{ip}".format(ip=ip),"",0.,subpadsYLow[ip],1.,subpadsYUp[ip])
        pad.SetFillColor(0);
        pad.SetBorderMode(0);
        pad.SetBorderSize(2);
        pad.SetTickx(1);
        pad.SetTicky(1);
        pad.SetLeftMargin(0.2);
        pad.SetRightMargin(0.05);
        pad.SetFrameFillStyle(0);
        pad.SetFrameBorderMode(0);
        pad.SetFrameFillStyle(0);
        pad.SetFrameBorderMode(0);
        pads.append(pad)
    return pads


def makeRatios(args):
    
    infiles = args[0:2]
    outfile = args[2]

    ## take the histograms from the files
    histos = {}
    for filein in infiles:
        tf = ROOT.TFile(filein)
        for k in tf.GetListOfKeys():
            name = k.GetName()
            h = tf.Get(name)
            if not h.InheritsFrom("TH1"): continue
            h.SetDirectory(None)
            if name in histos:
                histos[name].append(h)
            else:
                histos[name] = [h]
            
    ## make ratios of each variable
    for var,histoarr in histos.iteritems():
        for h in histoarr:
            h.Sumw2()
            if h.GetNbinsX()==1000:
                if var=='lptDressOverPreFSR': h.Rebin(4)
                else: h.Rebin(10)
            elif h.GetNbinsX()==5000: h.Rebin(50)
            else: h.Rebin(1)
            h.Scale(1./h.Integral())
            h.SetTitle('')
            if var=='lptDressOverPreFSR': h.GetXaxis().SetRangeUser(0.95,1.5)
            if var.endswith('lpt'): h.GetXaxis().SetRangeUser(0,100)
            
        photos,pythia = histoarr
        photos.SetName(var+'_photos')
        pythia.SetName(var+'_pythia')
        ratio = photos.Clone(var+'_ratio')
        ratio.SetDirectory(None)
        ratio.Divide(pythia)
        histoarr.append(ratio)

    ## save everything in one file
    toutfile = ROOT.TFile(outfile,'recreate')
    for var,histoarr in histos.iteritems():
        for h in histoarr:
            h.Write()

    toutfile.Close()

    return histos

def fitRatioDeltaR(photos,pythia,ratio,varName):
    ## error function
    erf = ROOT.TF1("modErf","[0]*TMath::Erf((x-[1])/[2])+[3]*(x>[4])*x",0.0,0.1)
    erf.SetParameter(0,1.2)
    erf.SetParameter(1,0.05)
    erf.SetParameter(2,0.01)
    erf.SetParLimits(3,1.1,1.9)
    erf.SetParLimits(4,0.006,0.015)
    erf.SetParLimits(2,0.001,0.02)
    ratio.Fit("modErf","R")

    pars = [erf.GetParameter(ipar) for ipar in xrange(0,5)]
    errs = [erf.GetParError(ipar) for ipar in xrange(0,5)]
    return (pars,errs)
    
def plotRatios(histos,axislabels,name):
    
    subpadsYLow = [0.4,0.06]
    subpadsYUp = [0.95,0.4]

    binnedVar = (type(histos.keys()[0])!=str)

    if binnedVar:
        photos_rwgt = ROOT.TFile("photos_rwgt_{name}.root".format(name=name),'recreate')
        th3pars = makeParametersHisto(name)
    
    for var,histoarr in histos.iteritems():
        print "var = ",var," histoarr = ",histoarr
        photos,pythia,ratio = histoarr

        ## two use cases: keys are either the var name in the 1d plotting case of all the variables
        ## or keys are (ieta,ipt) in the case of binned variable plotting
        if not binnedVar: varName = var
        else:
            ieta,ipt = var
            varName = 'ieta{ieta}_ipt{ipt}'.format(ieta=ieta,ipt=ipt)
            # these are the no-radiation cases 
            th3pars.SetBinContent(ieta,ipt,1,photos.GetBinContent(1)/pythia.GetBinContent(1))
            # now the overflows, i.e. the cases where the radiation is outside the cone
            ovfBin = pythia.GetNbinsX()+1
            th3pars.SetBinContent(ieta,ipt,7,photos.GetBinContent(ovfBin)/pythia.GetBinContent(ovfBin))
            
        canv = makeCanvas()
        pads = makePads(subpadsYLow,subpadsYUp)
        canv.cd()
        ## draw the single ones
        pads[0].Draw(); pads[0].cd(); ROOT.gPad.SetBottomMargin(0);
        if re.match('fsrpt.*',varName) or varName=='lptDressOverPreFSR':
            pads[0].SetLogy()

        photos.SetMarkerColor(ROOT.kRed)
        photos.SetLineColor(ROOT.kRed)
        pythia.SetLineColor(ROOT.kBlack)
        photos.SetMaximum(1.1*max(photos.GetMaximum(),pythia.GetMaximum()))
        photos.Draw('pe1')
        pythia.Draw('hist same')
        photos.GetYaxis().SetTitle('Normalized entries')
        photos.GetYaxis().SetTitleSize(0.07)
        photos.GetYaxis().SetDecimals()
        
        canv.cd()
        legend = ROOT.TLegend(0.25, 0.9, 0.90, 0.99);
        legend.SetFillStyle(0); legend.SetBorderSize(0)
        legend.SetNColumns(2)
        legend.AddEntry(photos,'PHOTOS','pl')
        legend.AddEntry(pythia,'PYTHIA 8','pl')
        legend.Draw()

        ## draw the ratio
        pads[1].Draw(); pads[1].cd(); ROOT.gPad.SetTopMargin(0); ROOT.gPad.SetBottomMargin(0.3); 
        ratio.SetMarkerColor(ROOT.kRed+1)
        ratio.SetLineColor(ROOT.kRed+1)
        #ratio.GetYaxis().SetRangeUser(0.95,1.05)
        ratio.GetYaxis().SetTitle('PHOTOS / PYTHIA 8')
        ratio.GetYaxis().SetTitleOffset(0.7)
        ratio.GetYaxis().SetTitleSize(0.08)
        ratio.GetYaxis().SetLabelSize(0.07)
        ratio.GetXaxis().SetLabelSize(0.09)
        ratio.GetYaxis().CenterTitle()
        ratio.GetYaxis().SetDecimals()
        ratio.Draw('pe1')

        # DeltaR(gamma-hard/lep) used for systematic reweighting
        if varName=='fsrdr_hard' or 'ieta' in varName:
            pars,errs = fitRatioDeltaR(photos,pythia,ratio,varName)
            if binnedVar:
                for ipar in xrange(len(pars)):
                    th3pars.SetBinContent(ieta,ipt,ipar+2,pars[ipar])
                    th3pars.SetBinError(ieta,ipt,ipar+2,errs[ipar])
            
        ratio.GetXaxis().SetTitleSize(0.12)
        ratio.GetXaxis().SetTitleOffset(0.9)
        if not binnedVar:
            if varName in axislabels:
                ratio.GetXaxis().SetTitle(axislabels[varName])
            else:
                ratio.GetXaxis().SetTitle(varName)
        else:
             ratio.GetXaxis().SetTitle('Delta R(l,#gamma_{hardest})') # ok, should be passed as argument...
        for ext in ['png','pdf']:
            canv.SaveAs("{var}_{name}.{ext}".format(var=varName,name=name,ext=ext))

    if binnedVar:
        photos_rwgt.cd()
        th3pars.Write()
        photos_rwgt.Close()

        
def getBinnedVar(varname,args):
    infiles = args[0:2]
    th3s = []
    for filein in infiles:
        tf = ROOT.TFile(filein)
        h = tf.Get(varname)
        if not h.InheritsFrom("TH1"):
            print "ERROR: Variable ",varname," is not a TH3F. "
            return
        h.SetDirectory(None)
        h.Sumw2()
        th3s.append(h)

    print "list of TH3F: ",th3s

    del BINSX[:]; del BINSY[:] # just to be sure
    for i in xrange(th3s[0].GetNbinsX()+1):
        BINSX.append(th3s[0].GetXaxis().GetXbins().At(i))
    for i in xrange(th3s[0].GetNbinsY()+1):
        BINSY.append(th3s[0].GetYaxis().GetXbins().At(i))

    print "binsy here = ",BINSY
    
    photos,pythia = th3s
    histos1d = {}
    nbinsz = th3s[0].GetNbinsZ()
    binsz = [th3s[0].GetZaxis().GetXbins().At(i) for i in xrange(nbinsz+1)]
    print "bins z = ",binsz
    for ieta in xrange(1,th3s[0].GetNbinsX()+1):
        for ipt in xrange(1,th3s[0].GetNbinsY()+1):
            photos1d = ROOT.TH1F('photos_ieta{ieta}_ipt{ipt}'.format(ieta=ieta,ipt=ipt),'',
                                 len(binsz)-1,array.array('f',binsz))
            pythia1d = ROOT.TH1F('pythia_ieta{ieta}_ipt{ipt}'.format(ieta=ieta,ipt=ipt),'',
                                 len(binsz)-1,array.array('f',binsz))
            photos1d.SetDirectory(None); pythia1d.SetDirectory(None)
            photos1d.GetXaxis().SetRangeUser(0,0.1); pythia1d.GetXaxis().SetRangeUser(0,0.1)
            ## normalization needs to be conserved in each preFSR bin
            photos1d.Sumw2(); pythia1d.Sumw2()
            # print "ieta= ",ieta,"  ipt= ",ipt
            for iz in xrange(1,nbinsz+2): # include the overflow bin !
                photos1d.SetBinContent(iz,photos.GetBinContent(ieta,ipt,iz))                
                photos1d.SetBinError(iz,photos.GetBinError(ieta,ipt,iz))                
                pythia1d.SetBinContent(iz,pythia.GetBinContent(ieta,ipt,iz))                
                pythia1d.SetBinError(iz,pythia.GetBinError(ieta,ipt,iz))                
                # if ieta==1 and ipt==1:
                #     print "Bin of dr histo = ",iz," has photos = ",photos.GetBinContent(ieta,ipt,iz)," pythia = ",pythia.GetBinContent(ieta,ipt,iz)
            # print "ph int = ",photos1d.Integral(),"   ",pythia1d.Integral()
            photos1d.Scale(1./photos1d.Integral()); pythia1d.Scale(1./pythia1d.Integral())
            ratio1d = photos1d.Clone('photos2pythia_ieta{ieta}_ipt{ipt}'.format(ieta=ieta,ipt=ipt))
            ratio1d.SetDirectory(None)
            ratio1d.Divide(pythia1d)
            ratio1d.GetXaxis().SetRangeUser(0,0.1)
            histos1d[(ieta,ipt)] = [photos1d]
            histos1d[(ieta,ipt)].append(pythia1d)
            histos1d[(ieta,ipt)].append(ratio1d)

    print "list of histogram arrays",histos1d
    return histos1d



if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] photos.root pythia8.root photosOverPythia8.root ")
    parser.add_option('', '--make'   , type='string'       , default='all' , help='run all (default) or only parts (plots, rwgt)')
    parser.add_option("-n", "--name",     dest="name", type="string", default="wp_mu", help="Channel naming (should be wp_mu', 'wm_mu', 'wp_el', 'wm_el'");
    (options, args) = parser.parse_args()
    if len(args)<3:
        print "Need photos.root pythia8.root photosOverPythia8.root"
        exit(1)

    possible_names = ['wp_mu', 'wm_mu', 'wp_el', 'wm_el']
    if options.name not in possible_names:
        print "name should be one among ",possible_names
        exit(1)

    axis_labels = {'leta':'lepton #eta', 'lpt': 'lepton p_{T} [GeV]',
                   'wy':'W rapidity', 'wpt': 'W p_{T} [GeV]', 'wmass': 'W mass',
                   'fsrdr_close': '#Delta R(l,#gamma_{closest})', 'fsrpt_close': 'p_{T}(#gamma_{closest})',
                   'fsrdr_hard': '#Delta R(l,#gamma_{hardest})', 'fsrpt_hard': 'p_{T}(#gamma_{hardest})', 'fsrptfrac_hard': 'p_{T}(#gamma_{hardest})/p_{T}^{l}',
                   'nfsr': 'number of FSR #gamma',
                   'lptDressOverPreFSR': 'p_{T}^{l} / p_{T}^{pre-FSR l}'}

    if options.make in ['all','plots']:
        histos = makeRatios(args)
        plotRatios(histos,axis_labels,options.name)

    if options.make in ['all','rwgt']:
        binnedVarRatios = getBinnedVar("h3d_fsrdr_hard",args)
        plotRatios(binnedVarRatios,axis_labels,options.name)
