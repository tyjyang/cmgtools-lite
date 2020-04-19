import ROOT
from math import *
import re, os, glob
from array import array
from CMGTools.WMass.plotter.tree2yield import scalarToVector
from CMGTools.TTHAnalysis.tools.plotDecorations import doSpam

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

if "/functions_cc.so" not in ROOT.gSystem.GetLibraries(): 
    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/plotter/functions.cc+" % os.environ['CMSSW_BASE']);
if "/w-mass-13TeV/functionsWMass_cc.so" not in ROOT.gSystem.GetLibraries(): 
    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/plotter/w-mass-13TeV/functionsWMass.cc+" % os.environ['CMSSW_BASE']);

def setAxisTH2(h2,options,xtitle=None,ytitle=None):

    charge = ""
    if options.selectCharge == "plus": charge = "positive"
    elif options.selectCharge == "minus": charge = "negative"

    h2.SetTitle("")
    h2.SetMarkerSize(1.2)            
    # X
    h2.GetXaxis().SetTitle(("%s muon #eta" % charge) if xtitle==None else xtitle)
    h2.GetXaxis().SetTitleSize(0.05)
    h2.GetXaxis().SetLabelSize(0.04)
    h2.GetXaxis().SetTitleOffset(0.95) 
    # Y
    h2.GetYaxis().SetTitle(("%s muon p_{T} (GeV)" % charge) if ytitle==None else ytitle)    
    h2.GetYaxis().SetTitleSize(0.05)
    h2.GetYaxis().SetLabelSize(0.04)
    h2.GetYaxis().SetTitleOffset(0.95)
    # Z
    h2.GetZaxis().SetTitleSize(0.05)
    h2.GetZaxis().SetLabelSize(0.04)
    h2.GetZaxis().SetTitleOffset(1.2) # 1.4


def getBinsFromName(name):

    # name = "hist_pt_%.1f_%.1f_eta_%.3f_%.3f"
    if not all([x in name for x in ["eta","pt"]]):
        #print "Error in getBinsFromName(): 'eta' or 'pt' not found in %s. Exit" % name
        return -999.,-999.,-999.,-999.
    tokens = name.split('_')
    for i,tkn in enumerate(tokens):
        #print "%d %s" % (i, tkn)                                                                                      
        if tkn == "eta": 
            etalow  = float(tokens[i + 1].replace('m','-').replace('p','.'))
            etahigh = float(tokens[i + 2].replace('m','-').replace('p','.'))
        if tkn == "pt":  
            ptlow  = float(tokens[i + 1].replace('m','-').replace('p','.'))
            pthigh = float(tokens[i + 2].replace('m','-').replace('p','.'))
    return ptlow,pthigh,etalow,etahigh

def makeSignalModel(model, w, useRightTailCB=False):
    if model == "Z-Voit":
        w.factory("Voigtian::signal(mass,sum::(dm[0,-10,10],MZ[91.1876]), GammaZ[2.495],sigma)")
        return (w.pdf("signal"), ["dm", "sigma"])
    if model == "Z-CB":
        w.factory("BreitWigner::zBW(mass,MZ[91.1876], GammaZ[2.495])")
        if useRightTailCB:
            w.factory("CBShape::resol(mass,dm[0,-10,10], sigma, alpha[-3., -8.0, -0.5], n[1, 0.1, 100.])")
        else:
            w.factory("CBShape::resol(mass,dm[0,-10,10], sigma, alpha[3., 0.5, 8], n[1, 0.1, 100.])")            
        w.factory("FCONV::signal(mass,zBW,resol)")
        return (w.pdf("signal"), ["dm", "sigma","alpha", "n"] )
    if model == "Z-DCB":
        #ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        #w.factory("BreitWigner::zBW(mass,MZ[91.1876], GammaZ[2.495])")
        ROOT.gSystem.Load("%s/src/CMGTools/WMass/python/plotter/testDoubleCB/my_double_CB_cc.so" % os.environ['CMSSW_BASE'])
        w.factory("BreitWigner::zBW(mass,MZ[91.1876], GammaZ[2.495])")
        w.factory("My_double_CB::resol(mass,dm[0,-10,10], sigma, alpha[3., 0.5, 8], n[1, 0.1, 100.], alpha2[3., 0.5, 8], n2[1, 0.1, 100.])")
        w.factory("FCONV::signal(mass,zBW,resol)")
        return (w.pdf("signal"), ["dm", "sigma","alpha", "n","alpha2", "n2"] )
    raise RuntimeError, "missing signal model: "+model

def makeSignalModel1D(model, w, useRightTailCB=False):
    w.factory("sigma[1.5,0.3,10]")
    return makeSignalModel(model,w, useRightTailCB=useRightTailCB)

def makeSignalModel2D(model, w):
    #w.factory('expr::sigma("@0*exp(@1)",lambda[1,0.5,2.0], logdm)')
    w.factory('prod::sigma(lambda[1,0.5,2.0], massErr)')
    (pdf,params) = makeSignalModel(model,w)
    params[params.index("sigma")] = "lambda"
    return (pdf,params)


def makeBackgroundModel(model, w):
    if model == "Expo":
        w.factory("Exponential::background(mass,slope[-0.1,-1,0.5])")
        return (w.pdf("background"), [])
    if model == "None":
        return (None, [])
    raise RuntimeError, "missing signal model: "+model

def makeSumModel(spdf,bpdf,w):
    if not bpdf: return spdf
    w.factory("SUM::sumModel(f_signal[0.9,0.5,1]*%s, %s)" % (spdf.GetName(), bpdf.GetName()))
    return w.pdf("sumModel")

def printPlot(frame, name, text, options, spam=[], canvas=None, legend=None):
    c1 = canvas if canvas != None else ROOT.TCanvas("c1","c1",1200,900)
    c1.SetTopMargin(0.08)
    c1.SetTickx(1)
    c1.SetTicky(1)
    xoffs = 0;
    frame.SetMaximum(frame.GetMaximum()*(1.05+0.07*len(spam)))
    frame.SetTitle("")
    c1.SetLeftMargin(0.12)
    if frame.GetMaximum() >= 100*1000:
        xoffs = 0.155
        frame.Draw();
        frame.GetYaxis().SetTitleOffset(1.7)
    elif frame.GetMaximum() >= 1000:
        xoffs = 0.02
        frame.Draw()
        frame.GetYaxis().SetTitleOffset(1.45)
    else:        
        frame.Draw()
        frame.GetYaxis().SetTitleOffset(1.25)
    c1.SetRightMargin(0.05); 
    frame.GetXaxis().SetNdivisions(505);   
    frame.GetXaxis().SetDecimals(1);   
    frame.GetXaxis().SetTitle("dilepton mass (GeV)");   
    frame.GetXaxis().SetTitleSize(0.05);   
    frame.GetXaxis().SetTitleOffset(1.02);       
    if legend != None: legend.Draw("same")
    printCanvas(c1, name, text, options, xoffs=xoffs, spam=spam)

def printCanvas(c1, name, text, options, xoffs=0, spam=[]):
    if   options.lumi > 3.54e+1: lumitext = "%.1f fb^{-1}" % options.lumi
    elif options.lumi > 3.54e+0: lumitext = "%.1f fb^{-1}" % options.lumi
    elif options.lumi > 3.54e-1: lumitext = "%.2f fb^{-1}" % options.lumi
    elif options.lumi > 3.54e-2: lumitext = "%.0f pb^{-1}" % (options.lumi*1000)
    elif options.lumi > 3.54e-3: lumitext = "%.1f pb^{-1}" % (options.lumi*1000)
    else                       : lumitext = "%.2f pb^{-1}" % (options.lumi*1000)

    doSpam("#bf{CMS} #it{%s}" % ("Simulation" if "_ref" in name else "Preliminary"), 0.15+xoffs, 0.94, 0.55+xoffs, 0.99, align=12, textSize=options.textSize)
    doSpam(lumitext+" (13 TeV)", .55+xoffs, .94, .97, .99, align=32, textSize=options.textSize)

    if spam:
        y0 = 0.88 - options.textSize*1.2
        for line in spam:
            niceline = re.sub(r"(\s+)-(\d+)",r"\1#minus\2", line)
            doSpam(niceline, 0.55, y0, 0.95, y0 + options.textSize*1.2, textSize=options.textSize, align=22)
            y0 -= options.textSize * 1.2

    lat = ROOT.TLatex()
    lat.SetTextSize(0.04)
    lat.SetTextFont(42)
    lat.SetTextColor(ROOT.kBlack)
    lat.SetNDC()
    ptlow,pthigh,etalow,etahigh = getBinsFromName(name)
    if ptlow > 0.0:
        lat.DrawLatex(0.15, 0.5,  "{ptlow:3g} < p_{{T}}(GeV) < {pthigh:3g}".format(ptlow=ptlow,pthigh=pthigh))
        lat.DrawLatex(0.15, 0.44, "{etalow:3g} < #eta < {etahigh:3g}".format(etalow=etalow,etahigh=etahigh))

    c1.RedrawAxis("sameaxis")
    exts = [x for x in options.plotExtension.split(',')]
    for ext in exts:
        c1.Print("%s/%s.%s" % (options.printDir, name, ext))
    log = open("%s/%s.%s"  % (options.printDir, name, "txt"), "w")
    for line in text: log.write(line +"\n")
    log.close()


def makePlot1D(w, data, pdf, params, result, name, options, histRoot=None):

    xmin = float(options.xvar[1].split(",")[1])
    xmax = float(options.xvar[1].split(",")[2])
    if options.fitRange:
        xmin = float(max(options.fitRange[0],xmin))
        xmax = float(min(options.fitRange[1],xmax))

    frame = w.var("mass").frame(ROOT.RooFit.Bins(w.var("mass").getBins()/options.rebin))
    data.plotOn(frame, ROOT.RooFit.Name("points"))
    if pdf.GetName() == "sumModel":
        pdf.plotOn(frame, ROOT.RooFit.Components("background"), ROOT.RooFit.LineColor(ROOT.kRed+1),ROOT.RooFit.LineStyle(2), ROOT.RooFit.Name("pdf_bkg"))
    pdf.plotOn(frame, ROOT.RooFit.Name("pdf_model"))

    chi2 = frame.chiSquare("pdf_model","points",result.floatParsFinal().getSize()) # this might be wrong if range is not set for the curve, but it should be already set since the fit is made in the desired range
    #
    # alternate Chi2, but it is different from the one above
    chi2_v2 = ROOT.RooChi2Var("chi2","chi2 var",pdf,data, 0, "fit");
    ndof = -1
    nbinsInRange = -1
    if histRoot:
        # Ndof = number of bins in range - number of parameters - 1
        # add epsilon = 0.0001 to value when looking for bin.
        # e.g. 3 bins from 70 to 73 means than 70.0001 and 73.0001 are bin number 1 and 4,so the difference gives 3
        ndof = histRoot.GetXaxis().FindFixBin(xmax+0.0001) - histRoot.GetXaxis().FindFixBin(xmin+0.0001)  - result.floatParsFinal().getSize() - 1
        nbinsInRange = histRoot.GetXaxis().FindFixBin(xmax+0.0001) - histRoot.GetXaxis().FindFixBin(xmin+0.0001)

    c1 = ROOT.TCanvas("c1","c1",1200,900)
    legfit = ROOT.TLegend(0.15, 0.72, 0.5, 0.9)
    legfit.SetFillColor(0)
    legfit.SetFillStyle(0)
    legfit.SetBorderSize(0)
    dummy1 = ROOT.TH1F("hdummy1","",1,0,1)
    dummy2 = ROOT.TH1F("hdummy2","",1,0,1)
    dummy3 = ROOT.TH1F("hdummy3","",1,0,1)
    dummy1.SetMarkerStyle(20)
    dummy2.SetLineColor(ROOT.kBlue);    dummy2.SetLineWidth(2);
    dummy3.SetLineColor(ROOT.kRed+1);    dummy3.SetLineWidth(2);  dummy3.SetLineStyle(2)
    if "ref" in name : legfit.AddEntry(dummy1,"MC data","LPE")
    else:              legfit.AddEntry(dummy1,"data","LPE")
    legfit.AddEntry(dummy2,"model","L")
    if pdf.GetName() == "sumModel":
        legfit.AddEntry(dummy3,"background","L")

    text = []
    for param in params:
        var = result.floatParsFinal().find(param)
        text.append("%s: %.3f +- %.3f" % (param, var.getVal(), var.getError()))
    if "dm" in params and "sigma" in params:
        dm = result.floatParsFinal().find("dm")
        sm = result.floatParsFinal().find("sigma")
        spam = [ ("#Delta = %.3f #pm %.3f GeV" % (dm.getVal(), dm.getError())),
                 ("#sigma = %.3f #pm %.3f GeV" % (sm.getVal(), sm.getError())) ]
    else:
        spam = []
    spam.append("#chi^{2} = %.2f   N^{model}_{pars} = %d" % (chi2, result.floatParsFinal().getSize()))
    #spam.append("#chi^{2}_{alt} = %.2f   ndf = %d" % (chi2_v2.getVal(), ndof))
    #spam.append("N_{bins}(%3.1f,%3.1f) = %d" % (xmin,xmax,nbinsInRange))
    if w.function("resol"):
        foms = getResolutionFOMs(w)
        for k,v in foms.iteritems(): 
            text.append("%s: %.3f" % (k,v))
    printPlot(frame,name,text,options, spam = spam, canvas=c1, legend=legfit)

class GioProjWData:
    def __init__(self, pdf, xvar, data, rebinFactor=1):
        self.xvar = xvar
        self.pdfOriginal = pdf
        self.dataOriginal = data
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING);
        vars = ROOT.RooArgSet(self.dataOriginal.get());  
        vars.remove(xvar,False,True);
        pars = pdf.getParameters(ROOT.RooArgSet(xvar));
        hist = None;
        if data.ClassName() == "RooDataSet":
            others = ROOT.RooDataSet("","dm",data,vars);
        else:
            others = ROOT.RooDataHist("","dm",vars,data);
        print "Will have to average over %d entries " % others.numEntries();
        for i in xrange(others.numEntries()):
            pars.assignValueOnly(others.get(i))
            h = pdf.createHistogram("htemp",xvar,ROOT.RooFit.Binning(rebinFactor*xvar.getBins()));
            h.SetDirectory(0);
            h.Scale(others.weight()/h.Integral());
            if hist == None: hist = h 
            else:            hist.Add(h); 
        self.hist = hist
        self.rdh  = ROOT.RooDataHist("","",ROOT.RooArgList(xvar),self.hist);
        self.pdf  = ROOT.RooHistPdf("","",ROOT.RooArgSet(xvar),self.rdh,2);
    def plotOn(self, frame, *args):
        return self.pdf.plotOn(frame, *args)

def makePlot2D(w, data, pdf, params, result, name, options):
    frame = w.var("mass").frame(ROOT.RooFit.Bins(w.var("mass").getBins()/options.rebin))
    data.plotOn(frame)
    #if pdf.GetName() == "sumModel":
    #    pdf.plotOn(frame, ROOT.RooFit.Components("background"), ROOT.RooFit.LineStyle(2))
    proj = GioProjWData(pdf, w.var("mass"), data)
    proj.plotOn(frame)
    text = []
    for param in params:
        var = result.floatParsFinal().find(param)
        text.append("%s: %.3f +- %.3f" % (param, var.getVal(), var.getError()))
    if "dm" in params and "lambda" in params:
        dm = result.floatParsFinal().find("dm")
        sm = result.floatParsFinal().find("lambda")
        spam = [ ("#Delta = %.2f #pm %.2f GeV" % (dm.getVal(), dm.getError())),
                 ("#lambda = %.2f #pm %.2f GeV" % (sm.getVal(), sm.getError())) ]
    else:
        spam = []
    printPlot(frame,name,text,options, spam = spam)

def makePlot1DRef(w, data, pdf, pdfref, params, result, refresult, name, options):
    frame = w.var("mass").frame()
    data.plotOn(frame, ROOT.RooFit.Name("points_data"))
    if pdfref.GetName() == "sumModel" and pdf.GetName() == "sumModel":
        pdf.plotOn(frame, ROOT.RooFit.Components("background"), ROOT.RooFit.LineColor(ROOT.kRed+1), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Name("pdf_bkg_data"))
        pdf.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.Name("pdf_model_data"))
        pdf.plotOn(frame, ROOT.RooFit.Components("signal"), ROOT.RooFit.LineColor(ROOT.kGreen+2), ROOT.RooFit.Name("pdf_signal_data"))
        pdfref.plotOn(frame, ROOT.RooFit.Components("signal"), ROOT.RooFit.LineColor(ROOT.kGray+1),ROOT.RooFit.Name("pdf_signal_mc"))
    else:
        pdf.plotOn(frame, ROOT.RooFit.Components("signal"), ROOT.RooFit.LineColor(ROOT.kGreen+2), ROOT.RooFit.Name("pdf_signal_data"))
        pdfref.plotOn(frame, ROOT.RooFit.Components("signal"), ROOT.RooFit.LineColor(ROOT.kGray+1),ROOT.RooFit.Name("pdf_signal_mc"))
        
    # dummy legend
    c1 = ROOT.TCanvas("c1","c1",1200,900)
    legfit = ROOT.TLegend(0.15, 0.6, 0.5, 0.9)
    legfit.SetFillColor(0)
    legfit.SetFillStyle(0)
    legfit.SetBorderSize(0)
    dummy1 = ROOT.TH1F("hdummy1","",1,0,1)
    dummy2 = ROOT.TH1F("hdummy2","",1,0,1)
    dummy3 = ROOT.TH1F("hdummy3","",1,0,1)
    dummy4 = ROOT.TH1F("hdummy4","",1,0,1)
    dummy5 = ROOT.TH1F("hdummy5","",1,0,1)
    dummy1.SetMarkerStyle(20)
    dummy2.SetLineColor(ROOT.kBlue);    dummy2.SetLineWidth(2);
    dummy3.SetLineColor(ROOT.kGreen+2);    dummy3.SetLineWidth(2);
    dummy4.SetLineColor(ROOT.kRed+1);    dummy4.SetLineWidth(2);  dummy4.SetLineStyle(2)
    dummy5.SetLineColor(ROOT.kGray+1);    dummy5.SetLineWidth(2);
    if pdfref.GetName() == "sumModel" and pdf.GetName() == "sumModel":
        legfit.AddEntry(dummy1,"data","LPE")
        legfit.AddEntry(dummy2,"model (data)","L")
        legfit.AddEntry(dummy3,"signal (data)","L")
        legfit.AddEntry(dummy4,"background (data)","L")
        legfit.AddEntry(dummy5,"signal (MC)","L")
    else:
        legfit.AddEntry(dummy1,"data","LPE")
        legfit.AddEntry(dummy3,"signal (data)","L")
        legfit.AddEntry(dummy5,"signal (MC)","L")
    
    lat = ROOT.TLatex()
    lat.SetTextSize(0.025)
    lat.SetTextFont(42)
    ptlow,pthigh,etalow,etahigh = getBinsFromName(name)
    lat.DrawLatex(0.18, 0.5,  "{ptlow:3g} < p_{{T}}(GeV) < {pthigh:3g}".format(ptlow=ptlow,pthigh=pthigh))
    lat.DrawLatex(0.18, 0.45, "{etalow:3g} < #eta < {etahigh:3g}".format(etalow=etalow,etahigh=etahigh))

    text = []
    for param in params:
        var = result.floatParsFinal().find(param)
        ref = refresult.floatParsFinal().find(param)
        text.append("%s: fit %.3f +- %.3f ref %.3f +- %.3f  diff %.3f +- %.3f" % (param,
                        var.getVal(), var.getError(), ref.getVal(), ref.getError(),
                        var.getVal()-ref.getVal(), hypot(var.getError(), ref.getError())))
    if "dm" in params and "sigma" in params:
        fdm =    result.floatParsFinal().find("dm")
        fsm =    result.floatParsFinal().find("sigma")
        rdm = refresult.floatParsFinal().find("dm")
        rsm = refresult.floatParsFinal().find("sigma")
        spam = [ "#Delta - #Delta_{MC} = %.3f #pm %.3f GeV" % (fdm.getVal() - rdm.getVal(),
                                                               hypot(fdm.getError(), rdm.getError()))
        ]
        spam.append("#Delta = %.3f #pm %.3f GeV" % (fdm.getVal(),fdm.getError()))
        spam.append("#Delta_{MC} = %.3f #pm %.3f GeV" % (rdm.getVal(),rdm.getError()))
        spam.append("#sigma/#sigma_{MC} = %.3f #pm %.3f" % (fsm.getVal()/rsm.getVal(), (fsm.getVal()/rsm.getVal())*hypot(fsm.getError()/fsm.getVal(), rsm.getError()/rsm.getVal())) )
    else:
        spam = []
    printPlot(frame,name,text,options, spam=spam, canvas=c1, legend=legfit)


def fit1D(hist, options, modelOverride=False, etaptbin=None):

    useRightTailCB = False  # default is to have tail of single-CB to the left
    # manage fit range
    xmin = float(options.xvar[1].split(",")[1])
    xmax = float(options.xvar[1].split(",")[2])
    if options.fitRange:
        xmin = max(options.fitRange[0],xmin)
        xmax = min(options.fitRange[1],xmax)
    if etaptbin:
        ptbin = etaptbin[0] -1
        etabin = etaptbin[1] -1
        ptedges = [float(x) for x in options.ptbins.split(',')]
        etaedges = [float(x) for x in options.etabins.split(',')]
        print "\n"*5
        print "="*40
        print ">>> Fitting bin: %.1f < pt < %.1f   %.1f < eta < %.1f" % (ptedges[ptbin],ptedges[ptbin+1],
                                                                         etaedges[etabin],etaedges[etabin+1])
        print "="*40
        print "\n"

        if options.signalModel == "Z-CB":
            if ptedges[ptbin] >= 45.0:
                useRightTailCB = True
                xmin = 80.0
                xmax = 110.0
            else:
                xmin = 70.0
                xmax = 100.0        
        #if ptedges[ptbin] >= 45.0 and ptedges[ptbin+1] <= 50 and etaedges[etabin] >= 1.6:
            #bpdf, bparams = (None, []) # use no bkg when fit range is very narrow (not enough leverage)

    w = ROOT.RooWorkspace("w")
    w.factory("mass[%g,%g]" % (hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax()))
    w.var("mass").setBins(hist.GetNbinsX())
    w.var("mass").SetTitle(hist.GetXaxis().GetTitle())
    data = ROOT.RooDataHist("data","data", ROOT.RooArgList(w.var("mass")), hist)
    spdf, sparams = makeSignalModel1D(options.signalModel if not modelOverride else modelOverride, w, useRightTailCB=useRightTailCB)
    bpdf, bparams = makeBackgroundModel(options.backgroundModel, w) if options.backgroundModel else (None, [])
    pdf = makeSumModel(spdf, bpdf, w)
    result = pdf.fitTo(data, ROOT.RooFit.Save(True),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Range(float(xmin), float(xmax)),ROOT.RooFit.Strategy(options.fitStrategy),ROOT.RooFit.PrintLevel(options.roofitPrintLevel))
    return (w,data,pdf,sparams+bparams,result)

def getResolutionFOMs(w,pdf="resol",xvar="mass",xrange=[-10,10]):
    w2 = ROOT.RooWorkspace("w2")
    w2.factory("%s[0,%g,%g]" % (xvar, xrange[0], xrange[1]))
    getattr(w2,'import')(w.pdf(pdf), ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence(True))
    w2.var(xvar).setBins(2000)
    hist = w2.pdf(pdf).createHistogram("",w2.var(xvar)); hist.SetDirectory(None)
    imax = hist.GetMaximumBin()
    ymax = hist.GetBinContent(imax)
    ymid = 0.5*ymax
    ysum = ymax; y68 = 0.6827*hist.Integral()
    ihi = ilo = imax
    nbins = hist.GetNbinsX()
    hist.SetBinContent(0,0); hist.SetBinContent(nbins+1,0)
    sigmaEff = -1; fwhm = -1
    while ilo >= 1 or ihi <= nbins:
        ileft, iright = ilo-1, ihi+1
        yleft, yright = map(hist.GetBinContent, (ileft,iright))
        if ilo >= 1 and yleft >= yright:
            ysum += yleft
            ilo   = ileft
        else:
            ysum += yright
            ihi   = iright
        if ysum > y68 and sigmaEff < 0:
            sigmaEff = 0.5*abs(hist.GetXaxis().GetBinLowEdge(ilo) - hist.GetXaxis().GetBinUpEdge(ihi))
        if max(yleft,yright) < ymid and fwhm < 0:
            fwhm = abs(hist.GetXaxis().GetBinLowEdge(ilo) - hist.GetXaxis().GetBinUpEdge(ihi))/(2.*1.177)
    return { 'fwhm23':fwhm, 'sigmaEff':sigmaEff }
    
    
def fit2D(hist, options):
    w = ROOT.RooWorkspace("w")
    w.factory("mass[%g,%g]" % (hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax()))
    w.factory("massErr[%g,%g]" % (hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax()))
    w.var("mass").setBins(hist.GetNbinsX())
    w.var("massErr").setBins(hist.GetNbinsY())
    w.var("mass").SetTitle(hist.GetXaxis().GetTitle())
    w.var("massErr").SetTitle(hist.GetYaxis().GetTitle())
    data = ROOT.RooDataHist("data","data", ROOT.RooArgList(w.var("mass"),w.var("massErr")), hist)
    spdf, sparams = makeSignalModel2D(options.signalModel, w)
    bpdf, bparams = makeBackgroundModel(options.backgroundModel, w) if options.backgroundModel else (None, [])
    pdf = makeSumModel(spdf, bpdf, w)
    conditional = ROOT.RooFit.ConditionalObservables(ROOT.RooArgSet(w.var("massErr")))
    result = pdf.fitTo(data, ROOT.RooFit.Save(True), ROOT.RooFit.SumW2Error(True), conditional)
    return (w,data,pdf,sparams+bparams,result)



def makeCut(cut):
    if cut == "Zee-BB":
        return "abs(LepGood_pdgId[0]) == 11 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) < 1.479"
    elif cut == "Zee-NonBB":
        return "abs(LepGood_pdgId[0]) == 11 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) > 1.479"
    elif cut == "Zee-BE":
        return "abs(LepGood_pdgId[0]) == 11 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) > 1.479 && min(abs(LepGood_eta[0]),abs(LepGood_eta[1])) < 1.479"
    elif cut == "Zee-EE":
        return "abs(LepGood_pdgId[0]) == 11 && min(abs(LepGood_eta[0]),abs(LepGood_eta[1])) > 1.479"
    elif cut == "Zee":
        return "abs(LepGood_pdgId[0]) == 11"
    elif cut == "Zee-dm132":
        return "abs(LepGood_pdgId[0]) == 11 && 1.3 < 0.5*mZ1*TMath::Hypot(LepGood_ptErr[0]/LepGood_pt[0], LepGood_ptErr[1]/LepGood_pt[1]) && 0.5*mZ1*TMath::Hypot(LepGood_ptErr[0]/LepGood_pt[0], LepGood_ptErr[1]/LepGood_pt[1]) < 2.00"
    elif cut == "Zee-wmass":
        return "abs(LepGood_pdgId[0]) == 11 && min(LepGood_eleMVAId[0],LepGood_eleMVAId[1]) > 1 && max(LepGood_relIso04[0],LepGood_relIso04[1]) < 0.4 && min(LepGood_convVetoFull[0],LepGood_convVetoFull[1])==1"
    elif cut == "Zee-binning":
        return "1==1"
    elif cut == "Zmm-BB":
        return "abs(LepGood_pdgId[0]) == 13 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) < 1.2"
    elif cut == "Zmm-NonBB":
        return "abs(LepGood_pdgId[0]) == 13 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) > 1.2"
    elif cut == "Zmm-BE":
        return "abs(LepGood_pdgId[0]) == 13 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) > 1.2 && min(abs(LepGood_eta[0]),abs(LepGood_eta[1])) < 1.2"
    elif cut == "Zmm-EE":
        return "abs(LepGood_pdgId[0]) == 13 && min(abs(LepGood_eta[0]),abs(LepGood_eta[1])) > 1.2"
    elif cut == "Zmm":
        return "abs(LepGood_pdgId[0]) == 13"
    elif cut == "Zmm-B0E1":
        return "abs(LepGood_pdgId[0]) == 13 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) > 1.2 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) < 1.5 && min(abs(LepGood_eta[0]),abs(LepGood_eta[1])) < 1.0"
    elif cut == "Zmm-B0E2":
        return "abs(LepGood_pdgId[0]) == 13 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) > 1.5 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) < 2.0 && min(abs(LepGood_eta[0]),abs(LepGood_eta[1])) < 1.0"
    elif cut == "Zmm-B0E3":
        return "abs(LepGood_pdgId[0]) == 13 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) > 2.0 && max(abs(LepGood_eta[0]),abs(LepGood_eta[1])) < 2.4 && min(abs(LepGood_eta[0]),abs(LepGood_eta[1])) < 1.0"
    elif cut == "Zmm-full":
        return "(passVertexPreSel == 1) && (HLT_BIT_HLT_IsoMu24_v > 0 || HLT_BIT_HLT_IsoTkMu24_v > 0 ) && nLepGood == 2 && (LepGood_pdgId[0] * LepGood_pdgId[1]) == -169 && LepGood_mediumMuonId[0]  > 0 && LepGood_mediumMuonId[1] > 0 && LepGood_relIso04[0] < 0.15 && LepGood_relIso04[1] < 0.15"
    else:
        return cut

def makeHist1D(tree, options, weightedMC = False):
    print "Filling tree for cut %s" % options.cut
    finalCut =  makeCut(options.cut)
    if weightedMC: 
        # add trigger SF for triggering object
        # add also prefiring, but should decide which lepton to use (unless using jet function)
        funcTrigSF = ""
        if options.useAllEvents:
            flag = 0  # need to pass 1 for positive charge, 0 for negative, because function triggerSFforChargedLeptonMatchingTriggerWlike takes as first argument a bool, which is 1 (0) for positive (negative) charge
            # the flag is typically defined using function isOddEvent(evt), which will be true (false) for odd (even) events, for which we generally want charge plus (minus)
            if options.selectCharge == "plus":
                flag = 1
            elif options.selectCharge == "minus":
                flag = 0
            funcTrigSF = "triggerSFforChargedLeptonMatchingTriggerWlike({f},(LepGood_matchedTrgObjMuDR[0]>0.0 || LepGood_matchedTrgObjTkMuDR[0] > 0.0),(LepGood_matchedTrgObjMuDR[1] >0.0 || LepGood_matchedTrgObjTkMuDR[1]>0.0),LepGood_pdgId[0],LepGood_pdgId[1],{pt}[0],{pt}[1],LepGood_eta[0],LepGood_eta[1])".format(f=flag,pt=options.ptvar)
        else:
            funcTrigSF = "triggerSFforChargedLeptonMatchingTriggerWlike(isOddEvent(evt),(LepGood_matchedTrgObjMuDR[0]>0.0 || LepGood_matchedTrgObjTkMuDR[0] > 0.0),(LepGood_matchedTrgObjMuDR[1] >0.0 || LepGood_matchedTrgObjTkMuDR[1]>0.0),LepGood_pdgId[0],LepGood_pdgId[1],{pt}[0],{pt}[1],LepGood_eta[0],LepGood_eta[1])".format(pt=options.ptvar)

        wgt = "puw2016_nTrueInt_36fb(nTrueInt)*_get_muonSF_recoToSelection(LepGood_pdgId[0],{pt}[0],LepGood_eta[1])*_get_muonSF_recoToSelection(LepGood_pdgId[1],{pt}[1],LepGood_eta[1])* {trigSF} * jetPrefireSF_JetAll_Clean".format(pt=options.ptvar, trigSF=funcTrigSF)
        finalCut = "(%s) * %s" % (finalCut, wgt)
    if options.selectCharge in ["plus", "minus"] and not options.useAllEvents:
        finalCut += " * {func}(evt)".format(func="isOddEvent" if options.selectCharge=="plus" else "isEvenEvent") 
    finalCut = scalarToVector(finalCut)
    tree.Draw("%s>>hist(%s)" % (scalarToVector(options.xvar[0]), options.xvar[1]), finalCut, "goff", options.maxEntries)
    print "Done, fitting."
    hist = ROOT.gROOT.FindObject("hist")
    if weightedMC:
        for b in xrange(hist.GetNbinsX()+2): 
            hist.SetBinContent(b, max(0, hist.GetBinContent(b)))
    hist.GetXaxis().SetTitle(options.xtitle)
    return hist

def makeHist2D(tree, options, weightedMC = False, 
               dmvar = "0.5*mZ1*TMath::Hypot(LepGood_ptErrTk[0]/LepGood_pt[0], LepGood_ptErrTk[1]/LepGood_pt[1])",
               dmbins = 50, dmrange=[0.3,5]):
    print "Filling tree for cut %s" % options.cut
    finalCut =  makeCut(options.cut)
    if weightedMC: 
        # add trigger SF for triggering object
        # add also prefiring, but should decide which lepton to use (unless using jet function)
        funcTrigSF = ""
        if options.useAllEvents:
            flag = 0  # need to pass 1 for positive charge, 0 for negative, because function triggerSFforChargedLeptonMatchingTriggerWlike takes as first argument a bool, which is 1 (0) for positive (negative) charge
            # the flag is typically defined using function isOddEvent(evt), which will be true (false) for odd (even) events, for which we generally want charge plus (minus)
            if options.selectCharge == "plus":
                flag = 1
            elif options.selectCharge == "minus":
                flag = 0
            funcTrigSF = "triggerSFforChargedLeptonMatchingTriggerWlike({f},(LepGood_matchedTrgObjMuDR[0]>0.0 || LepGood_matchedTrgObjTkMuDR[0] > 0.0),(LepGood_matchedTrgObjMuDR[1] >0.0 || LepGood_matchedTrgObjTkMuDR[1]>0.0),LepGood_pdgId[0],LepGood_pdgId[1],{pt}[0],{pt}[1],LepGood_eta[0],LepGood_eta[1])".format(f=flag,pt=options.ptvar)
        else:
            funcTrigSF = "triggerSFforChargedLeptonMatchingTriggerWlike(isOddEvent(evt),(LepGood_matchedTrgObjMuDR[0]>0.0 || LepGood_matchedTrgObjTkMuDR[0] > 0.0),(LepGood_matchedTrgObjMuDR[1] >0.0 || LepGood_matchedTrgObjTkMuDR[1]>0.0),LepGood_pdgId[0],LepGood_pdgId[1],{pt}[0],{pt}[1],LepGood_eta[0],LepGood_eta[1])".format(pt=options.ptvar)

        wgt = "puw2016_nTrueInt_36fb(nTrueInt)*_get_muonSF_recoToSelection(LepGood_pdgId[0],{pt}[0],LepGood_eta[1])*_get_muonSF_recoToSelection(LepGood_pdgId[1],{pt}[1],LepGood_eta[1])* {trigSF} * jetPrefireSF_JetAll_Clean".format(pt=options.ptvar, trgiSF=funcTrigSF)
        finalCut = "(%s) * %s" % (finalCut, wgt)
    if options.selectCharge in ["plus", "minus"] and not options.useAllEvents:
        finalCut += " * {func}(evt)".format(func="isOddEvent" if options.selectCharge=="plus" else "isEvenEvent") 
    finalCut = scalarToVector(finalCut)
    tree.Draw("%s:%s>>hist(%s,%d,%g,%g)" % (scalarToVector(dmvar), scalarToVector(options.xvar[0]), scalarToVector(options.xvar[1]), dmbins,dmrange[0],dmrange[1]), finalCut, "goff", options.maxEntries)
    hist = ROOT.gROOT.FindObject("hist")
    if weightedMC:
        for bx in xrange(hist.GetNbinsX()+2): 
          for by in xrange(hist.GetNbinsY()+2): 
            hist.SetBinContent(bx,by, max(0, hist.GetBinContent(bx,by)))
    return hist

def makeHist1D_DMSlices(tree, options, weightedMC = False, 
        dmvar = "0.5*mZ1*TMath::Hypot(LepGood_ptErrTk[0]/LepGood_pt[0], LepGood_ptErrTk[1]/LepGood_pt[1])",
        dmbins = 25, dmrange=[0.5,3.5]):
    hist = makeHist2D(tree,options,weightedMC=weightedMC, dmvar=dmvar, dmbins=dmbins, dmrange=dmrange)
    print "Splitting."
    ret = []; xbins, xmin, xmax = map(float, options.xvar[1].split(",")) 
    for b in xrange(1,dmbins+1):
        hist1D = ROOT.TH1D("hist%d" % b,"hist",int(xbins),xmin,xmax)
        hist1D.SetDirectory(None)
        for bx in xrange(1,hist.GetNbinsX()+1):
            hist1D.SetBinContent(bx, hist.GetBinContent(bx,b))
        hist1D.dm = [ hist.GetYaxis().GetBinLowEdge(b), hist.GetYaxis().GetBinUpEdge(b) ]
        if options.xcut and (hist1D.dm[1] < options.xcut[0] or hist1D.dm[0] > options.xcut[1]):
            continue
        ret.append(hist1D)
    return ret


def processOneDMSlice(args):
    i,hist,dm,refhist,myparams,options = args
    res = {}; refres = {}
    name = options.name+"_dm_%.2f_%.2f" % (dm[0], dm[1])
    (w,data,pdf, params, result) = fit1D(hist, options)
    makePlot1D(w, data, pdf, params, result, name, options, hist)
    for x in myparams:
        var = result.floatParsFinal().find(x)
        res[x] = ( 0.5*(dm[0] + dm[1]), 0.5*(dm[1] - dm[0]), 
                   var.getVal(), var.getError() )
    if w.function("resol"):
        foms = getResolutionFOMs(w)
        for k,v in foms.iteritems(): res[k] = v
    if refhist:
        (wref,dataref, pdfref, params, refresult) = fit1D(refhist, options)
        makePlot1D(wref, dataref, pdfref, params, refresult, name+"_ref", options, refhist)
        makePlot1DRef(w, data, pdf, pdfref, params, result, refresult, name+"_comp", options)
        for x in myparams:
            var = refresult.floatParsFinal().find(x)
            refres[x] = ( 0.5*(dm[0] + dm[1]), 0.5*(dm[1] - dm[0]), 
                          var.getVal(), var.getError() )
        if w.function("resol"):
            foms = getResolutionFOMs(wref)
            for k,v in foms.iteritems(): refres[k] = v
    return (i,res,refres)

def processOnePtEtaBin(args):
    bin,hist,refhist,myparams,options = args
    res = {}; refres = {}
    name = hist.GetName().replace("hist_",options.name+"_")
    (w,data,pdf, params, result) = fit1D(hist, options, etaptbin=bin)
    makePlot1D(w, data, pdf, params, result, name, options, hist, )
    for x in myparams:
        var = result.floatParsFinal().find(x)
        res[x] = ( var.getVal(), var.getError() )
    if refhist:
        (wref,dataref, pdfref, params, refresult) = fit1D(refhist, options, etaptbin=bin)
        makePlot1D(wref, dataref, pdfref, params, refresult, name+"_ref", options, refhist)
        makePlot1DRef(w, data, pdf, pdfref, params, result, refresult, name+"_comp", options)
        for x in myparams:
            var = refresult.floatParsFinal().find(x)
            refres[x] = ( var.getVal(), var.getError() )
    return (bin,res,refres)

def makeHist3D_(tree, y, z, options, weightedMC = False):
    print "Filling tree for cut %s" % options.cut
    ## cut maker    
    finalCut =  makeCut(options.cut)
    flag = 0  # need to pass 1 for positive charge, 0 for negative, because function triggerSFforChargedLeptonMatchingTriggerWlike takes as first argument a bool, which is 1 (0) for positive (negative) charge
    # the flag is typically defined using function isOddEvent(evt), which will be true (false) for odd (even) events, for which we generally want charge plus (minus)
    if options.useAllEvents:
        if options.selectCharge == "plus":
            flag = 1
        elif options.selectCharge == "minus":
            flag = 0

    if weightedMC: 
        funcTrigSF = ""
        if options.useAllEvents:
            funcTrigSF = "triggerSFforChargedLeptonMatchingTriggerWlike({f},(LepGood_matchedTrgObjMuDR[0]>0.0 || LepGood_matchedTrgObjTkMuDR[0] > 0.0),(LepGood_matchedTrgObjMuDR[1] >0.0 || LepGood_matchedTrgObjTkMuDR[1]>0.0),LepGood_pdgId[0],LepGood_pdgId[1],{pt}[0],{pt}[1],LepGood_eta[0],LepGood_eta[1])".format(f=flag,pt=options.ptvar)
        else:
            funcTrigSF = "triggerSFforChargedLeptonMatchingTriggerWlike(isOddEvent(evt),(LepGood_matchedTrgObjMuDR[0]>0.0 || LepGood_matchedTrgObjTkMuDR[0] > 0.0),(LepGood_matchedTrgObjMuDR[1] >0.0 || LepGood_matchedTrgObjTkMuDR[1]>0.0),LepGood_pdgId[0],LepGood_pdgId[1],{pt}[0],{pt}[1],LepGood_eta[0],LepGood_eta[1])".format(pt=options.ptvar)

        # MC weight = 36300*1999.0*2346.27/84696508770.880005        
        # i.e. luminosity*xsec*genWeight/SumWeights
        # sumWeights now is hardcoded, should be taken summing integral of histogram from each root file
        # 2346.27 is the genWeight for pilot DY sample with photos, but it has a sign
        # so it is 36300*1999.0/84696508770.880005 = 0.000856749599872
        wgt = "genWeight * 0.000856749599872 * puw2016_nTrueInt_36fb(nTrueInt)*_get_muonSF_recoToSelection(LepGood_pdgId[0],{pt}[0],LepGood_eta[1])*_get_muonSF_recoToSelection(LepGood_pdgId[1],{pt}[1],LepGood_eta[1])* {trigSF} * jetPrefireSF_JetAll_Clean".format(pt=options.ptvar, trigSF=funcTrigSF)
        finalCut = "(%s) * %s" % (finalCut, wgt)
    # trigger match and event parity if needed
    if options.selectCharge in ["plus", "minus"]:
        if options.useAllEvents:            
            finalCut += " * triggerMatchWlike({f},LepGood_matchedTrgObjMuDR[0],LepGood_matchedTrgObjMuDR[1],LepGood_matchedTrgObjTkMuDR[0],LepGood_matchedTrgObjTkMuDR[1],LepGood_pdgId[0],LepGood_pdgId[1])".format(f=flag)
        else:
            finalCut += " * triggerMatchWlike(isOddEvent(evt),LepGood_matchedTrgObjMuDR[0],LepGood_matchedTrgObjMuDR[1],LepGood_matchedTrgObjTkMuDR[0],LepGood_matchedTrgObjTkMuDR[1],LepGood_pdgId[0],LepGood_pdgId[1]) * {func}(evt)".format(func="isOddEvent" if options.selectCharge=="plus" else "isEvenEvent") 
    ## end of cut maker    
    
    #finalCut += " * ({pt}[0]<{ptmax} && {pt}[1]<{ptmax} )".format(pt=options.ptvar,ptmax=options.ptbins.split(',')[-1])
    finalCut = scalarToVector(finalCut)
    yexpr, yedges = y
    zexpr, zedges = z
    xbins, xmin, xmax = map(float, options.xvar[1].split(","))
    xedges = [ xmin + i * (xmax-xmin)/int(xbins) for i in xrange(0,int(xbins)+1) ]
    hist = ROOT.TH3D("hist","hist", len(xedges)-1,array('f',xedges),  len(yedges)-1,array('f',yedges),  len(zedges)-1,array('f',zedges))
    nent = tree.Draw("%s:%s:%s>>hist" % (scalarToVector(zexpr),scalarToVector(yexpr),scalarToVector(options.xvar[0])), finalCut, "goff", options.maxEntries)
    if nent < 0:
        print "Warning: something fishy happenend with TTree::Draw() in makeHist3D_(). Abort"
        quit()
    hist = ROOT.gROOT.FindObject("hist")
    if weightedMC:
        for bx in xrange(hist.GetNbinsX()+2): 
          for by in xrange(hist.GetNbinsY()+2): 
             for bz in xrange(hist.GetNbinsZ()+2): 
                hist.SetBinContent(bx,by,bz, max(0, hist.GetBinContent(bx,by,bz)))
    #if options.saveHistoInFile:
    # always save histogram, why not, unless it was loaded from file
    if not options.loadHistoFromFile:
        histToSave = hist.Clone("%s_histo_mass_pt_eta" % ("mc" if weightedMC else "data"))
        outf = ROOT.TFile.Open(options.printDir + "/histo3D_mass_pt_eta.root","UPDATE")
        outf.cd()
        histToSave.Write(histToSave.GetName(),ROOT.TObject.kOverwrite)
        outf.Close()

    return hist

def makeHistsMPtEta(tree, ptedges, etaedges, options,  weightedMC = False):
    hl1 = None
    # get or make histogram
    if options.loadHistoFromFile:
        fname = options.loadHistoFromFile
        hname = "%s_histo_mass_pt_eta" % ("mc" if weightedMC else "data")
        infile = ROOT.TFile.Open(fname,"READ")
        hl1 = infile.Get(hname)
        if not hl1:
            print "Warning: I could not get histogram %s from file %s. Please check" % (hname,fname)
            quit()
        else:
            hl1.SetDirectory(0)
        infile.Close()
    else:   
        useSignedEta = False
        if float(options.etabins.split(',')[0]) < 0.0:
            useSignedEta = True
        if options.selectCharge in ["plus", "minus"]:
            if options.useAllEvents:
                selCharge = 0;
                if options.selectCharge == "plus": 
                    selCharge = 1
                else:
                    selCharge = -1
                funcPt = "returnChargeValAllEvt({c},{pt}[0],LepGood_charge[0],{pt}[1],LepGood_charge[1],0.0)".format(c=selCharge,pt=options.ptvar)
                funcEta = "returnChargeValAllEvt({c},LepGood_eta[0],LepGood_charge[0],LepGood_eta[1],LepGood_charge[1],-999.0)".format(c=selCharge)
            else:
                funcPt = "returnChargeVal({pt}[0],LepGood_charge[0],{pt}[1],LepGood_charge[1],evt)".format(pt=options.ptvar)
                funcEta = "returnChargeVal(LepGood_eta[0],LepGood_charge[0],LepGood_eta[1],LepGood_charge[1],evt)"

            funcEta = "{feta}".format(feta=funcEta) if useSignedEta else "abs({feta})".format(feta=funcEta)
            hl1 = makeHist3D_(tree, (funcPt, ptedges), ("{feta}".format(feta=funcEta), etaedges), options, weightedMC=weightedMC)
        else:
            funcEta = "LepGood_eta[0]" if useSignedEta else "abs(LepGood_eta[0])"
            hl1 = makeHist3D_(tree, ("{pt}[0]".format(pt=options.ptvar), ptedges), ("{feta}".format(feta=funcEta), etaedges), options, weightedMC=weightedMC)
            hl1 = hl1.Clone("hist_l1");
            funcEta = "LepGood_eta[1]" if useSignedEta else "abs(LepGood_eta[1])"
            hl2 = makeHist3D_(tree, ("{pt}[1]".format(pt=options.ptvar), ptedges), ("{feta}".format(feta=funcEta), etaedges), options, weightedMC=weightedMC)
            hl2 = hl2.Clone("hist_l2");
            hl1.Scale(0.5);
            hl2.Scale(0.5);
            hl1.Add(hl2)
    ##
    # now we have the histogram

    xbins, xmin, xmax = map(float, options.xvar[1].split(","))
    ret = []
    for ipt in xrange(1,len(ptedges)):
        for ieta in xrange(1,len(etaedges)):
            name = "hist_pt_%.1f_%.1f_eta_%.3f_%.3f" % (ptedges[ipt-1],ptedges[ipt],etaedges[ieta-1],etaedges[ieta])
            name = name.replace(".","p")
            name = name.replace("-","m")
            hist1D = ROOT.TH1D(name,name, int(xbins),xmin,xmax)
            hist1D.SetDirectory(None)
            for bx in xrange(1,hist1D.GetNbinsX()+1):
                hist1D.SetBinContent(bx, hl1.GetBinContent(bx,ipt,ieta))
            hist1D.bin = (ipt,ieta)
            ret.append(hist1D)
    h2d = hl1.Project3D("yz")
    if not h2d:
        print "Error in makeHistsMPtEta(): h2d is Null"
        quit()
    return (h2d,ret)
            
def styleScatterData(gdata,gmc=None):
        gdata.SetMarkerStyle(ROOT.kFullCircle)
        gdata.SetMarkerColor(ROOT.kGreen+2)
        gdata.SetLineColor(ROOT.kGreen+2)
def styleScatterMC(gmc):
        gmc.SetMarkerStyle(ROOT.kOpenCircle)
        gmc.SetMarkerColor(ROOT.kGray+3)
        gmc.SetLineColor(ROOT.kGray+3)
 
def addZFitterOptions(parser):
    parser.add_option("-n", "--name",   dest="name", default='plot', help="name");
    parser.add_option("-r", "--refmc",   dest="refmc", default=None, help="path to mc samples for chain");
    parser.add_option("-d", "--data",   dest="data", default=None, help="path to data samples for chain");
    parser.add_option("--data-regexp",   dest="dataRegexp", type="string", default="SingleMuon.*", help="regular expression to use some data samples")
    parser.add_option("--mc-regexp",   dest="mcRegexp", type="string", default="ZJToMuMu_powhegMiNNLO_pythia8_photos.*", help="regular expression to use some MC samples")
    parser.add_option("-m", "--mode",   dest="mode", default='1D_PtEtaSlices', help="mode");
    parser.add_option("-s", "--signalModel",   dest="signalModel", default='Z-CB', help="Signal model");
    parser.add_option("-b", "--backgroundModel",   dest="backgroundModel", default='Expo', help="Background model");
    parser.add_option("-t", "--tree",    dest="tree", default='tree', help="Tree name");
    parser.add_option("-f", "--friends", dest="friends", default='Friends', help="Friend tree name (pass empty string to disable friends");
    parser.add_option("-c", "--cut",     dest="cut", type="string", default="Zee", help="cut")
    parser.add_option("--xcut",     dest="xcut", type="float", nargs=2, default=None, help="x axis cut")
    parser.add_option("-x", "--x-var",   dest="xvar", type="string", default=(" mass_2(LepGood1_rocPt,LepGood1_eta,LepGood1_phi,0.1057,LepGood2_rocPt,LepGood2_eta,LepGood2_phi,0.1057)","80,70,110"), nargs=2, help="X var and bin")
    parser.add_option("--ptvar",   dest="ptvar", type="string", default="LepGood_pt", help="pt variable in tree to use (default is uncorrected one)")
    parser.add_option("--fit-range",    dest="fitRange", type="float", nargs=2, default=None, help="Range for fit on variable specified with -x")
    parser.add_option("--setRangeClosure",   dest="setRangeClosure", type="float", default=0.0, help="Range for final closure 2D plot in GeV. Pass positive value X, range is set to [-X,X]. By default, use histogram range extended by 20%")
    parser.add_option("--fit-strategy",  dest="fitStrategy", type="int", default="1", help="Roofit fit strategy passed to RooFit::Strategy(XX): XX can be 0,1,2 where 2 is the slowest but most accurate strategy")
    parser.add_option("--select-charge",   dest="selectCharge", type="string", default="all", help="Use only one charge (plus|minus) or all (all): in latter case both leptons are used to fill the distribution")
    parser.add_option("--ptbins",   dest="ptbins", type="string", default="23,30,35,40,45,50,55,60", help="Comma separated list of pt bin edges")
    parser.add_option("--etabins",   dest="etabins", type="string", default="-2.4,-1.6,-0.8,0,0.8,1.6,2.4", help="Comma separated list of eta bin edges")
    parser.add_option("--useAllEvents",    dest="useAllEvents", default=False, action='store_true' , help="When selecting charge, use all events instead of using only event with given parity to select the charged lepton (this will imply statistical correlations among the two charges)");
    parser.add_option("-u", "--no-weigth-mc",   dest="noWeightMC", default=False, action='store_true' , help="Do not use weigths for MC (scale factors and PU weight)");
    #parser.add_option("--chi2fit",   dest="chi2fit", default=False, action='store_true' , help="Use Chi^2 fit (default is maximum likelihood fit)");
    parser.add_option(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (let it be odd, so the central strip is white if not using --abs-value and the range is symmetric)')
    parser.add_option(     '--palette'  , dest='palette',      default=0, type=int, help='Set palette: default is a built-in one, 1 is colore one , 55 is kRainbow')
    parser.add_option("--xtitle",   dest="xtitle", type="string", default="mass (GeV)", help="X title")
    parser.add_option("--textSize",   dest="textSize", type="float", default=0.04, help="Text size")
    parser.add_option("-l","--lumi",   dest="lumi", type="float", default=35.9, help="Text size")
    parser.add_option("--pdir", "--print-dir", dest="printDir", type="string", default="plots", help="print out plots in this directory");
    parser.add_option("--max-entries",     dest="maxEntries", default=1000000000, type="int", help="Max entries to process in each tree") 
    parser.add_option("-j", "--jobs",    dest="jobs",      type="int",    default=0, help="Use N threads");
    parser.add_option("--rebin",    dest="rebin",      type="int",    default=1, help="Rebin mass plots by N");
    #parser.add_option("--saveHistoInFile",   dest="saveHistoInFile", default=False, action='store_true', help="Save histogram in file for later usage (can be used as input with --histoFromFile to avoid rerunning on ntuples to change fits)")
    parser.add_option("--loadHistoFromFile",   dest="loadHistoFromFile", type="string", default="", help="Load histogram from file instead of running on ntuples. Pass file name")
    parser.add_option("--plot-extension",   dest="plotExtension", type="string", default="png,pdf", help="Comma separated list of extensions to save plots (can just use png to save space)")
    parser.add_option("--roofitPrintLevel",    dest="roofitPrintLevel",      type="int",    default=0, help="RooFit::PrintLevel flag in fitTo: can use -1 to suppress everything, 0 for minimal output, 1 is default for RooFit");

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] tree [tree2] ")
    addZFitterOptions(parser)
    (options, args) = parser.parse_args()
    if not options.data:
        print "You must specify at least one data tree to fit"
        quit()
    ROOT.gROOT.SetBatch(True)
    if not os.path.exists(options.printDir):
        os.system("mkdir -p "+options.printDir)
        if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/m/mciprian/public/index.php "+options.printDir)


    #ROOT.TH1.SetDefaultSumw2()  ## keep commented, otherwise for some mysterious reason the mass distribution for MC has all the data points with error equal to the content. Anyway, the histograms are still filled with weigths, so the errors should be correct (and when fitting one only needs to set ROOT.RooFit.SumW2Error() to use the actual MC stat to compute the uncertainty)
    
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    # built-in palette from light blue to orange
    ROOT.TColor.CreateGradientColorTable(3,
                                         array ("d", [0.00, 0.50, 1.00]),
                                         ##array ("d", [1.00, 1.00, 0.00]),
                                         ##array ("d", [0.70, 1.00, 0.34]),
                                         ##array ("d", [0.00, 1.00, 0.82]),
                                         array ("d", [0.00, 1.00, 1.00]),
                                         array ("d", [0.34, 1.00, 0.65]),
                                         array ("d", [0.82, 1.00, 0.00]),
                                         255,  0.95)



    # tree = None
    # if "*" or "?" in options.data:
    #     tree = ROOT.TChain(options.tree)
    #     files = glob.glob(options.data)
    #     for fname in files: 
    #         print fname
    #         tree.Add(fname)
    # else:
    #     datafile = ROOT.TFile.Open(options.data)
    #     tree = datafile.Get(options.tree)

    # if options.refmc:
    #     if "*" or "?" in options.refmc:
    #         reftree = ROOT.TChain(options.tree)
    #         files = glob.glob(options.refmc)
    #         for fname in files: reftree.Add(fname)
    #     else:
    #         reffile = ROOT.TFile.Open(options.refmc)
    #         reftree = reffile.Get(options.tree)
    # some hardcoded stuff for tests

    ## prepare trees
    tree = None
    reftree = None
    if not options.loadHistoFromFile:
        print "Preparing trees to run on ..."
        # data
        dataMatch = re.compile(options.dataRegexp)
        datapath = options.data
        if not datapath.endswith("/"): datapath += "/"
        folders = [x for x in os.listdir(datapath) if dataMatch.match(x)]
        tree = ROOT.TChain(options.tree)
        ftree = ROOT.TChain(options.friends)
        for i,f in enumerate(folders):
            tree.Add(datapath + f + "/treeProducerWMass/tree.root")
            ftree.Add(datapath + "friends/tree_Friend_" + f + ".root")
        if tree.GetEntries() != ftree.GetEntries():
            print "Error: data trees and friends have different number of entries"
            quit()
        print "Data chain made of %d trees: number of events (unweighted) = %d" % (tree.GetNtrees(),tree.GetEntries())
        tree.AddFriend(ftree)
        # and now mc
        mcMatch = re.compile(options.mcRegexp)
        mcpath = options.refmc
        if not mcpath.endswith("/"): mcpath += "/"
        folders = [x for x in os.listdir(mcpath) if mcMatch.match(x)]
        reftree = ROOT.TChain(options.tree)
        refftree = ROOT.TChain(options.friends)
        for i,f in enumerate(folders):
            reftree.Add(mcpath + f + "/treeProducerWMass/tree.root")
            refftree.Add(mcpath + "friends/tree_Friend_" + f + ".root")
        if reftree.GetEntries() != refftree.GetEntries():
            print "Error: mc trees and friends have different number of entries"
            quit()
        print "MC chain made of %d trees: number of events (unweighted) = %d" % (reftree.GetNtrees(),reftree.GetEntries())
        reftree.AddFriend(refftree)
        #### Done with trees

    useWeightMC = False if options.noWeightMC else True

    if options.mode == "1D":
        hist =   makeHist1D(tree, options)
        (w,data,pdf, params, result) = fit1D(hist, options)
        makePlot1D(w, data, pdf, params, result, options.name, options,hist)
        if options.refmc:
            refhist = makeHist1D(reftree, options, weightedMC=useWeightMC)
            (wref,dataref, pdfref, params, refresult) = fit1D(refhist, options)
            makePlot1D(wref, dataref, pdfref, params, refresult, options.name+"_ref", options, refhist)
            makePlot1DRef(w, data, pdf, pdfref, params, result, refresult, options.name+"_comp", options)
    elif options.mode == "2D":
        hist =   makeHist2D(tree, options)
        (w,data,pdf, params, result) = fit2D(hist, options)
        makePlot2D(w, data, pdf, params, result, options.name, options)
        if options.refmc:
            refhist = makeHist2D(reftree, options, weightedMC=useWeightMC)
            (wref,dataref, pdfref, params, refresult) = fit2D(refhist, options)
            makePlot2D(wref, dataref, pdfref, params, refresult, options.name+"_ref", options)
        #    makePlot1DRef(w, data, pdf, pdfref, params, result, refresult, options.name+"_comp", options)
    elif options.mode == "1D_DMSlices":
        hists =   makeHist1D_DMSlices(tree, options)
        if options.refmc:
            refhists = makeHist1D_DMSlices(reftree, options, weightedMC=useWeightMC)
        else: 
            refhists = [None for h in hists]
        myparams = [ "dm", "sigma" ]
        tasks = [ (i,hist,hist.dm,refhist,myparams,options) for i,(hist, refhist) in enumerate(zip(hists,refhists)) ]
        if options.jobs > 0:
            from multiprocessing import Pool
            fits = Pool(options.jobs).map(processOneDMSlice, tasks)
        else:
            fits = map(processOneDMSlice, tasks)
        # strip clearly bad fits
        fits = [ (i,res,refres) for (i,res,refres) in fits if abs(res["dm"][2]) < 3 ]
        if len(refres):
            fits = [ (i,res,refres) for (i,res,refres) in fits if abs(refres["dm"][2]) < 3 ]
        gdata = dict([ (x,ROOT.TGraphErrors(len(fits))) for x in myparams ])
        gref  = dict([ (x,ROOT.TGraphErrors(len(fits))) for x in myparams ])
        for (i, res, refres) in fits:
           for x in myparams:
                gdata[x].SetPoint(i,      res[x][0], res[x][2] )
                gdata[x].SetPointError(i, res[x][1], res[x][3] )
           if len(refres):
               for x in myparams:
                    gref[x].SetPoint(i,      refres[x][0], refres[x][2] )
                    gref[x].SetPointError(i, refres[x][1], refres[x][3] )
        c1 = ROOT.TCanvas("c1","c1")
        xmax = max(hist.dm[1] for hist in hists) 
        framedm = ROOT.TH2D("frame",";#sigma(pred);#sigma(fit)",100,0,xmax*1.3,100,0,xmax*1.3)
        framedm.Draw(); ROOT.gStyle.SetOptStat(0)
        line = ROOT.TLine(0,0,xmax*1.3,xmax*1.3)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.DrawLine(0,0,xmax*1.3,xmax*1.3)
        line.SetLineColor(ROOT.kRed+1)
        line.DrawLine(0,0,xmax*1.3/1.2,xmax*1.3    )
        line.DrawLine(0,0,xmax*1.3    ,xmax*1.3/1.2)
        line.SetLineWidth(1)
        line.SetLineColor(ROOT.kOrange+7)
        line.DrawLine(0,0,xmax*1.3/1.1,xmax*1.3    )
        line.DrawLine(0,0,xmax*1.3    ,xmax*1.3/1.1)
        if refhists[0] != None:
            styleScatterMC(gref["sigma"])
            gref["sigma"].Draw("P SAME")
        styleScatterData(gdata["sigma"])
        gdata["sigma"].Draw("P SAME")
        printCanvas(c1, options.name+"_summary", [], options)
        if "sigmaEff" in fits[0][1]:
            framedm.GetYaxis().SetTitle("#sigma_{eff}(fit)")
            for (i, res, refres) in fits:
                gdata["sigma"].SetPoint(i,      res["sigma"][0], res["sigmaEff"] )
                gdata["sigma"].SetPointError(i, res["sigma"][1], res["sigmaEff"]*res["sigma"][3]/res["sigma"][2] )
                if refhists[0] != None:
                    gref["sigma"].SetPoint(i,      refres["sigma"][0], refres["sigmaEff"] )
                    gref["sigma"].SetPointError(i, refres["sigma"][1], refres["sigmaEff"]*refres["sigma"][3]/refres["sigma"][2] )
            framedm.Draw(); ROOT.gStyle.SetOptStat(0)
            line = ROOT.TLine(0,0,xmax*1.3,xmax*1.3)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.DrawLine(0,0,xmax*1.3,xmax*1.3)
            line.SetLineColor(ROOT.kRed+1)
            line.DrawLine(0,0,xmax*1.3/1.2,xmax*1.3    )
            line.DrawLine(0,0,xmax*1.3    ,xmax*1.3/1.2)
            line.SetLineWidth(1)
            line.SetLineColor(ROOT.kOrange+7)
            line.DrawLine(0,0,xmax*1.3/1.1,xmax*1.3    )
            line.DrawLine(0,0,xmax*1.3    ,xmax*1.3/1.1)
            if refhists[0] != None:
                gref["sigma"].Draw("P SAME")
            gdata["sigma"].Draw("P SAME")
            printCanvas(c1, options.name+"_summary_eff", [], options)
    elif options.mode == "1D_PtEtaSlices":
        #ptbins = [23,30,35,40,45,50,60,70]; # depends on skim in the ntuples
        #etabins = [0, 1.0, 1.5, 2.1, 2.5] if "Zee" in options.cut else [0, 0.8, 1.6, 2.4]
        ptbins = [float(x) for x in options.ptbins.split(',')]
        etabins = [float(x) for x in options.etabins.split(',')]
        frame2D, hists = makeHistsMPtEta(tree, ptbins, etabins, options)
        if options.refmc:
            _, refhists = makeHistsMPtEta(reftree, ptbins, etabins, options, weightedMC=useWeightMC)
        else: 
            refhists = [None for h in hists]
        myparams = [ "dm", "sigma" ]
        tasks = [ (hist.bin,hist,refhist,myparams,options) for (hist, refhist) in zip(hists,refhists) ]
        if options.jobs > 0:
            from multiprocessing import Pool
            fits = Pool(options.jobs).map(processOnePtEtaBin, tasks)
        else:
            fits = map(processOnePtEtaBin, tasks)
        # strip clearly bad fits
        fits = [ (i,res,refres) for (i,res,refres) in fits if abs(res["dm"][0]) < 3 ]
        if options.refmc:
            fits = [ (i,res,refres) for (i,res,refres) in fits if abs(refres["dm"][0]) < 3 ]
        # get bin labels from hist names
        binLabel = {}; binIndex = {}
        for i,hist in enumerate(hists):
            ptmin,ptmax,etamin,etamax = ptbins[hist.bin[0]-1],ptbins[hist.bin[0]],etabins[hist.bin[1]-1],etabins[hist.bin[1]]
            binIndex[hist.bin] = i
            #binLabel[hist.bin] = "#splitline{%g < p_{T} < %g}{%g < |#eta| < %g}" % (ptmin,ptmax,etamin,etamax)
            binLabel[hist.bin] = "%g < p_{T} < %g   %g < #eta < %g" % (ptmin,ptmax,etamin,etamax)
        frame2D.GetXaxis().SetNdivisions(505)
        frame2D.SetMarkerSize(2.0)
        # make graphs
        gdata = dict([ (x,ROOT.TGraphErrors(len(fits))) for x in myparams ])
        gref  = dict([ (x,ROOT.TGraphErrors(len(fits))) for x in myparams ])
        gdiff = dict([ (x,ROOT.TGraphErrors(len(fits))) for x in myparams ])
        # make 2Ds
        hdata = dict([ (x,frame2D.Clone(options.name+"_"+x        )) for x in myparams ])
        href  = dict([ (x,frame2D.Clone(options.name+"_"+x+"_ref" )) for x in myparams ])
        hdiff = dict([ (x,frame2D.Clone(options.name+"_"+x+"_diff")) for x in myparams ])
        text = dict([ (x,[]) for x in myparams ])
        for i,(bin, res, refres) in enumerate(fits):
           header = hists[binIndex[bin]].GetName().replace("hist_","").replace("_"," ")
           for x in myparams:
                gdata[x].SetPoint(     i, res[x][0], binIndex[bin]+0.6 )
                gdata[x].SetPointError(i, res[x][1],  0                )
                hdata[x].SetBinContent(bin[1],bin[0], res[x][0])
                hdata[x].SetBinError(  bin[1],bin[0], res[x][1])
           if options.refmc:
               for x in myparams:
                    gref[x].SetPoint(     i, refres[x][0], binIndex[bin]+0.4 )
                    gref[x].SetPointError(i, refres[x][1],  0                )
                    href[x].SetBinContent(bin[1],bin[0], refres[x][0] )
                    href[x].SetBinError(  bin[1],bin[0], refres[x][1] )
                    diff = res[x][0] - refres[x][0] if x != "sigma" else res[x][0]/refres[x][0]
                    err  = hypot(res[x][1],refres[x][1]) if x != "sigma" else diff*hypot(res[x][1]/res[x][0],refres[x][1]/refres[x][0])
                    gdiff[x].SetPoint(i,      diff, binIndex[bin]+0.5)
                    gdiff[x].SetPointError(i, err , 0.0)
                    hdiff[x].SetBinContent(bin[1],bin[0], diff )
                    hdiff[x].SetBinError(  bin[1],bin[0], err  )
                    text[x].append("%s    %+5.3f +- %5.3f     %+5.3f +- %5.3f     %+5.3f +- %5.3f" % ( 
                                        header, res[x][0], res[x][1], refres[x][0], refres[x][1], diff, err))
           else:
              for x in myparams:
                   text[x].append("%s    %+5.3f +- %5.3f " % (header, res[x][0], res[x][1]))
        c1 = ROOT.TCanvas("c1","c1",1800,1200)
        c1.SetTickx(1)
        c1.SetTicky(1)
        #c1.SetGrid()
        c1.SetRightMargin(0.2)
        ROOT.gStyle.SetOptStat(0)
        if options.palette > 0: ROOT.gStyle.SetPalette(options.palette);
        ROOT.gStyle.SetNumberContours(options.nContours)
        ROOT.gStyle.SetPaintTextFormat(".3f")
        def normalizeZ(h,pivot=0):
            swing = max(abs(h.GetMaximum()-pivot), abs(h.GetMinimum()-pivot))
            h.GetZaxis().SetRangeUser(pivot-1.5*swing,pivot+1.5*swing)
        for x in myparams:
            if x == "dm": normalizeZ(hdata[x])
            #hdata[x].SetContour(100)
            hdata[x].Draw("COLZ0 TEXTE")
            setAxisTH2(hdata[x],options)
            hdata[x].GetZaxis().SetTitle({'dm':"#Deltam  (GeV)",'sigma':'#sigma(m)  (GeV)'}[x])
            printCanvas(c1, hdata[x].GetName(), text[x], options, xoffs=-0.1)
            if options.refmc:
                if x == "dm": normalizeZ(href[x])
                #href[x].SetContour(100)
                setAxisTH2(href[x],options)
                href[x].Draw("COLZ0 TEXTE")
                href[x].GetZaxis().SetTitle({'dm':"#Deltam  (GeV)",'sigma':'#sigma(m)  (GeV)'}[x])
                printCanvas(c1, href[x].GetName(), text[x], options, xoffs=-0.1)
                #normalizeZ(hdiff[x], pivot = 1 if x == "sigma" else 0)
                maxz = float(max(abs(hdiff[x].GetBinContent(hdiff[x].GetMinimumBin())),abs(hdiff[x].GetBinContent(hdiff[x].GetMaximumBin()))))
                if options.setRangeClosure and x == "dm":
                    hdiff[x].GetZaxis().SetRangeUser(-1.0*options.setRangeClosure,options.setRangeClosure)
                else:
                    hdiff[x].GetZaxis().SetRangeUser(-1.2*maxz,1.2*maxz)
                #hdiff[x].SetContour(100)
                ROOT.gStyle.SetPaintTextFormat(".3f")
                setAxisTH2(hdiff[x],options)
                hdiff[x].Draw("COLZ0 TEXTE")
                hdiff[x].GetZaxis().SetTitle({'dm'    :"#Deltam - #Deltam_{MC}  (GeV)",
                                              'sigma' :'#sigma/#sigma_{MC}'}[x]
                )
                printCanvas(c1, hdiff[x].GetName(), text[x], options, xoffs=-0.1)
                hdiff[x].SaveAs("%s/%s.root" % (options.printDir, hdiff[x].GetName()))
        c1 = ROOT.TCanvas("c1","c1",1000,1000)
        c1.SetTickx(1)
        c1.SetTicky(1)
        c1.SetGrid()
        c1.SetRightMargin(0.04)
        c1.SetLeftMargin(0.3)
        for x in myparams:
            xmax = max(g[x].GetX()[i] + 1.3*g[x].GetErrorX(i) for g in (gdata,gref) for i in xrange(len(fits)))
            xmin = min(g[x].GetX()[i] - 1.3*g[x].GetErrorX(i) for g in (gdata,gref) for i in xrange(len(fits)))
            dx = 0.1*(xmax-xmin)
            frame = ROOT.TH2D("frame","", 100, xmin-dx, xmax+dx, len(hists), 0., len(hists))
            frame.GetXaxis().SetTitle({'dm':"#Deltam  (GeV)",'sigma':'#sigma(m)  (GeV)'}[x]);
            frame.GetXaxis().SetNdivisions(505)
            for bin,label in binLabel.iteritems():
                frame.GetYaxis().SetBinLabel(binIndex[bin]+1, label)
            frame.Draw()
            if options.refmc:
                styleScatterMC(gref[x])
                gref[x].Draw("PZ SAME")
            styleScatterData(gdata[x])
            gdata[x].Draw("PZ SAME")
            frame.GetXaxis().SetRangeUser(xmin-dx,xmax+dx)
            gdata[x].GetXaxis().SetRangeUser(xmin-dx,xmax+dx)
            gref[x].GetXaxis().SetRangeUser(xmin-dx,xmax+dx)
            printCanvas(c1, options.name+"_"+x+"_summary", text[x], options)
            if options.refmc:
                xmax = max(g[x].GetX()[i] + 1.3*g[x].GetErrorX(i) for g in (gdiff,) for i in xrange(len(fits)))
                xmin = min(g[x].GetX()[i] - 1.3*g[x].GetErrorX(i) for g in (gdiff,) for i in xrange(len(fits)))
                dx = 0.1*(xmax-xmin)
                frame = ROOT.TH2D("frame","", 100, xmin-dx, xmax+dx, len(hists), 0., len(hists))
                frame.GetXaxis().SetTitle({'dm'    :"#Deltam - #Deltam_{MC}  (GeV)",
                                           'sigma' :'#sigma/#sigma_{MC}'} [x]);
                frame.GetXaxis().SetNdivisions(505)
                for bin,label in binLabel.iteritems():
                    frame.GetYaxis().SetBinLabel(binIndex[bin]+1, label)
                frame.Draw()
                styleScatterData(gdiff[x])
                gdiff[x].Draw("PZ SAME")
                if options.setRangeClosure and x == "dm":
                    frame.GetXaxis().SetRangeUser(-1.0*options.setRangeClosure,options.setRangeClosure)
                    gdiff[x].GetXaxis().SetRangeUser(-1.0*options.setRangeClosure,options.setRangeClosure)
                else:
                    frame.GetXaxis().SetRangeUser(xmin-dx,xmax+dx)
                    gdiff[x].GetXaxis().SetRangeUser(xmin-dx,xmax+dx)
                printCanvas(c1, options.name+"_"+x+"_diff_summary", text[x], options)
                #gdiff[x].SaveAs("%s/%s.root" % (options.printDir, gdiff[x].GetName()))

 
        
        

