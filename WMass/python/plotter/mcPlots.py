#!/usr/bin/env python

from mcAnalysis import *
import CMS_lumi as CMS_lumi
import itertools
import shutil

CMS_lumi.writeExtraText = 1

_global_workspaces=[] # avoid crash in 80X, to be investigated

if "/bin2Dto1Dlib_cc.so" not in ROOT.gSystem.GetLibraries():
    compileMacro("ccFiles/bin2Dto1Dlib.cc")

if "/fakeRate_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/fakeRate.cc")

SAFE_COLOR_LIST=[
ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+2, ROOT.kViolet+5, ROOT.kSpring+5, ROOT.kAzure+1, ROOT.kPink+7, ROOT.kOrange+3, ROOT.kBlue+3, ROOT.kMagenta+3, ROOT.kRed+2,
]+list(range(11,40))

def _unTLatex(string):
    return string.replace("#chi","x").replace("#mu","m").replace("#rightarrow","->")

class PlotFile:
    def __init__(self,fileName,options):
        self._options = options
        self._plots = []
        defaults = {}
        infile = open(fileName,'r')
        for line in infile:
            if re.match("\s*#.*", line) or len(line.strip())==0: continue
            while line.strip()[-1] == "\\":
                line = line.strip()[:-1] + infile.next()
            extra = {}
            if ";" in line:
                (line,more) = line.split(";")[:2]
                more = more.replace("\\,",";")
                for setting in [f.strip().replace(";",",") for f in more.split(',')]:
                    if "=" in setting: 
                        # in following line, if setting has more than one '=' an error will occur, e.g., if you have XTitle="muon isolation (#DeltaR=0.4)"
                        # therefore, split only on the first occurrence of '='
                        (key,val) = [f.strip() for f in setting.split("=",1)]  
                        extra[key] = eval(val)
                    else: extra[setting] = True
            line = re.sub("#.*","",line) 
            field = [f.strip().replace(";",":") for f in line.replace("::",";;").replace("\\:",";").split(':')]
            if len(field) == 1 and field[0] == "*":
                if len(self._plots): raise RuntimeError("PlotFile defaults ('*') can be specified only before all plots")
                logging.info("Setting the following defaults for all plots: ")
                for k,v in extra.items():
                    logging.info("\t%s: %r" % (k,v))
                    defaults[k] = v
                continue
            else:
                for k,v in defaults.items():
                    if k not in extra: extra[k] = v
            if len(field) <= 2: continue
            if len(options.plotselect):
                skipMe = True
                for p0 in options.plotselect:
                    for p in p0.split(","):
                        if re.match(p+"$", field[0]): skipMe = False
                if skipMe: continue
            if len(options.plotexclude):
                skipMe = False
                for p0 in options.plotexclude:
                    for p in p0.split(","):
                        if re.match(p+"$", field[0]): skipMe = True
                if skipMe: continue
            if options.globalRebin: extra['rebinFactor'] = options.globalRebin
            self._plots.append(PlotSpec(field[0],field[1],field[2],extra))
    def plots(self):
        return self._plots[:]

def getDataPoissonErrors(h, drawZeroBins=False, drawXbars=False):
    xaxis = h.GetXaxis()
    q=(1-0.6827)/2.;
    points = []
    errors = []
    for i in range(h.GetNbinsX()):
        N = h.GetBinContent(i+1);
        dN = h.GetBinError(i+1);
        if drawZeroBins or N > 0:
            if N > 0 and dN > 0 and abs(dN**2/N-1) > 1e-4: 
                logging.debug("This is not Poisson to begin with! %.2f, %.2f, neff = %.2f, yscale = %.5g" % (N, dN, (N/dN)**2, (dN**2/N)))
                yscale = (dN**2/N)
                N = (N/dN)**2
            else:
                yscale = 1
            x = xaxis.GetBinCenter(i+1);
            points.append( (x,yscale*N) )
            EYlow  = (N-ROOT.ROOT.Math.chisquared_quantile_c(1-q,2*N)/2.) if N > 0 else 0
            EYhigh = ROOT.ROOT.Math.chisquared_quantile_c(q,2*(N+1))/2.-N;
            EXhigh, EXlow = (xaxis.GetBinUpEdge(i+1)-x, x-xaxis.GetBinLowEdge(i+1)) if drawXbars else (0,0)
            errors.append( (EXlow,EXhigh,yscale*EYlow,yscale*EYhigh) )
    ret = ROOT.TGraphAsymmErrors(len(points))
    ret.SetName(h.GetName()+"_graph")
    for i,((x,y),(EXlow,EXhigh,EYlow,EYhigh)) in enumerate(zip(points,errors)):
        ret.SetPoint(i, x, y)
        ret.SetPointError(i, EXlow,EXhigh,EYlow,EYhigh)
    ret.SetLineWidth(h.GetLineWidth())
    ret.SetLineColor(h.GetLineColor())
    ret.SetLineStyle(h.GetLineStyle())
    ret.SetMarkerSize(h.GetMarkerSize())
    ret.SetMarkerColor(h.GetMarkerColor())
    ret.SetMarkerStyle(h.GetMarkerStyle())
    return ret

def PrintHisto(h):
    if not h:
        return
    logging.info("Hist name is %s" % h.GetName())
    c=[]
    if "TH1" in h.ClassName():
        for i in range(h.GetNbinsX()):
            c.append((h.GetBinContent(i+1),h.GetBinError(i+1)))
    elif "TH2" in h.ClassName():
        for i in range(h.GetNbinsX()):
            for j in range(h.GetNbinsY()):
                c.append((h.GetBinContent(i+1,j+1),h.GetBinError(i+1,j+1)))
    else:
        logging.info('Hist is not th1 or th2')
    logging.info(c)

def doSpam(text,x1,y1,x2,y2,align=12,fill=False,textSize=0.033,_noDelete={}):
    cmsprel = ROOT.TPaveText(x1,y1,x2,y2,"NDC");
    cmsprel.SetTextSize(textSize);
    cmsprel.SetFillColor(0);
    cmsprel.SetFillStyle(1001 if fill else 0);
    cmsprel.SetLineStyle(2);
    cmsprel.SetLineColor(0);
    cmsprel.SetTextAlign(align);
    cmsprel.SetTextFont(42);
    cmsprel.AddText(text);
    cmsprel.Draw("same");
    _noDelete[text] = cmsprel; ## so it doesn't get deleted by PyROOT
    return cmsprel


def doTinyCmsPrelim(textLeft="_default_",textRight="_default_",hasExpo=False,textSize=0.033,lumi=None, xoffs=0, options=None, doWide=False):
    if textLeft  == "_default_": textLeft  = options.lspam
    if textRight == "_default_": textRight = options.rspam
    if lumi      == None       : lumi      = options.lumi
    if   lumi > 3.54e+1: lumitext = "%.0f fb^{-1}" % lumi
    elif lumi > 3.54e+0: lumitext = "%.1f fb^{-1}" % lumi
    elif lumi > 3.54e-1: lumitext = "%.2f fb^{-1}" % lumi
    elif lumi > 3.54e-2: lumitext = "%.0f pb^{-1}" % (lumi*1000)
    elif lumi > 3.54e-3: lumitext = "%.1f pb^{-1}" % (lumi*1000)
    else               : lumitext = "%.2f pb^{-1}" % (lumi*1000)
    lumitext = "%.1f fb^{-1}" % lumi
    textLeft = textLeft.replace("%(lumi)",lumitext)
    textRight = textRight.replace("%(lumi)",lumitext)
    if textLeft not in ['', None]:
        doSpam(textLeft, (.30 if hasExpo else 0.07 if doWide else .19)+xoffs, .955, .60+xoffs, .995, align=12, textSize=textSize)
    if textRight not in ['', None]:
        doSpam(textRight,(0.6 if doWide else .68)+xoffs, .955, .99+xoffs, .995, align=32, textSize=textSize)

def reMax(hist,hist2,islog,factorLin=1.3,factorLog=2.0,doWide=False):
    if  hist.ClassName() == 'THStack':
        hist = hist.GetHistogram()  # better to use sum of all components, GetHistogram() returns a fake TH1 used to build the axis (could have no relation with the plot)
        #hist = hist.GetStack().Last() # this is the sum of all components
    #max0 = hist.GetBinContent(hist.GetMaximumBin())*(factorLog if islog else factorLin)
    #max2 = hist2.GetBinContent(hist2.GetMaximumBin())*(factorLog if islog else factorLin)
    max0 = hist.GetMaximum()
    max2 = hist2.GetMaximum()*(factorLog if islog else factorLin)
    # Below, use a protection against cases where uncertainty is much bigger than value (might happen with weird situations or QCD MC)
    if hasattr(hist2,'poissonGraph'):
       for i in range(hist2.poissonGraph.GetN()):
          if (hist2.poissonGraph.GetErrorYhigh(i) > hist2.poissonGraph.GetY()[i]):
              tmpvalue = (factorLog if islog else factorLin) * hist2.poissonGraph.GetY()[i]
          else:
              tmpvalue = (hist2.poissonGraph.GetY()[i] + 1.3*hist2.poissonGraph.GetErrorYhigh(i))
          max2 = max(max2, tmpvalue*(factorLog if islog else factorLin))
    elif "TH1" in hist2.ClassName():
       for b in range(1,hist2.GetNbinsX()+1):
          if (hist2.GetBinError(b) > hist2.GetBinContent(b)):
              tmpvalue = (factorLog if islog else factorLin) * hist2.GetBinContent(b)
          else:
              tmpvalue = hist2.GetBinContent(b) + 1.3*hist2.GetBinError(b)
          max2 = max(max2, tmpvalue*(factorLog if islog else factorLin))
    if max2 > max0 and not "TH2" in hist.ClassName():
        max0 = max2;
        if islog: hist.GetYaxis().SetRangeUser(0.1 if doWide else 0.9, max0)
        else:     hist.GetYaxis().SetRangeUser(0,max0)

def doShadedUncertainty(h):
    xaxis = h.GetXaxis()
    points = []; errors = []
    for i in range(h.GetNbinsX()):
        N = h.GetBinContent(i+1); dN = h.GetBinError(i+1);
        if N == 0 and dN == 0: continue
        x = xaxis.GetBinCenter(i+1);
        points.append( (x,N) )
        EYlow, EYhigh  = dN, min(dN,N);
        EXhigh, EXlow = (xaxis.GetBinUpEdge(i+1)-x, x-xaxis.GetBinLowEdge(i+1))
        errors.append( (EXlow,EXhigh,EYlow,EYhigh) )
    ret = ROOT.TGraphAsymmErrors(len(points))
    ret.SetName(h.GetName()+"_errors")
    for i,((x,y),(EXlow,EXhigh,EYlow,EYhigh)) in enumerate(zip(points,errors)):
        ret.SetPoint(i, x, y)
        ret.SetPointError(i, EXlow,EXhigh,EYlow,EYhigh)
    ret.SetFillStyle(3244);
    ret.SetFillColor(ROOT.kGray+2)
    ret.SetMarkerStyle(0)
    ret.Draw("PE2 SAME")
    return ret

def doDataNorm(pspec,pmap):
    if "data" not in pmap: return None
    total = sum([v.Integral() for k,v in pmap.items() if k != 'data' and not hasattr(v,'summary')])
    sig = pmap["data"].Clone(pspec.name+"_data_norm")
    sig.SetFillStyle(0)
    sig.SetLineColor(1)
    sig.SetLineWidth(3)
    sig.SetLineStyle(2)
    if sig.Integral() > 0:
        sig.Scale(total/sig.Integral())
    sig.Draw("HIST SAME")
    return sig

def doStackSignalNorm(pspec,pmap,individuals,extrascale=1.0,norm=True):
    total = sum([v.Integral() for k,v in pmap.items() if k != 'data' and not hasattr(v,'summary')])
    if options.noStackSig:
        total = sum([v.Integral() for k,v in pmap.items() if not hasattr(v,'summary') and mca.isBackground(k) ])
    if individuals:
        sigs = []
        for sig in [pmap[x] for x in mca.listSignals() if x in pmap and pmap[x].Integral() > 0]:
            sig = sig.Clone(sig.GetName()+"_norm")
            sig.SetFillStyle(0)
            sig.SetLineColor(sig.GetFillColor())
            sig.SetLineWidth(4)
            if norm: sig.Scale(total*extrascale/sig.Integral())
            sig.Draw("HIST SAME")
            sigs.append(sig)
        return sigs
    else:
        sig = None
        if "signal" in pmap: sig = pmap["signal"].Clone(pspec.name+"_signal_norm")
        else: 
            sigs = [pmap[x] for x in mca.listBackgrounds() if x in pmap and pmap[x].Integral() > 0]
            sig = sigs[0].Clone(sigs.GetName()+"_norm")
        sig.SetFillStyle(0)
        sig.SetLineColor(206)
        sig.SetLineWidth(4)
        if norm and sig.Integral() > 0:
            sig.Scale(total*extrascale/sig.Integral())
        sig.Draw("HIST SAME")
        return [sig]

def doStackSigScaledNormData(pspec,pmap):
    if "data"       not in pmap: return (None,-1.0)
    if "signal"     not in pmap: return (None,-1.0)
    data = pmap["data"]
    sig = pmap["signal"].Clone("sig_refloat")
    bkg = None
    if "background" in pmap:
        bkg = pmap["background"]
    else:
        bkg = sig.Clone(); bkg.Reset()
    sf = (data.Integral()-bkg.Integral())/sig.Integral()
    sig.Scale(sf)
    sig.Add(bkg)
    sig.SetFillStyle(0)
    sig.SetLineColor(206)
    sig.SetLineWidth(3)
    sig.SetLineStyle(2)
    sig.Draw("HIST SAME")
    return (sig,sf)

def doScaleSigNormData(pspec,pmap,mca):
    if "data"       not in pmap: return -1.0
    if "signal"     not in pmap: return -1.0
    data = pmap["data"]
    sig = pmap["signal"].Clone("sig_refloat")
    bkg = None
    if "background" in pmap:
        bkg = pmap["background"]
    else:
        bkg = sig.Clone(); bkg.Reset()
    sf = (data.Integral()-bkg.Integral())/sig.Integral()
    signals = [ "signal" ] + mca.listSignals()
    for p,h in pmap.items():
        if p in signals: h.Scale(sf)
    pspec.setLog("ScaleSig", [ "Signal processes scaled by %g" % sf ] )
    return sf

def doScaleBkgNormData(pspec,pmap,mca,list = []):
    if "data"       not in pmap: return -1.0
    if "background" not in pmap: return -1.0
    if any([l not in pmap for l in list]): return -1.0
    data = pmap["data"]
    bkg  = pmap["background"]
    int = sum([pmap[l].Integral() for l in list])
    rm = bkg.Integral() - int
    sf = (data.Integral() - rm) / int
    bkgs = ["background"] + list
    for p,h in pmap.items():
        if p in bkgs: h.Scale(sf)
    return sf


def doNormFit(pspec,pmap,mca,saveScales=False):
    global _global_workspaces
    if "data" not in pmap: return -1.0
    gKill = ROOT.RooMsgService.instance().globalKillBelow()
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    data = pmap["data"]
    w = ROOT.RooWorkspace("w","w")
    _global_workspaces.append(w)
    x = w.factory("x[%g,%g]" % (data.GetXaxis().GetXmin(), data.GetXaxis().GetXmax()))
    x.setBins(data.GetNbinsX())
    obs = ROOT.RooArgList(w.var("x"))
    hdata = pmap['data']; hdata.killbins = False
    hmc = mergePlots('htemp', [v for (k,v) in pmap.items() if k != 'data'])
    for b in range(1,hmc.GetNbinsX()+2):
        if hdata.GetBinContent(b) > 0 and hmc.GetBinContent(b) == 0:
            if not hdata.killbins:
                hdata = hdata.Clone()
                hdata.killbins = True
            for b2 in range(b-1,0,-1):
                if hmc.GetBinContent(b2) > 0:
                    hdata.SetBinContent(b2, hdata.GetBinContent(b2) + hdata.GetBinContent(b))
                    hdata.SetBinContent(b, 0)
                    break
            if hdata.GetBinContent(b) > 0:
                for b2 in range(b+1,hmc.GetNbinsX()+2):
                    if hmc.GetBinContent(b2) > 0:
                        hdata.SetBinContent(b2, hdata.GetBinContent(b2) + hdata.GetBinContent(b))
                        hdata.SetBinContent(b, 0)
                        break
            if hdata.GetBinContent(b) > 0: hdata.SetBinContent(b, 0)
    rdhs = {};
    w.imp = getattr(w, 'import')
    for p,h in pmap.items():
        rdhs[p] = ROOT.RooDataHist("hist_"+p,"",obs,h if p != "data" else hdata)
        w.imp(rdhs[p])
    pdfs   = ROOT.RooArgList()
    coeffs = ROOT.RooArgList()
    constraints = ROOT.RooArgList()
    dontDelete = []
    procNormMap = {}
    for p in mca.listBackgrounds(allProcs=True) + mca.listSignals(allProcs=True):
        if p not in pmap: continue
        if pmap[p].Integral() == 0: continue
        hpdf = ROOT.RooHistPdf("pdf_"+p,"",ROOT.RooArgSet(x), rdhs[p])
        pdfs.add(hpdf); dontDelete.append(hpdf)
        if mca.getProcessOption(p,'FreeFloat',False):
            normTermName = mca.getProcessOption(p,'PegNormToProcess',p)
            logging.info("%s scale as %s" % (p, normTermName))
            normterm = w.factory('prod::norm_%s(%g,syst_%s[1,%g,%g])' % (p, pmap[p].Integral(), normTermName, 0.2, 5))
            dontDelete.append((normterm,))
            coeffs.add(normterm)
            procNormMap[p] = normterm
        elif mca.getProcessOption(p,'NormSystematic',0.0) > 0:
            syst = mca.getProcessOption(p,'NormSystematic',0.0)
            normTermName = mca.getProcessOption(p,'PegNormToProcess',p)
            logging.info("%s scale as %s with %s constraint" % (p, normTermName, syst))
            normterm = w.factory('expr::norm_%s("%g*pow(%g,@0)",syst_%s[-5,5])' % (p, pmap[p].Integral(), 1+syst, normTermName))
            if not w.pdf("systpdf_%s" % normTermName): 
                constterm = w.factory('Gaussian::systpdf_%s(syst_%s,0,1)' % (normTermName,normTermName))
                constraints.add(constterm)
                dontDelete.append((normterm,constterm))
            else:
                dontDelete.append((normterm))
            coeffs.add(normterm)
            procNormMap[p] = normterm
        else:    
            logging.info("%s is fixed" % p)
            normterm = w.factory('norm_%s[%g]' % (p, pmap[p].Integral()))
            dontDelete.append((normterm,))
            coeffs.add(normterm)
    pdfs.Print("V")
    coeffs.Print("V")
    addpdf = ROOT.RooAddPdf("tot","",pdfs,coeffs,False)
    model  = addpdf
    if constraints.getSize() > 0:
        constraints.add(addpdf)
        model = ROOT.RooProdPdf("prod","",constraints)
    result = model.fitTo( rdhs["data"], ROOT.RooFit.Save(1) )
    totsig, totbkg = None, None
    if "signal" in pmap and "signal" not in mca.listSignals(): 
        totsig = pmap["signal"]; totsig.Reset()
    if "background" in pmap and "background" not in mca.listBackgrounds(): 
        totbkg = pmap["background"]; totbkg.Reset()
    fitlog = []
    for p in mca.listBackgrounds(allProcs=True) + mca.listSignals(allProcs=True):
        normSystematic = mca.getProcessOption(p,'NormSystematic', 0.0)
        if p in pmap and p in procNormMap:
           # setthe scale
           newscale = procNormMap[p].getVal()/pmap[p].Integral()
           pmap[p].Scale(newscale)
           # now get the 1 sigma
           normTermName = mca.getProcessOption(p,'PegNormToProcess',p)
           nuis = w.var("syst_"+normTermName);
           val,err = (nuis.getVal(), nuis.getError())
           v0 =  procNormMap[p].getVal()
           nuis.setVal(val+err)
           v1 =  procNormMap[p].getVal()
           nuis.setVal(val)
           fitlog.append("Process %s scaled by %.3f +/- %.3f" % (p,newscale,newscale*(v1-v0)/v0))
           if saveScales:
              logging.info("Scaling process %s by the extracted scale factor %.3f with rel. syst uncertainty %.3f" % (p,newscale,(v1-v0)/v0))
              mca.setProcessOption(p,'NormSystematic', (v1-v0)/v0);
              mca.scaleUpProcess(p,newscale)
        # recompute totals
        if p in pmap:
            htot = totsig if mca.isSignal(p) else totbkg
            if htot != None:
                htot.Add(pmap[p])
                syst = normSystematic
                if syst > 0:
                    for b in range(1,htot.GetNbinsX()+1):
                        htot.SetBinError(b, hypot(htot.GetBinError(b), pmap[p].GetBinContent(b)*syst))
    pspec.setLog("Fitting", fitlog)
    ROOT.RooMsgService.instance().setGlobalKillBelow(gKill)
    

def doRatioHists(pspec,pmap,total,totalSyst,marange,firange=False,fitRatio=None,errorsOnRef=True,onlyStatErrorsOnRef=False,ratioNums="signal",ratioDen="background",ylabel="Data/pred.",doWide=False,showStatTotLegend=False,errorBarsOnRatio=True,ratioYLabelSize=0.06,ratioNumsWithData=""):
    numkeys =[] if ratioDen == "data" else  [ "data" ]
    if len(ratioNumsWithData): 
        for p in pmap.keys():                
            for s in ratioNumsWithData.split(","):
                #print "p, s : %s,%s" % (p,s)
                # do we want a match or equality? If I have QCD in numerator but I have processes QCD and QCD_1, I will have 2 matches, and this is not what I want
                # if re.match(s,p): 
                if s==p: 
                    numkeys.append(p)
                    break

    if "data" not in pmap or ratioDen == "data": 
        #print str(pmap)
        # >= 3 instead of 4 because I might have no signal process, 
        # while I always have background as sum of everything but data (minimum two processes to make ratio of them)
        if len(pmap) >= 3 and ratioDen in pmap:   
            numkeys = []
            for p in pmap.keys():                
                for s in ratioNums.split(","):
                    #print "p, s : %s,%s" % (p,s)
                    # do we want a match or equality? If I have QCD in numerator but I have processes QCD and QCD_1, I will have 2 matches, and this is not what I want
                    # if re.match(s,p): 
                    if s==p: 
                        numkeys.append(p)
                        break
            if len(numkeys) == 0:
                return (None,None,None,None)
            # do this first
            total.GetXaxis().SetLabelOffset(999) ## send them away
            total.GetXaxis().SetTitleOffset(999) ## in outer space
            total.GetYaxis().SetTitleSize(0.06)
            total.GetYaxis().SetTitleOffset(0.75 if doWide else 1.48)
            if options.setTitleYoffset > 0: total.GetYaxis().SetTitleOffset(options.setTitleYoffset)
            total.GetYaxis().SetLabelSize(0.05)
            total.GetYaxis().SetLabelOffset(0.007)
            # then we can overwrite total with background
            numkey = 'signal'
            total     = pmap[ratioDen]
            totalSyst = pmap[ratioDen]
        else:    
            return (None,None,None,None)
    ratios = [] #None
    for numkey in numkeys:
        if hasattr(pmap[numkey], 'poissonGraph'):
            ratio = pmap[numkey].poissonGraph.Clone(numkey+"_div"); 
            for i in range(ratio.GetN()):
                x    = ratio.GetX()[i]
                div  = total.GetBinContent(total.GetXaxis().FindBin(x))
                ratio.SetPoint(i, x, ratio.GetY()[i]/div if div > 0 else 0)
                ratio.SetPointError(i, ratio.GetErrorXlow(i), ratio.GetErrorXhigh(i), 
                                       ratio.GetErrorYlow(i)/div  if div > 0 else 0, 
                                       ratio.GetErrorYhigh(i)/div if div > 0 else 0) 
        else:
            ratio = pmap[numkey].Clone(numkey+"_div"); 
            ratio.Divide(total)
            if not errorBarsOnRatio:
                for i in range(ratio.GetNbinsX()+2):
                    ratio.SetBinError(i,0)
        ratios.append(ratio)
    unity  = totalSyst.Clone("sim_div");
    unity0 = total.Clone("sim_div");
    rmin, rmax =  1,1
    for b in range(1,unity.GetNbinsX()+1):
        e,e0,n = unity.GetBinError(b), unity0.GetBinError(b), unity.GetBinContent(b)
        unity.SetBinContent(b, 1 if n > 0 else 0)
        unity0.SetBinContent(b,  1 if n > 0 else 0)
        if errorsOnRef:
            if onlyStatErrorsOnRef:
                unity.SetBinError(b, 0)
            else:
                unity.SetBinError(b, e/n if n > 0 else 0)
            unity0.SetBinError(b, e0/n if n > 0 else 0)
        else:
            unity.SetBinError(b, 0)
            unity0.SetBinError(b, 0)
        rmin = min([ rmin, 1-2*e/n if n > 0 else 1])
        rmax = max([ rmax, 1+2*e/n if n > 0 else 1])
    for ratio in ratios:
        if ratio.ClassName() != "TGraphAsymmErrors":
            for b in range(1,unity.GetNbinsX()+1):
                if ratio.GetBinContent(b) == 0: continue
                rmin = min([ rmin, ratio.GetBinContent(b) - 2*ratio.GetBinError(b) ]) 
                rmax = max([ rmax, ratio.GetBinContent(b) + 2*ratio.GetBinError(b) ])  
        else:
            for i in range(ratio.GetN()):
                rmin = min([ rmin, ratio.GetY()[i] - 2*ratio.GetErrorYlow(i)  ]) 
                rmax = max([ rmax, ratio.GetY()[i] + 2*ratio.GetErrorYhigh(i) ])  
    if rmin < marange[0] or firange: rmin = marange[0]; 
    if rmax > marange[1] or firange: rmax = marange[1];
    if (rmax > 3 and rmax <= 3.4): rmax = 3.4
    if (rmax > 2 and rmax <= 2.4): rmax = 2.4
    unity.SetFillStyle(1001);
    unity.SetFillColor(ROOT.kCyan);
    unity.SetMarkerStyle(1);
    unity.SetMarkerColor(ROOT.kCyan);
    unity0.SetFillStyle(1001);
    unity0.SetFillColor(ROOT.kBlue-7);
    unity0.SetMarkerStyle(1);
    unity0.SetMarkerColor(ROOT.kBlue-7);
    ROOT.gStyle.SetErrorX(0.5);
    if errorsOnRef and not onlyStatErrorsOnRef:
        unity.Draw("E2");
    else:
        unity.Draw("AXIS");
    if fitRatio != None and len(ratios) == 1:
        from CMGTools.TTHAnalysis.tools.plotDecorations import fitTGraph
        fitTGraph(ratio,order=fitRatio)
        unity.SetFillStyle(3013);
        unity0.SetFillStyle(3013);
        if errorsOnRef:
            if not onlyStatErrorsOnRef: unity.Draw("AXIS SAME");
            unity0.Draw("E2 SAME");
    else:
        if (total != totalSyst or onlyStatErrorsOnRef) and errorsOnRef:
            unity0.Draw("E2 SAME");
    rmin = float(pspec.getOption("RMin",rmin))
    rmax = float(pspec.getOption("RMax",rmax))
    unity.GetYaxis().SetRangeUser(rmin,rmax);
    unity.GetXaxis().SetTitleFont(42)
    unity.GetXaxis().SetTitleSize(0.14)
    unity.GetXaxis().SetTitleOffset(options.setTitleXoffset)
    unity.GetXaxis().SetLabelFont(42)
    unity.GetXaxis().SetLabelSize(0.1)
    unity.GetXaxis().SetLabelOffset(0.007)
    unity.GetYaxis().SetNdivisions(505)
    unity.GetYaxis().SetTitleFont(42)
    unity.GetYaxis().SetTitleSize(ratioYLabelSize) # 0.14
    offset = 0.32 if doWide else 0.62
    unity.GetYaxis().SetTitleOffset(offset)
    if options.setTitleYoffset > 0: unity.GetYaxis().SetTitleOffset(offset * options.setTitleYoffset / 1.48)
    unity.GetYaxis().SetLabelFont(42)
    unity.GetYaxis().SetLabelSize(0.11)
    unity.GetYaxis().SetLabelOffset(0.007)
    unity.GetYaxis().SetDecimals(True) 
    unity.GetYaxis().SetTitle(ylabel)
    total.GetXaxis().SetLabelOffset(999) ## send them away
    total.GetXaxis().SetTitleOffset(999) ## in outer space
    total.GetYaxis().SetTitleSize(0.065)
    total.GetYaxis().SetTitleOffset(0.75 if doWide else 1.48)
    if options.setTitleYoffset > 0: total.GetYaxis().SetTitleOffset(options.setTitleYoffset)
    total.GetYaxis().SetLabelSize(0.05)
    total.GetYaxis().SetLabelOffset(0.007)
    binlabels = pspec.getOption("xBinLabels","")
    if binlabels != "" and len(binlabels.split(",")) == unity.GetNbinsX():
        blist = binlabels.split(",")
        for i in range(1,unity.GetNbinsX()+1): 
            unity.GetXaxis().SetBinLabel(i,blist[i-1]) 
        unity.GetXaxis().SetLabelSize(0.15)
    #ratio.SetMarkerSize(0.7*ratio.GetMarkerSize()) # no it is confusing
    binlabels = pspec.getOption("xBinLabels","")
    if binlabels != "" and len(binlabels.split(",")) == unity.GetNbinsX():
        blist = binlabels.split(",")
        for i in range(1,unity.GetNbinsX()+1): 
            unity.GetXaxis().SetBinLabel(i,blist[i-1]) 
    #$ROOT.gStyle.SetErrorX(0.0);
    line = ROOT.TLine(unity.GetXaxis().GetXmin(),1,unity.GetXaxis().GetXmax(),1)
    line.SetLineWidth(2);
    line.SetLineColor(58);
    line.Draw("L")
    for ratio in ratios:
        ratio.Draw("E SAME" if ratio.ClassName() != "TGraphAsymmErrors" else "PZ SAME");
        if len(ratioNumsWithData) or ratioDen == "data":
            ratio.SetMarkerColor(ratio.GetLineColor())
            if ratioDen == "data":
                ratio.SetMarkerStyle(0)
                ratio.SetLineColor(ratio.GetFillColor())
                ratio.SetLineWidth(2)
    leg0 = ROOT.TLegend(0.12 if doWide else 0.2, 0.8, 0.25 if doWide else 0.45, 0.9)
    leg0.SetFillColor(0)
    leg0.SetShadowColor(0)
    leg0.SetLineColor(0)
    leg0.SetTextFont(42)
    leg0.SetTextSize(0.035*0.7/0.3)
    leg0.AddEntry(unity0, "stat. bkg. unc.", "F")
    if showStatTotLegend: leg0.Draw()
    leg1 = ROOT.TLegend(0.25 if doWide else 0.45, 0.8, 0.38 if doWide else 0.7, 0.9)
    leg1.SetFillColor(0)
    leg1.SetShadowColor(0)
    leg1.SetLineColor(0)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.035*0.7/0.3)
    leg1.AddEntry(unity, "total bkg. unc.", "F")
    if showStatTotLegend and not onlyStatErrorsOnRef: leg1.Draw()
    global legendratio0_, legendratio1_
    legendratio0_ = leg0
    legendratio1_ = leg1
    return (ratios, unity, unity0, line)

def doRatio2DHists(pspec,pmap,total,totalSyst,marange,firange=False,ratioNums="signal",ratioDen="background",ylabel="Data/pred.",ratioNumsWithData=""):
    numkeys =[] if ratioDen == "data" else  [ "data" ]
    if len(ratioNumsWithData):
        for p in pmap.keys():                
            for s in ratioNumsWithData.split(","):
                #print "p, s : %s,%s" % (p,s)
                # do we want a match or equality? If I have QCD in numerator but I have processes QCD and QCD_1, I will have 2 matches, and this is not what I want
                # if re.match(s,p): 
                if s==p: 
                    numkeys.append(p)
                    break

    if "data" not in pmap or ratioDen == "data": 
        #print str(pmap)
        # >= 3 instead of 4 because I might have no signal process, 
        # while I always have background as sum of everything but data (minimum two processes to make ratio of them)
        if len(pmap) >= 3 and ratioDen in pmap:   
            numkeys = []
            for p in pmap.keys():                
                for s in ratioNums.split(","):
                    #print "p, s : %s,%s" % (p,s)
                    # do we want a match or equality? If I have QCD in numerator but I have processes QCD and QCD_1, I will have 2 matches, and this is not what I want
                    # if re.match(s,p): 
                    if s==p: 
                        numkeys.append(p)
                        break
            if len(numkeys) == 0:
                return (None,None,None,None)
            # then we can overwrite total with background
            # numkey = 'signal' # I think this line should not have been here, it is not used anywhere
            total     = pmap[ratioDen]
            totalSyst = pmap[ratioDen]
        else:    
            return (None,None,None,None)
    rmin, rmax =  1,1
    if firange: rmin = marange[0] 
    if firange: rmax = marange[1]
    rmin = float(pspec.getOption("RMin",rmin))
    rmax = float(pspec.getOption("RMax",rmax))

    ratios = [] #None
    for numkey in numkeys:
        # clone and divide should be enough, but cloning keeps the palette and the zaxis title of the cloned
        # ratio = pmap[numkey].Clone("data_div"); 
        # build in the hard way
        #
        # the following doesn't work if the TH2D was defined with nx,min,max,ny,min,max (don't know why, boh)
        # create the bins looping of edges
        #xbins, ybins = pmap[numkey].GetXaxis().GetXbins(), pmap[numkey].GetYaxis().GetXbins()        
        xbins = [pmap[numkey].GetXaxis().GetBinLowEdge(i+1) for i in range(pmap[numkey].GetNbinsX()+1)] # nbins+1 because I also need the low edge of overflow bin
        ybins = [pmap[numkey].GetYaxis().GetBinLowEdge(i+1) for i in range(pmap[numkey].GetNbinsY()+1)]
        #print "xbins = " + ",".join(str(x) for x in xbins)
        #print "ybins = " + ",".join(str(x) for x in ybins)
        #print "Defining TH2 for ratio: %s / %s" % (numkey, ratioDen)
        ratioName = "ratio__{d}__{n}".format(d=ratioDen,n=numkey)
        if ROOT.gROOT.FindObject(ratioName) != None:
            ROOT.gROOT.FindObject(ratioName).Delete()
        ratio = ROOT.TH2D(ratioName,ratioName,len(xbins)-1,array('d',xbins),len(ybins)-1,array('d',ybins))
        ratio.GetXaxis().SetTitle(pmap[numkey].GetXaxis().GetTitle())
        ratio.GetYaxis().SetTitle(pmap[numkey].GetYaxis().GetTitle())
        for ix in range(1,ratio.GetNbinsX()+1):
            for iy in range(1,ratio.GetNbinsY()+1):
                r = 0 if total.GetBinContent(ix, iy)==0 else pmap[numkey].GetBinContent(ix, iy)/total.GetBinContent(ix, iy)
                ratio.SetBinContent(ix, iy, r)
                ratio.SetBinError  (ix, iy, r*hypot(pmap[numkey].GetBinError(ix, iy)/pmap[numkey].GetBinContent(ix, iy) if pmap[numkey].GetBinContent(ix, iy) else 0,
                                                    total.GetBinError(ix, iy)/total.GetBinContent(ix, iy) if total.GetBinContent(ix, iy) else 0))
        ratio.GetZaxis().SetRangeUser(rmin,rmax)
        ratio.GetZaxis().SetTitle(ylabel)
        ratio.GetZaxis().SetTitleFont(42)
        ratio.GetZaxis().SetTitleSize(0.055)
        ratio.GetZaxis().SetTitleOffset(1.3)
        ratio.GetZaxis().SetLabelFont(42)
        ratio.GetZaxis().SetLabelSize(0.05)
        ratio.GetZaxis().SetLabelOffset(0.007)
        ratio.SetContour(100)
        ratios.append(ratio)

    #print 'this is ratios', ratios

    return ratios

def doStatTests(total,data,test,legendCorner):
    #print "Stat tests for %s:" % total.GetName()
    #ksprob = data.KolmogorovTest(total,"XN")
    #print "\tKS  %.4f" % ksprob
    chi2l, chi2p, chi2gq, chi2lp, nb = 0, 0, 0, 0, 0
    for b in range(1,data.GetNbinsX()+1):
        oi = data.GetBinContent(b)
        ei = total.GetBinContent(b)
        dei = total.GetBinError(b)
        if ei <= 0: continue
        nb += 1
        chi2l += - 2*(oi*log(ei/oi)+(oi-ei) if oi > 0 else -ei)
        chi2p += (oi-ei)**2 / ei
        chi2gq += (oi-ei)**2 /(ei+dei**2)
        #chi2lp +=
    logging.info("\tc2p  %.4f (%6.2f/%3d)" % (ROOT.TMath.Prob(chi2p,  nb), chi2p,  nb))
    logging.info("\tc2l  %.4f (%6.2f/%3d)" % (ROOT.TMath.Prob(chi2l,  nb), chi2l,  nb))
    logging.info("\tc2qg %.4f (%6.2f/%3d)" % (ROOT.TMath.Prob(chi2gq, nb), chi2gq, nb))
    #print "\tc2lp %.4f (%6.2f/%3d)" % (ROOT.TMath.Prob(chi2lp, nb), chi2lp, nb)
    chi2s = { "chi2l":chi2l, "chi2p":chi2p, "chi2gq":chi2gq, "chi2lp":chi2lp }
    if test in chi2s:
        chi2 = chi2s[test]
        pval = ROOT.TMath.Prob(chi2, nb)
        chi2fmt = ("%.1f" if nb < 10 else "%.0f") % chi2
        text = ("#chi^{2} %s/%d p-value %.3f" if pval < 0.02 else "#chi^{2} %s/%d p-value %.2f") % (chi2fmt, nb, pval)
    else:
        text = "Unknown test %s" % test
    if legendCorner == "TR":
        doSpam(text, .23, .85, .48, .93, align=12, textSize=0.05)
    elif legendCorner == "TL":
        doSpam(text, .75, .85, .93, .93, align=32, textSize=0.05)



legend_ = None;
def doLegend(pmap,mca,corner="TR",textSize=0.035,cutoff=1e-2,cutoffSignals=True,mcStyle="F",legWidth=0.18,legBorder=True,signalPlotScale=None,totalError=None,header="",doWide=False,overrideLegCoord=None,nColumns=1,allProcInLegend=False):
        if (corner == None): return
        total = sum([x.Integral() for x in pmap.values()])
        sigEntries = []; bgEntries = []
        cutoffTimesTotal = 0.0 if allProcInLegend else cutoff*total
        for p in mca.listSignals(allProcs=True):
            if mca.getProcessOption(p,'HideInLegend',False): continue
            if p in pmap and pmap[p].Integral() > (cutoff*total if cutoffSignals else 0): 
                lbl = mca.getProcessOption(p,'Label',p)
                if signalPlotScale and signalPlotScale!=1: 
                    lbl=lbl+" x "+("%d"%signalPlotScale if floor(signalPlotScale)==signalPlotScale else "%.2f"%signalPlotScale)
                myStyle = mcStyle if type(mcStyle) == str else mcStyle[1] if any(re.match(x,p) for x in options.forceFillColorNostackMode.split(",")) else mcStyle[0]
                sigEntries.append( (pmap[p],lbl,myStyle) )
        backgrounds = mca.listBackgrounds(allProcs=True)
        for p in backgrounds:
            if mca.getProcessOption(p,'HideInLegend',False): continue
            if p in pmap and pmap[p].Integral() >= cutoff*total: 
                lbl = mca.getProcessOption(p,'Label',p)
                myStyle = mcStyle if type(mcStyle) == str else mcStyle[1] if any(re.match(x,p) for x in options.forceFillColorNostackMode.split(",")) else mcStyle[0]
                bgEntries.append( (pmap[p],lbl,myStyle) )
        nentries = len(sigEntries) + len(bgEntries) + ('data' in pmap)

        (x1,y1,x2,y2) = (0.97-legWidth if doWide else .85-legWidth, .7 - textSize*max(nentries-3,0), .90, .91)
        if corner == "TR":
            (x1,y1,x2,y2) = (0.97-legWidth if doWide else 0.80 if options.veryWide else .85-legWidth, .7 - textSize*max(nentries-3,0), .90, .91)
        elif corner == "TC":
            (x1,y1,x2,y2) = (.5, .75 - textSize*max(nentries-3,0), .5+legWidth, .91)
        elif corner == "TL":
            (x1,y1,x2,y2) = (.2, .75 - textSize*max(nentries-3,0), .2+legWidth, .91)
        elif corner == "BR":
            (x1,y1,x2,y2) = (.85-legWidth, .33 + textSize*max(nentries-3,0), .90, .15)
        elif corner == "BC":
            (x1,y1,x2,y2) = (.5, .33 + textSize*max(nentries-3,0), .5+legWidth, .15)
        elif corner == "BL":
            (x1,y1,x2,y2) = (.2, .33 + textSize*max(nentries-3,0), .2+legWidth, .15)

        if overrideLegCoord != None:
            if len(overrideLegCoord.split(",")) == 4:
                [x1,y1,x2,y2] = [float(coord) for coord in overrideLegCoord.split(",")]
            else:
                logging.warning("In doLegend(): len(overrideLegCoord.split(",")) != 4 --> option will be ignored")

        #(x1,y1,x2,y2) = (0.2, 0.78,0.9,0.93)  # I just noticed that this line was added by Marc for monster plot (I was getting weird plots)
        # I think this is hardcoding stuff too much (we don't always plot monster plots ;)  )
        # and we already have a nice option to set coordinates (see option --setLegendCoordinates)
        # to set how many columns, use option --n-column-legend

        leg = ROOT.TLegend(x1,y1,x2,y2)
        leg.SetNColumns(nColumns)
        if header: leg.SetHeader(header.replace("\#", "#"))
        leg.SetFillColor(0)
        leg.SetFillColorAlpha(0,0.6)  # should make the legend semitransparent (second number is 0 for fully transparent, 1 for full opaque)
        #leg.SetFillStyle(0) # transparent legend, so it will not cover plots (markers of legend entries will cover it unless one changes the histogram FillStyle, but this has other effects on color, so better not touching the FillStyle)
        leg.SetShadowColor(0)
        if header: leg.SetHeader(header.replace("\#", "#"))       
        if not legBorder:
            leg.SetLineColor(0)
            leg.SetBorderSize(0)  # remove border  (otherwise it is drawn with a white line, visible if it overlaps with plots
        leg.SetTextFont(42)
        leg.SetTextSize(textSize)
        if 'data' in pmap: 
            leg.AddEntry(pmap['data'], mca.getProcessOption('data','Label','Data', noThrow=True), 'LPE')
        total = sum([x.Integral() for x in pmap.values()])
        for (plot,label,style) in sigEntries: leg.AddEntry(plot,label,style)
        for (plot,label,style) in  bgEntries: leg.AddEntry(plot,label,style)
        if totalError: leg.AddEntry(totalError,"total bkg. unc.","F") 
        leg.Draw()
        ## assign it to a global variable so it's not deleted
        global legend_
        legend_ = leg 
        return leg

class PlotMaker:
    def __init__(self,tdir,options):
        self._options = options
        self._dir = tdir
        ROOT.gROOT.ProcessLine(".x ccFiles/tdrstyle.cc") # keep here, otherwise plots are screwed up        
        if not options.drawStatBox:
            ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        if self._options.perBin and not "txt" in self._options.printPlots: raise RuntimeError("Error: cannot print yields per bin if txt option not given")
 
    def run(self,mca,cuts,plots,makeStack=True,makeCanvas=True):
        if self._options.wideplot: ROOT.gStyle.SetTitleYOffset(0.55)
        sets = [ (None, 'all cuts', cuts.allCuts()) ]
        for subname, title, cut in sets:
            logging.info(" cut set: %s" % title)
            cdir = self._dir
            if subname:
                if self._dir.Get(subname):
                    cdir = self._dir.Get(subname)
                else:
                    cdir = self._dir.mkdir(subname,title)
            cdir.cd()
            pspecs = plots.plots()
            if self._options.preFitData:
                matchspec = [ p for p in pspecs if p.name == self._options.preFitData ]
                if not matchspec: raise RuntimeError("Error: plot %s not found" % self._options.preFitData)
                pspecs = matchspec + [ p for p in pspecs if p.name != self._options.preFitData ]
            #print ' this is pspecs', pspecs
            pmaps = mca.getPlots(pspecs,cuts if self._options.printYieldsRDF else cut, makeSummary=True)
            if self._options.skipPlot:
                logging.info("Skip plotting: will save histograms in output file and exit")
                cdir.cd()
                for x in pmaps:
                    for p in list(x.keys()):
                        #print(">>>>> %s" % p)             
                        #print(x[p].GetName())
                        hname = x[p].GetName()
                        if hname.endswith("_background"): continue
                        logging.info("Writing histogram: %s" % hname)
                        cdir.WriteTObject(x[p],hname)
                return
            
            for ipspec,pspec in enumerate(pspecs):
                logging.info("    plot: %s" % pspec.name)
                #
                # blinding policy
                blind = pspec.getOption('Blinded','None') if 'data' in pmaps[ipspec] else 'None'
                if self._options.unblind == True: blind = 'None'
                xblind = [9e99,-9e99]
                if re.match(r'(bin|x)\s*([<>]?)\s*(\+|-)?\d+(\.\d+)?|(\+|-)?\d+(\.\d+)?\s*<\s*(bin|x)\s*<\s*(\+|-)?\d+(\.\d+)?', blind):
                    xfunc = (lambda h,b: b)             if 'bin' in blind else (lambda h,b : h.GetXaxis().GetBinCenter(b));
                    test  = eval("lambda bin : "+blind) if 'bin' in blind else eval("lambda x : "+blind) 
                    hdata = pmaps[ipspec]['data']
                    for b in range(1,hdata.GetNbinsX()+1):
                        if test(xfunc(hdata,b)):
                            #print "blinding bin %d, x = [%s, %s]" % (b, hdata.GetXaxis().GetBinLowEdge(b), hdata.GetXaxis().GetBinUpEdge(b))
                            hdata.SetBinContent(b,0)
                            hdata.SetBinError(b,0)
                            xblind[0] = min(xblind[0],hdata.GetXaxis().GetBinLowEdge(b))
                            xblind[1] = max(xblind[1],hdata.GetXaxis().GetBinUpEdge(b))
                    #print "final blinded range x = [%s, %s]" % (xblind[0],xblind[1])
                elif blind != "None":
                    raise RuntimeError("Unrecongnized value for 'Blinded' option, stopping here")
                #
                # Pseudo-data?
                if self._options.pseudoData:
                    if "data" in pmaps[ipspec]: raise RuntimeError("Can't use --pseudoData if there's also real data (maybe you want --xp data?)")
                    if "background" in self._options.pseudoData:
                        pdata = pmaps[ipspec]["background"]
                        pdata = pdata.Clone(str(pdata.GetName()).replace("_background","_data"))
                    elif "all" in self._options.pseudoData:
                        pdata = pmaps[ipspec]["background"]
                        pdata = pdata.Clone(str(pdata.GetName()).replace("_background","_data"))
                        if "signal" in pmaps[ipspec]: pdata.Add(pmaps[ipspec]["signal"])
                    else:
                        raise RuntimeError("Pseudo-data option %s not supported" % self._options.pseudoData)
                    if "asimov" not in self._options.pseudoData:
                        if "TH1" in pdata.ClassName():
                            for i in range(1,pdata.GetNbinsX()+1):
                                pdata.SetBinContent(i, ROOT.gRandom.Poisson(pdata.GetBinContent(i)))
                                pdata.SetBinError(i, sqrt(pdata.GetBinContent(i)))
                        elif "TH2" in pdata.ClassName():
                            for ix in range(1,pdata.GetNbinsX()+1):
                              for iy in range(1,pdata.GetNbinsY()+1):
                                pdata.SetBinContent(ix, iy, ROOT.gRandom.Poisson(pdata.GetBinContent(ix, iy)))
                                pdata.SetBinError(ix, iy, sqrt(pdata.GetBinContent(ix, iy)))
                        else:
                            raise RuntimeError("Can't make pseudo-data for %s" % pdata.ClassName())
                    pmaps[ipspec]["data"] = pdata
                #
                if not makeStack: 
                    for k,v in pmaps[ipspec].items():
                        if v.InheritsFrom("TH1"): v.SetDirectory(cdir) 
                        cdir.WriteTObject(v)
                    continue
                #
                stack = ROOT.THStack(pspec.name+"_stack",pspec.name)
                hists = [v for k,v in pmaps[ipspec].items() if k != 'data']
                total = hists[0].Clone(pspec.name+"_total"); total.Reset()
                totalSyst = hists[0].Clone(pspec.name+"_totalSyst"); totalSyst.Reset()
                if self._options.plotmode == "norm": 
                    if 'data' in pmaps[ipspec]:
                        total.GetYaxis().SetTitle(total.GetYaxis().GetTitle()+" (normalized)")
                    else:
                        total.GetYaxis().SetTitle("density/bin")
                    total.GetYaxis().SetDecimals(True)
                if self._options.scaleSigToData: self._sf = doScaleSigNormData(pspec,pmaps[ipspec],mca)
                if self._options.scaleBkgToData != []: self._sf = doScaleBkgNormData(pspec,pmaps[ipspec],mca,self._options.scaleBkgToData)
                elif self._options.fitData: doNormFit(pspec,pmaps[ipspec],mca)
                elif self._options.preFitData and pspec.name == self._options.preFitData:
                    doNormFit(pspec,pmaps[ipspec],mca,saveScales=True)
                #
                for k,v in pmaps[ipspec].items():
                    if v.InheritsFrom("TH1"): v.SetDirectory(cdir) 
                    cdir.WriteTObject(v)
                #
                self.printOnePlot(mca,pspec,pmaps[ipspec],
                                  xblind=xblind,
                                  makeCanvas=makeCanvas,
                                  outputDir=cdir,
                                  printDir=self._options.printDir+(("/"+subname) if subname else ""),
                                  drawBox=self._options.drawBox,
                                  contentAxisTitle=self._options.contentAxisTitle)

    def printOnePlot(self,mca,pspec,pmap,makeCanvas=True,outputDir=None,printDir=None,xblind=[9e99,-9e99],extraProcesses=[],plotmode="auto",outputName=None, drawBox=None, contentAxisTitle=None):
                options = self._options
                if printDir == None: printDir=self._options.printDir
                if outputDir == None: outputDir = self._dir
                if plotmode == "auto": plotmode = self._options.plotmode
                if outputName == None: outputName = pspec.name
                stack = ROOT.THStack(outputName+"_stack",outputName)
                hists = [v for k,v in pmap.items() if k != 'data']
                total = hists[0].Clone(outputName+"_total"); total.Reset("ICESM") # ICES is default, but does not reset maximum and minimum (need M as well)
                totalSyst = hists[0].Clone(outputName+"_totalSyst"); totalSyst.Reset("ICESM")

                if plotmode == "norm": 
                    if 'data' in pmap:
                        total.GetYaxis().SetTitle(total.GetYaxis().GetTitle()+" (normalized)")
                    else:
                        total.GetYaxis().SetTitle("density/bin")
                    total.GetYaxis().SetDecimals(True)

                for p in itertools.chain(reversed(mca.listBackgrounds(allProcs=True)), reversed(mca.listSignals(allProcs=True)), extraProcesses):
                    if p in pmap: 
                        plot = pmap[p]
                        #if plot.Integral() == 0:
                        #    print 'Warning: plotting histo %s with zero integral, there might be problems in the following'%p
                        if pspec.getOption('NormBinWidth',False):
                            # I hope it doesn't conflict with anything done below
                            plot.Scale(1.0, "width")
                            
                        if plot.Integral() < 0:
                            logging.warning('Plotting histo %s with negative integral (%f), the stack plot will probably be incorrect.'%(p,plot.Integral()))
                        if 'TH1' in plot.ClassName():
                            for b in range(1,plot.GetNbinsX()+1):
                                if plot.GetBinContent(b)<0: logging.warning('Histo %s has bin %d with negative content (%f), the stack plot will probably be incorrect.'%(p,b,plot.GetBinContent(b)))
                        elif 'TH2' in plot.ClassName():
                            for b1 in range(1,plot.GetNbinsX()+1):
                                for b2 in range(1,plot.GetNbinsY()+1):
                                    if plot.GetBinContent(b1,b2)<0: logging.warning('histo %s has bin %d,%d with negative content (%f), the stack plot will probably be incorrect.'%(p,b1,b2,plot.GetBinContent(b1,b2)))
#                        if plot.Integral() <= 0: continue
                        if mca.isSignal(p):
                            plot.Scale(options.signalPlotScale)
                        if mca.isSignal(p) and options.noStackSig == True: 
                            plot.SetLineWidth(3)
                            plot.SetLineColor(plot.GetFillColor())
                            continue 
                        #if mca.getProcessOption(p,'UseAsSystematic',False):
                        #    # will be used to draw an uncertainty band around nominal MC (works only when having a single MC shape)
                        #    # this is being implemented, check if this feature is ready
                        #    continue
                        if plotmode == "stack":
                            stack.Add(plot)
                            total.Add(plot)
                            totalSyst.Add(plot)
                            if mca.getProcessOption(p,'NormSystematic',0.0) > 0:
                                syst = mca.getProcessOption(p,'NormSystematic',0.0)
                                if "TH1" in plot.ClassName():
                                    for b in range(1,plot.GetNbinsX()+1):
                                        totalSyst.SetBinError(b, hypot(totalSyst.GetBinError(b), syst*plot.GetBinContent(b)))
                        else:
                            plot.SetLineColor(plot.GetFillColor())
                            plot.SetLineWidth(3)
                            if len(options.forceFillColorNostackMode) and any(re.match(x,p) for x in options.forceFillColorNostackMode.split(",")):
                                pass
                            else:
                                plot.SetFillStyle(0)
                            if plotmode == "norm" and (plot.ClassName()[:2] == "TH"):
                                ref = pmap['data'].Integral() if 'data' in pmap else 1.0
                                if (plot.Integral()): plot.Scale(ref/plot.Integral())
                            stack.Add(plot)
                            total.SetMaximum(max(total.GetMaximum(),1.3*plot.GetMaximum()))
                        if self._options.errors and plotmode != "stack":
                            plot.SetMarkerColor(plot.GetFillColor())
                            plot.SetMarkerStyle(21)
                            plot.SetMarkerSize(1.5)
                        else:
                            plot.SetMarkerStyle(0)

                binlabels = pspec.getOption("xBinLabels","")
                if binlabels != "" and len(binlabels.split(",")) == total.GetNbinsX():
                    blist = binlabels.split(",")
                    for i in range(1,total.GetNbinsX()+1): 
                        total.GetXaxis().SetBinLabel(i,blist[i-1]) 
                        total.GetYaxis().SetLabelSize(0.05)

                if not self._options.emptyStack and stack.GetNhists() == 0:
                    logging.error("For %s, all histograms are empty\n " % pspec.name)
                    #print 'this is stack', stack
                    return

                # define aspect ratio
                doWide = True if self._options.wideplot or pspec.getOption("Wide",False) else False
                plotformat = (1200,600) if doWide else (2400,600) if options.veryWide else (options.setCanvasSize[0],options.setCanvasSize[1])
                sf = 20./plotformat[0]
                # there are way too things set automatically depending on settings that the user might want to customize
                # should probably allow the user to set any of them, with some sensible default values
                if doWide or options.veryWide:
                    ROOT.gStyle.SetPadLeftMargin(600.*0.18/plotformat[0])
                else:
                    ROOT.gStyle.SetPadLeftMargin(0.18)
                    if "TH2" in total.ClassName():
                        ROOT.gStyle.SetPadLeftMargin(0.16)                        
                    if plotformat[0] > 1000: 
                        ROOT.gStyle.SetPadLeftMargin(0.05)
                        ROOT.gStyle.SetPadRightMargin(0.02)

                    

                stack.Draw("GOFF")
                ytitle = "Events" if not self._options.printBinning else "Events / %s" %(self._options.printBinning)
                if plotmode == "norm":
                    ytitle = "Arbitrary units"
                total.GetXaxis().SetTitleFont(42)
                total.GetXaxis().SetTitleSize(0.05)
                total.GetXaxis().SetTitleOffset(options.setTitleXoffset)
                total.GetXaxis().SetLabelFont(42)
                total.GetXaxis().SetLabelSize(0.05)
                total.GetXaxis().SetLabelOffset(0.007)
                total.GetYaxis().SetTitleFont(42)
                total.GetYaxis().SetTitleSize(0.05)
                total.GetYaxis().SetTitleOffset(0.90 if doWide else 1.7)
                if options.setTitleYoffset > 0: total.GetYaxis().SetTitleOffset(options.setTitleYoffset)
                total.GetYaxis().SetLabelFont(42)
                total.GetYaxis().SetLabelSize(0.05)
                total.GetYaxis().SetLabelOffset(0.007)
                total.GetYaxis().SetTitle(pspec.getOption('YTitle',ytitle))
                total.GetXaxis().SetTitle(pspec.getOption('XTitle',outputName))
                total.GetXaxis().SetNdivisions(pspec.getOption('XNDiv',510))
                if contentAxisTitle != None:
                    if "TH2" in total.ClassName() or "TProfile2D" in total.ClassName():
                        total.GetZaxis().SetTitle(contentAxisTitle)
                    else:
                        total.GetYaxis().SetTitle(contentAxisTitle)
                if outputDir: outputDir.WriteTObject(stack)
                # 
                if not makeCanvas and not self._options.printPlots: return
                #doRatio = self._options.showRatio and ('data' in pmap or (plotmode != "stack" ) and ("TH2" not in total.ClassName())
                doRatio = self._options.showRatio and ("TH2" not in total.ClassName())
                islog = pspec.hasOption('Logy'); 
                if doRatio: ROOT.gStyle.SetPaperSize(20.,sf*(plotformat[1]+150))
                else:       ROOT.gStyle.SetPaperSize(20.,sf*plotformat[1])
                # create canvas
                height = plotformat[1]+150 if doRatio else plotformat[1]
                width = plotformat[0]+150 if "TH2" in total.ClassName() else plotformat[0]
                c1 = ROOT.TCanvas(outputName+"_canvas", outputName, width, height)
                c1.SetTopMargin(c1.GetTopMargin()*options.topSpamSize);
                topsize = 0.12*600./height if doRatio else 0.06*600./height
                if self._options.doOfficialCMS: c1.SetTopMargin(topsize*1.2 if doWide else topsize)
                c1.Draw()
                p1, p2 = c1, None # high and low panes
                # set borders, if necessary create subpads
                c1.SetWindowSize(width + (width - c1.GetWw()), (height + (height - c1.GetWh())));
                if doRatio:
                    p1 = ROOT.TPad("pad1","pad1",0,0.30,1,1);
                    p1.SetTopMargin(p1.GetTopMargin()*options.topSpamSize);
                    p1.SetBottomMargin(0.025);
                    p1.Draw();
                    p2 = ROOT.TPad("pad2","pad2",0,0,1,0.30);
                    p2.SetTopMargin(0.06);
                    p2.SetBottomMargin(0.3);
                    p2.SetFillStyle(0);
                    p2.Draw();
                    p1.cd();

                p1.SetLogy(islog)
                p1.SetLogz(True if pspec.hasOption('Logz') else False)
                if pspec.hasOption('Logx'):
                    p1.SetLogx(True)
                    if p2: p2.SetLogx(True)
                    total.GetXaxis().SetNoExponent(True)
                    total.GetXaxis().SetMoreLogLabels(True)
                if islog: total.SetMaximum(2*total.GetMaximum())
                if not islog: total.SetMinimum(0)
                total.Draw("HIST")
                if plotmode == "stack":
                    stack.Draw("SAME HIST")
                    total.Draw("AXIS SAME")
                else: 
                    if self._options.errors:
                        ROOT.gStyle.SetErrorX(0.5)
                        stack.Draw("SAME E NOSTACK")
                    else:
                        stack.Draw("SAME HIST NOSTACK")
                if pspec.getOption('MoreY',1.0) > 1.0:
                    total.SetMaximum(pspec.getOption('MoreY',1.0)*total.GetMaximum())
                totalError=None
                if options.showMCError:
                    totalError = doShadedUncertainty(totalSyst)
                is2D = total.InheritsFrom("TH2")
                if 'data' in pmap: 
                    if options.poisson and not is2D:
                        pdata = getDataPoissonErrors(pmap['data'], False, True)
                        pdata.Draw("PZ SAME")
                        pmap['data'].poissonGraph = pdata ## attach it so it doesn't get deleted
                    else:
                        pmap['data'].Draw("E SAME")
                    reMax(total,pmap['data'],islog,doWide=doWide)
                    if xblind[0] < xblind[1]:
                        blindbox = ROOT.TBox(xblind[0],total.GetYaxis().GetXmin(),xblind[1],total.GetMaximum())
                        blindbox.SetFillColor(ROOT.kBlue+3)
                        blindbox.SetFillStyle(3944)
                        blindbox.Draw()
                        xblind.append(blindbox) # so it doesn't get deleted
                    if options.doStatTests:
                        doStatTests(totalSyst, pmap['data'], options.doStatTests, legendCorner=pspec.getOption('Legend','TR'))
                if pspec.hasOption('YMin') and pspec.hasOption('YMax'):
                    total.GetYaxis().SetRangeUser(pspec.getOption('YMin',1.0), pspec.getOption('YMax',1.0))
                if pspec.hasOption('ZMin') and pspec.hasOption('ZMax'):
                    total.GetZaxis().SetRangeUser(pspec.getOption('ZMin',1.0), pspec.getOption('ZMax',1.0))
                #if options.yrange: 
                #    total.GetYaxis().SetRangeUser(options.yrange[0], options.yrange[1])
                legendCutoff = pspec.getOption('LegendCutoff', 1e-5 if c1.GetLogy() else 1e-2)                
                if plotmode == "norm": legendCutoff = 0 
                if plotmode == "stack":
                    if options.noStackSig: mcStyle = ("L","F")
                    else:                  mcStyle = "F"
                else: mcStyle = ("L","F") if len(options.forceFillColorNostackMode) else "L"
                if options.allProcInLegend: legendCutoff = 0.0
                doLegend(pmap,mca,corner=pspec.getOption('Legend','TR'),
                         cutoff=legendCutoff, mcStyle=mcStyle,
                         cutoffSignals=not(options.showSigShape or options.showIndivSigShapes or options.showSFitShape), 
                         textSize=( (0.045 if doRatio else 0.022) if options.legendFontSize <= 0 else options.legendFontSize ),
                         legWidth=options.legendWidth, legBorder=options.legendBorder, signalPlotScale=options.signalPlotScale,
                         header=self._options.legendHeader if self._options.legendHeader else pspec.getOption("LegendHeader", ""),
                         doWide=doWide, totalError=totalError,
                         overrideLegCoord=self._options.setLegendCoordinates, nColumns=self._options.nColumnLegend)
                if self._options.doOfficialCMS:
                    CMS_lumi.lumi_13TeV = "%.1f fb^{-1}" % self._options.lumi
                    CMS_lumi.extraText  = self._options.cmsprel
                    CMS_lumi.lumi_sqrtS = self._options.cmssqrtS
                    CMS_lumi.CMS_lumi(ROOT.gPad, 4, 0, -0.005 if doWide and doRatio else 0.01 if doWide else 0.05)
                else: 
                    rightTextOffset = 0.04 - p1.GetRightMargin() -0.02 # 0.04 should be the default, then subtract other 0.02
                    #print "rightTextOffset = ", str(rightTextOffset)
                    # xoffs should be negative if the right text goes outside the right margin of the canvas
                    doTinyCmsPrelim(hasExpo = total.GetMaximum() > 9e4 or c1.GetLogy(),textSize=(0.040 if doRatio else 0.033)*options.topSpamSize, options=options,doWide=doWide, xoffs=rightTextOffset)  
                if options.addspam:
                    if pspec.getOption('Legend','TR')=="TL":
                        doSpam(options.addspam, .68, .855, .9, .895, align=32, textSize=(0.040 if doRatio else 0.033)*options.topSpamSize)
                    else:
                        doSpam(options.addspam, .23, .855, .6, .895, align=12, textSize=(0.040 if doRatio else 0.033)*options.topSpamSize)
                signorm = None; datnorm = None; sfitnorm = None
                if options.showSigShape or options.showIndivSigShapes or options.showIndivSigs: 
                    signorms = doStackSignalNorm(pspec,pmap,options.showIndivSigShapes or options.showIndivSigs,extrascale=options.signalPlotScale, norm=not options.showIndivSigs)
                    for signorm in signorms:
                        if outputDir: 
                            signorm.SetDirectory(outputDir); outputDir.WriteTObject(signorm)
                        reMax(total,signorm,islog,doWide=doWide)
                if options.showDatShape: 
                    datnorm = doDataNorm(pspec,pmap)
                    if datnorm != None:
                        if outputDir: 
                            datnorm.SetDirectory(outputDir); outputDir.WriteTObject(datnorm)
                        reMax(total,datnorm,islog,doWide=doWide)
                if options.showSFitShape: 
                    (sfitnorm,sf) = doStackSigScaledNormData(pspec,pmap)
                    if sfitnorm != None:
                        if outputDir: 
                            sfitnorm.SetDirectory(outputDir); outputDir.WriteTObject(sfitnorm)
                        reMax(total,sfitnorm,islog,doWide=doWide)
                if options.flagDifferences and len(pmap) == 4:
                    new = pmap['signal']
                    ref = pmap['background']
                    if "TH1" in new.ClassName():
                        for b in range(1,new.GetNbinsX()+1):
                            if abs(new.GetBinContent(b) - ref.GetBinContent(b)) > options.toleranceForDiff*ref.GetBinContent(b):
                                logging.info("Plot: difference found in %s, bin %d" % (outputName, b))
                                p1.SetFillColor(ROOT.kYellow-10)
                                if p2: p2.SetFillColor(ROOT.kYellow-10)
                                break
                if makeCanvas and outputDir: outputDir.WriteTObject(c1)
                rdata,rnorm,rnorm2,rline = (None,None,None,None)
                if doRatio:
                    p2.cd(); 
                    rdata,rnorm,rnorm2,rline = doRatioHists(pspec,pmap,total,totalSyst, marange=options.maxRatioRange, firange=options.fixRatioRange,
                                                            fitRatio=options.fitRatio, errorsOnRef=options.errorBandOnRatio, onlyStatErrorsOnRef=options.onlyStatErrorOnRatio,
                                                            ratioNums=options.ratioNums, ratioDen=options.ratioDen, ylabel=options.ratioYLabel, 
                                                            doWide=doWide, showStatTotLegend=(False if options.noLegendRatioPlot else True),
                                                            errorBarsOnRatio = options.errorBarsOnRatio,
                                                            ratioYLabelSize = options.ratioYLabelSize,
                                                            ratioNumsWithData = options.ratioNumsWithData)
                if self._options.printPlots:
                    for ext in self._options.printPlots.split(","):
                        fdir = printDir;
                        if not os.path.exists(fdir): 
                            os.makedirs(fdir); 
                            if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/m/mciprian/public/index.php "+fdir)
                            elif os.path.exists("/pool/ciencias/"): os.system("cp /pool/ciencias/HeppyTrees/RA7/additionalReferenceCode/index.php "+fdir)
                        if ext == "txt" and self._options.perBin:
                            dump = open("%s/%s_perBin.%s" % (fdir, outputName, ext), "w")
                            maxlen = max([len(mca.getProcessOption(p,'Label',p)) for p in mca.listSignals(allProcs=True) + mca.listBackgrounds(allProcs=True)]+[7])
                            step = (plot.GetXaxis().GetXmax()-plot.GetXaxis().GetXmin())/plot.GetNbinsX()
                            bins = [plot.GetXaxis().GetXmin()+i*step for i in range(plot.GetNbinsX())]
                            fmh    = "%%-%ds" % (maxlen+1)
                            fmt    = "%9.2f +/- %9.2f (stat)"
                            dump.write(fmh % pspec.expr + " " + " ".join("%d" % (i) for i in bins) + "\n")
                            dump.write(("-"*(maxlen+45))+"\n");
                            bkgsyst = [0 for i in range(pmap["background"].GetNbinsX())]; sigsyst = bkgsyst
                            for p in mca.listSignals(allProcs=True) + mca.listBackgrounds(allProcs=True) + ["signal", "background"]:
                                if p not in pmap: continue
                                plot = pmap[p]
                                if plot.Integral() <= 0: continue
                                norm = plot.Integral()
                                if p not in ["signal","background"] and mca.isSignal(p): norm /= options.signalPlotScale # un-scale what was scaled
                                if p == "signal": dump.write(("-"*(maxlen+45))+"\n");
                                dump.write(fmh % (_unTLatex(mca.getProcessOption(p,'Label',p) if p not in ["signal", "background"] else p.upper())))
                                bins = []
                                for b in range(1,plot.GetNbinsX()+1):
                                    syst = plot.GetBinContent(b) * mca.getProcessOption(p,'NormSystematic',0.0) if p not in ["signal", "background"] else 0;
                                    if p in mca.listBackgrounds(allProcs=True): bkgsyst[b-1] += syst*syst 
                                    if p in mca.listSignals(allProcs=True)    : sigsyst[b-1] += syst*syst
                                    line = fmt % (plot.GetBinContent(b), plot.GetBinError(b))
                                    #if syst: line += " +/- %9.2f (syst)"  % syst
                                    if   p == "signal"     and sigsyst[b-1]: line += " +/- %9.2f (syst)" % math.sqrt(sigsyst[b-1])
                                    elif p == "background" and bkgsyst[b-1]: line += " +/- %9.2f (syst)" % math.sqrt(bkgsyst[b-1])
                                    else: line += " +/- %9.2f (syst)"  % syst
                                    bins.append(line)
                                dump.write(" ".join(bins) + "\n")
                            if 'data' in pmap: 
                                dump.write(("-"*(maxlen+45))+"\n");
                                dump.write("%%%ds " % (maxlen+1) % ('DATA'))
                                plot = pmap['data']
                                bins = []
                                for b in range(1,plot.GetNbinsX()+1):
                                    bins.append("%7.0f" % plot.GetBinContent(b))
                                dump.write(" ".join(bins) + "\n")
                            for logname, loglines in pspec.allLogs():
                                dump.write("\n\n --- %s --- \n" % logname)
                                for line in loglines: dump.write("%s\n" % line)
                                dump.write("\n")
                            dump.close()
                        if ext == "txt":
                            dump = open("%s/%s.%s" % (fdir, outputName, ext), "w")
                            maxlen = max([len(mca.getProcessOption(p,'Label',p)) for p in mca.listSignals(allProcs=True) + mca.listBackgrounds(allProcs=True)]+[7])
                            fmt    = "%%-%ds %%9.2f +/- %%9.2f (stat)" % (maxlen+1)
                            for p in mca.listSignals(allProcs=True) + mca.listBackgrounds(allProcs=True) + ["signal", "background"]:
                                if p not in pmap: continue
                                plot = pmap[p]
                                if plot.Integral() <= 0: continue
                                norm = plot.Integral()
                                if p not in ["signal","background"] and mca.isSignal(p): norm /= options.signalPlotScale # un-scale what was scaled
                                stat = sqrt(sum([plot.GetBinError(b)**2 for b in range(1,plot.GetNbinsX()+1)]))
                                syst = norm * mca.getProcessOption(p,'NormSystematic',0.0) if p not in ["signal", "background"] else 0;
                                if p == "signal": dump.write(("-"*(maxlen+45))+"\n");
                                dump.write(fmt % (_unTLatex(mca.getProcessOption(p,'Label',p) if p not in ["signal", "background"] else p.upper()), norm, stat))
                                if syst: dump.write(" +/- %9.2f (syst)"  % syst)
                                dump.write("\n")
                            if 'data' in pmap: 
                                dump.write(("-"*(maxlen+45))+"\n");
                                dump.write(("%%%ds %%7.0f\n" % (maxlen+1)) % ('DATA', pmap['data'].Integral()))
                            for logname, loglines in pspec.allLogs():
                                dump.write("\n\n --- %s --- \n" % logname)
                                for line in loglines: dump.write("%s\n" % line)
                                dump.write("\n")
                            dump.close()
                        else:
                            savErrorLevel = ROOT.gErrorIgnoreLevel; ROOT.gErrorIgnoreLevel = ROOT.kWarning;
                            if "TH2" in total.ClassName() or "TProfile2D" in total.ClassName():
                                pmap["total"] = total
                                #stuffToPlot = mca.listSignals(allProcs=True) + mca.listBackgrounds(allProcs=True)                          
                                for p in mca.listSignals(allProcs=True) + mca.listBackgrounds(allProcs=True) + ["signal", "background", "data", "total"]:
                                    if p not in pmap: continue
                                    plot = pmap[p]
                                    if "TGraph" in plot.ClassName(): continue
                                    c1.SetRightMargin(0.20)
                                    plot.SetContour(100)
                                    ROOT.gStyle.SetPalette(options.palette)
                                    ROOT.gStyle.SetPaintTextFormat(pspec.getOption("PaintTextFormat","g"))
                                    plot.SetMarkerSize(pspec.getOption("MarkerSize",1))
                                    if pspec.hasOption('ZMin') and pspec.hasOption('ZMax'):
                                        plot.GetZaxis().SetRangeUser(pspec.getOption('ZMin',1.0), pspec.getOption('ZMax',1.0))
                                    # plot.SetMarkerStyle(mca.getProcessOption(p,'MarkerStyle',1))
                                    # plot.SetMarkerColor(mca.getProcessOption(p,'FillColor',ROOT.kBlack))
                                    #plot.Draw(pspec.getOption("PlotMode","COLZ TEXT45"))
                                    plot.Draw(pspec.getOption("PlotMode","COLZ"))
                                    # Z axis setting ######################
                                    plot.GetZaxis().SetTitle(pspec.getOption('ZTitle',ytitle)) # use same content of default ytitle defined above, Events or Events/XX
                                    plot.GetZaxis().SetTitleFont(42)
                                    plot.GetZaxis().SetTitleSize(0.055)
                                    plot.GetZaxis().SetTitleOffset(0.90 if doWide else 1.3) 
                                    plot.GetZaxis().SetLabelFont(42)
                                    plot.GetZaxis().SetLabelSize(0.05)
                                    plot.GetZaxis().SetLabelOffset(0.007)
                                    #plot.GetZaxis().SetMaxDigits(2)  # avoid too many digits: with 2, for Nz > 99 the scientific notation is used, does not work in 8_0_25
                                    if contentAxisTitle != None: plot.GetZaxis().SetTitle(contentAxisTitle)
                                    #########################################
                                    if drawBox != None:
                                        liner = ROOT.TLine()
                                        liner.SetLineColor(options.boxColor)
                                        liner.SetLineWidth(3)
                                        x1,x2,y1,y2 = (float(x) for x in drawBox.split(','))
                                        liner.DrawLine(x1, y1, x1, y2)  # vertical at x=x1, y in [y1,y2]
                                        liner.DrawLine(x2, y1, x2, y2)  # vertical at x=x2, y in [y1,y2]
                                        liner.DrawLine(x1, y1, x2, y1)  # horizontal at y=y1, x in [x1,x2]
                                        liner.DrawLine(x1, y2, x2, y2)  # horizontal at y=y2, x in [x1,x2]
                                    c1.Print("%s/%s_%s.%s" % (fdir, outputName, p, ext))
                                if "data" in pmap and "TGraph" in pmap["data"].ClassName():
                                    pmap["data"].SetMarkerStyle(mca.getProcessOption('data','MarkerStyle',1))
                                    pmap["data"].SetMarkerSize(pspec.getOption("MarkerSize",1.6))
                                    for p in ["signal", "background", "total"]:
                                        if p not in pmap: continue
                                        plot = pmap[p]
                                        c1.SetRightMargin(0.20)
                                        plot.SetContour(100)
                                        #plot.Draw(pspec.getOption("PlotMode","COLZ TEXT45"))
                                        plot.Draw(pspec.getOption("PlotMode","COLZ"))
                                        pmap["data"].Draw("P SAME")
                                        c1.Print("%s/%s_data_%s.%s" % (fdir, outputName, p, ext))
                                if self._options.showRatio and ("TH2" in total.ClassName()):
                                    # following function is called twice, because we are looping on printplots = pdf and png
                                    # inside the function some histograms are created, which means they are being redefined
                                    # this will issue a warning that a histogram with same name is being replaced, but should be harmless
                                    rdata = doRatio2DHists(pspec,pmap,total,totalSyst, marange=options.maxRatioRange, firange=options.fixRatioRange,
                                                           ratioNums=options.ratioNums, ratioDen=options.ratioDen, ylabel=options.ratioYLabel,ratioNumsWithData=options.ratioNumsWithData)
                                    for r in rdata:
                                        if r == None:
                                            continue
                                        r.Draw(pspec.getOption("PlotMode","COLZ0"))
                                        if drawBox != None:
                                            # liner = ROOT.TLine()
                                            # liner.SetLineColor(options.boxColor)
                                            # liner.SetLineWidth(3)
                                            # x1,x2,y1,y2 = (float(x) for x in drawBox.split(','))
                                            liner.DrawLine(x1, y1, x1, y2)  # vertical at x=x1, y in [y1,y2]
                                            liner.DrawLine(x2, y1, x2, y2)  # vertical at x=x2, y in [y1,y2]
                                            liner.DrawLine(x1, y1, x2, y1)  # horizontal at y=y1, x in [x1,x2]
                                            liner.DrawLine(x1, y2, x2, y2)  # horizontal at y=y2, x in [x1,x2]
                                        # no log scale for ratio
                                        if pspec.getOption('Logz',False):
                                            c1.SetLogz(0) 
                                        c1.Print("%s/%s_%s.%s" % (fdir, outputName, r.GetName(), ext))
                                        # reset to log scale if modified above, even though at this point the plot is made and for the next one the option is specifically evaluated again
                                        if pspec.getOption('Logz',False):
                                            c1.SetLogz() 
                                        r.Write(r.GetName())
                            else:
                                c1.Print("%s/%s.%s" % (fdir, outputName, ext))
                            ROOT.gErrorIgnoreLevel = savErrorLevel;
                c1.Close()


def addPlotMakerOptions(parser, addAlsoMCAnalysis=True):
    if addAlsoMCAnalysis: addMCAnalysisOptions(parser)
    parser.add_argument("-l", "--lumi", type=float, default="1", help="Luminosity (in 1/fb). Only for plots, not as weight")
    parser.add_argument("--ss",  "--scale-signal", dest="signalPlotScale", default=1.0, type=float, help="scale the signal in the plots by this amount");
    parser.add_argument("--lspam", type=str, default="#bf{CMS} #it{Preliminary}", help="Spam text on the left hand side");
    parser.add_argument("--rspam", type=str, default="%(lumi) (13 TeV)", help="Spam text on the right hand side");
    parser.add_argument("--addspam", type=str, help="Additional spam text on the top left side, in the frame");
    parser.add_argument("--topSpamSize", type=float, default=1.2, help="Zoom factor for the top spam");
    parser.add_argument("--print", dest="printPlots", type=str, default="png,pdf,txt", help="print out plots in this format or formats (e.g. 'png,pdf,txt')");
    parser.add_argument("--pdir", "--print-dir", dest="printDir", type=str, default="plots", help="print out plots in this directory");
    parser.add_argument("--showSigShape", action="store_true", help="Superimpose a normalized signal shape")
    parser.add_argument("--showIndivSigShapes", action="store_true", help="Superimpose normalized shapes for each signal individually")
    parser.add_argument("--showIndivSigs", action="store_true", help="Superimpose shapes for each signal individually (normalized to their expected event yield)")
    parser.add_argument("--noStackSig", action="store_true", help="Don't add the signal shape to the stack (useful with --showSigShape)")
    parser.add_argument("--showDatShape", action="store_true", help="Stack a normalized data shape")
    parser.add_argument("--showSFitShape", action="store_true", help="Stack a shape of background + scaled signal normalized to total data")
    parser.add_argument("--showMCError", action="store_true", help="Show a shaded area for MC uncertainty")
    parser.add_argument("--showRatio", action="store_true", help="Add a data/sim ratio plot at the bottom")
    parser.add_argument("--ratioDen", type=str, default="background", help="Denominator of the ratio, when comparing MCs. Default is %(default)s")
    parser.add_argument("--ratioNumsWithData", type=str, default="", help="When plotting data and MC, use also these processes as numerators to make ratio with total/background (useful to plot ratios of unstacked components). Need a comma separated list of processes");
    parser.add_argument("--ratioNums", type=str, default="signal", help="Numerator(s) of the ratio, when comparing MCs (comma separated list of regexps)")
    parser.add_argument("--ratioYLabel", type=str, default="Data/pred.", help="Y axis label of the ratio histogram.")
    parser.add_argument("--noErrorBandOnRatio", dest="errorBandOnRatio", action="store_false", help="Do not show the error band on the reference in the ratio plots")
    parser.add_argument("--onlyStatErrorOnRatio", dest="onlyStatErrorOnRatio", action="store_true", help="Show only stat error on reference in ratio plots (when --noErrorBandOnRatio is False)")
    parser.add_argument("--noErrorBarsOnRatio", dest="errorBarsOnRatio", action="store_false", help="Do not show the error bars on the ratio plots (affect each numerator, unlike --noErrorBandOnRatio")
    parser.add_argument("--fitRatio", type=int, help="Fit the ratio with a polynomial of the specified order")
    parser.add_argument("--scaleSigToData", action="store_true", help="Scale all signal processes so that the overall event yield matches the observed one")
    parser.add_argument("--scaleBkgToData", action="append", default=[], help="Scale all background processes so that the overall event yield matches the observed one")
    parser.add_argument("--showSF", action="store_true", help="Show scale factor extracted from either --scaleSigToData or --scaleBkgToData on the plot")
    parser.add_argument("--fitData", action="store_true", help="Perform a fit to the data")
    parser.add_argument("--preFitData", type=str, help="Perform a pre-fit to the data using the specified distribution, then plot the rest")
    parser.add_argument("--maxRatioRange", type=float, nargs=2, default=(0.0, 5.0), help="Min and max for the ratio")
    parser.add_argument("--fixRatioRange", action="store_true", help="Fix the ratio range to --maxRatioRange")
    parser.add_argument("--doStatTests", type=str, help="Do this stat test: chi2p (Pearson chi2), chi2l (binned likelihood equivalent of chi2)")
    parser.add_argument("--plotmode", type=str, default="stack", help="Show as stacked plot (stack), a non-stacked comparison (nostack) and a non-stacked comparison of normalized shapes (norm)")
    parser.add_argument("--rebin", dest="globalRebin", type=int, default=0, help="Rebin all plots by this factor")
    parser.add_argument("--no-poisson", dest="poisson", action="store_false", help="Don't draw Poisson error bars")
    parser.add_argument("--unblind", action="store_true", help="Unblind plots irrespectively of plot file")
    parser.add_argument("--select-plot", "--sP", dest="plotselect", action="append", default=[], help="Select only these plots out of the full file")
    parser.add_argument("--exclude-plot", "--xP", dest="plotexclude", action="append", default=[], help="Exclude these plots from the full file")
    parser.add_argument("--legendWidth", type=float, default=0.25, help="Width of the legend")
    parser.add_argument("--legendBorder", type=int, default=0, help="Use a border in the legend (1=yes, 0=no)")
    parser.add_argument("--legendFontSize", type=float, default=0.055, help="Font size in the legend (if <=0, use the default)")
    parser.add_argument("--flagDifferences", action="store_true", default=False, help="Flag plots that are different (when using only two processes, and plotmode nostack")
    parser.add_argument("--toleranceForDiff", default=0.0, type=float, help="set numerical tollerance to define when two histogram bins are considered different");
    parser.add_argument("--pseudoData", type=str, default=None, help="If set to 'background' or 'all', it will plot also a pseudo-dataset made from background (or signal+background) with Poisson fluctuations in each bin.")
    parser.add_argument("--wide", dest="wideplot", action="store_true", help="Draw a wide canvas")
    parser.add_argument("--verywide", dest="veryWide", action="store_true", help="Draw a very wide canvas")
    parser.add_argument("--canvasSize", dest="setCanvasSize", type=int, nargs=2, default=(600, 600), help="Set canvas height and width")
    parser.add_argument("--setTitleYoffset", type=float, default=-1.0, help="Set Y axis offset for title (must be >0, if <0 the default is used)")
    parser.add_argument("--setTitleXoffset", type=float, default=0.90, help="Set X axis offset for title. The default is 0.9, which is fine unless there are superscripts, in that case 1.1 is suggested. It should be tuned based on the canvas size and presence of ratio plot.")
    parser.add_argument("--yrange", nargs=2, type=float, help="Y axis range");
    parser.add_argument("--emptyStack", action="store_true", help="Allow empty stack in order to plot, for example, only signals but no backgrounds.")
    parser.add_argument("--noLegendRatioPlot", action="store_true", help="Remove legend in ratio plot (by default it is drawn)");
    parser.add_argument("--perBin", action="store_true", help="Print the contents of every bin in another txt file");
    parser.add_argument("--legendHeader", type=str, help="Put a header to the legend")
    parser.add_argument("--ratioOffset", type=float, help="Put an offset between ratio and main pad")
    parser.add_argument("--ratioYLabelSize", type=float, default=0.14, help="Set size of title for y axis in ratio (note that it might require adjusting the offset")
    parser.add_argument("--noCms", dest="doOfficialCMS", action="store_false", help="Use official tool to write CMS spam")
    parser.add_argument("--cmsprel", type=str, default="Preliminary", help="Additional text (Simulation, Preliminary, Internal)")
    parser.add_argument("--cmssqrtS", type=str, default="13 TeV", help="Sqrt of s to be written in the official CMS text.")
    parser.add_argument("--printBin", dest="printBinning", type=str, help="Write 'Events/xx' instead of 'Events' on the y axis")
    parser.add_argument("--drawBox", type=str, help="For TH2: draw box with passing comma separated list of coordinates x1,x2,y1,y2. Example: --drawBox 'x1,x2,y1,y2'")
    parser.add_argument("--box-color", dest="boxColor", type=int, default=1, help="Set line color for box (works with --drawBox). Default is black")
    parser.add_argument("--contentAxisTitle", type=str, help="Set name of axis with bin content (Y for TH1, Z for TH2), overriding the one set by default or in the MCA file")
    parser.add_argument("--setLegendCoordinates", type=str, help="Pass a comma separated list of 4 numbers (x1,y1,x2,y2) used to build the legend. It overrides options to set the legend position in MCA, so this option is more suited to create a single plot")
    parser.add_argument("--n-column-legend", dest="nColumnLegend", type=int, default=1, help="Number of columns for legend (for monster plot, 3 is better)")
    parser.add_argument("--updateRootFile", action="store_true", help="Open the root file in UPDATE more (useful when you want to add a new histogram without running all the others)");
    parser.add_argument("--palette", type=int, default=57, help="Set color palette (only for 2D histograms)")
    parser.add_argument("--allProcInLegend", action="store_true", help="Put all processes in legend, regardless their integral.")
    parser.add_argument("--forceFillColorNostackMode", type=str, default="", help="Use fill color and style defined in MCA file when using --plotmode nostack|norm (comma separated list of regexps, by default only lines are used).")
    parser.add_argument("--drawStatBox", action="store_true", help="Draw stat box");
    parser.add_argument("-o", "--out", type=str, help="Output file name. by default equal to plots -'.txt' +'.root'");
    parser.add_argument("--no-rdf-report", dest="printYieldsRDF", action="store_false", help="Use RDF Report functionality to print yields per process (requires multiple filters, one for each line in cut file)")
    parser.add_argument("--yields-outfile", dest="yieldsOutfile", type=str, help="Output file name for yields. By default equal to PlotFile + '_yields' with extension set automatically by option --tf/--text-format");
    parser.add_argument("--skipPlot", action="store_true", default=False, help="After making the histograms save them in output file and exit, skipping plotting (usually when making histograms for datacards)")
    parser.add_argument("plotFile", type=str, help="Text file with plot format specifications")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    addPlotMakerOptions(parser)
    args = parser.parse_args()
    # TODO: The fact that options is in the global scope is kind of abused. 
    # Would be better to pass this to the functions that need it. For now, just leaving it
    options = args
    setLogging(args.verbose)
    mca  = MCAnalysis(args.sampleFile, args)
    cuts = CutsFile(args.cutFile, args)
    plots = PlotFile(args.plotFile, args)
    outname  = args.out if args.out else (args.plotFile.replace(".txt","")+".root")
    if (not args.out) and args.printDir:
        outname = args.printDir + "/"+os.path.basename(args.plotFile.replace(".txt","")+".root")
    outdir = os.path.dirname(outname) 
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
        htmlpath = "./templates/index.php"
        shutil.copy(htmlpath, outdir)
    logging.info("Will save plots to %s " % outname)
    fcmd = open(re.sub("\.root$","",outname)+"_command.txt","w")
    fcmd.write("%s\n\n" % " ".join(sys.argv))
    fcmd.close()
    shutil.copy(args.plotFile, re.sub("\.root$","",outname)+"_plots.txt")
    shutil.copy(args.sampleFile, re.sub("\.root$","",outname)+"_mca.txt")
    if args.yieldsOutfile:
        args.yieldsOutfile = args.yieldsOutfile + args.txtfmt
    else:
        args.yieldsOutfile = re.sub("\.root$","",outname)+"_yields." + args.txtfmt
    fcut = open(re.sub("\.root$","",outname)+"_cuts.txt","w")
    fcut.write("%s\n" % cuts);
    if args.rdfDefineFile or len(args.rdfDefine) or len(args.rdfAlias):        
        #frdfdefine = open(re.sub("\.root$","",outname)+"_rdfdefine.txt","w")
        fcut.write("\n\n")
        fcut.write("## Defines\n")
        lines = []
        if args.rdfDefineFile:
            with open(args.rdfDefineFile) as f:
                lines = [x.strip() for x in f if not x.startswith("#") and len(x) > 0]
        defs = [x.split(":")[0].strip() for x in lines]
        for l in args.rdfDefine:
            defname = l.split(":")[0].strip()
            for i,d in enumerate(defs):
                if defname == d:
                    lines[i] = l
                else:
                    if defname not in defs:
                        lines.append(l)
        for l in lines:
            fcut.write("%s\n" % l)
        fcut.write("## Aliases\n")
        for l in args.rdfAlias:
            fcut.write("%s\n" % l)
        #frdfdefine.close()
    fcut.close()
    rootFileOpenMode = "UPDATE" if args.updateRootFile else "RECREATE"
    outfile  = ROOT.TFile(outname,rootFileOpenMode)
    plotter = PlotMaker(outfile, args)
    plotter.run(mca,cuts,plots)
    outfile.Close()
    print(f"Histograms saved in file {outname}\n")

