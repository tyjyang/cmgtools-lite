#!/bin/env python

# Path of some root files
# /afs/cern.ch/user/m/mdunser/public/triggerTnP_muons_fullData.root
# official_muon_ScaleFactors_2016/
# /afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/data/scaleFactors/muons/

################################
# Exaples
################################

# muon trigger scale factors
#    python w-helicity-13TeV/smoothLeptonScaleFactors.py -i /afs/cern.ch/user/m/mdunser/public/triggerTnP_muons_fullData.root -o ~/www/wmass/13TeV/scaleFactors/muon/trigger/ -c mu -n smoothEfficiency.root -t

# muon isolation (replace ISO with ID for ID scale factors)
#    python w-helicity-13TeV/smoothLeptonScaleFactors.py -i /afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/data/scaleFactors/muons/efficienciesISO.root -o ~/www/wmass/13TeV/scaleFactors/muon/efficiency_ISO_GH/ -c mu -n smoothEfficiency.root -e GH -v ISO

# electron trigger
#    python w-helicity-13TeV/smoothLeptonScaleFactors.py -i ~mdunser/public/electronTriggerEfficiencies_allEras.root -o ~/www/wmass/13TeV/scaleFactors/electron/trigger/ -c el -n smoothEfficiency.root -t

################################
################################


import ROOT, os, sys, re, array, math

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gROOT.SetBatch(True)

import utilities
utilities = utilities.util()

def makeSFweightedAverage(f1, f2, w1, w2, hname, hnewname, uncertaintyRule="max"):

    # uncertaintyRule can be max, weight, diff
    # max takes the largest of the two uncertainties as the uncertainty on the mean
    # weight compute the usual uncertainty on the weighted mean
    # diff takes 1/2 of the difference as uncertainty (so the difference is the full error bars)

    tf = ROOT.TFile.Open(f1)        
    h1 =   tf.Get(hname).Clone(hname+"_1")
    if (h1 == 0):
        print "Error in makeSFweightedAverage(): could not retrieve histogram from input file %s. Exit" % f1
        quit()
    h1.SetDirectory(0)
    tf.Close()

    tf = ROOT.TFile.Open(f2)        
    h2 =   tf.Get(hname).Clone(hname+"_2")
    if (h2 == 0):
        print "Error in makeSFweightedAverage(): could not retrieve histogram from input file %s. Exit" % f2
        quit()
    h2.SetDirectory(0)
    tf.Close()

    #if not ROOT.TH1.CheckConsistency(h1,h2):
    #    print "Error in makeSFweightedAverage(): input histograms are not consistent. Exit"
    #    quit()
    xbins = h1.GetXaxis().GetXbins()
    ybins = [h1.GetYaxis().GetBinLowEdge(i) for i in range(1,h1.GetNbinsY()+2)]
    #ybins = h1.GetYaxis().GetXbins() # does not work, maybe because it has 500 elements and the maximum allowed is 255
    #print "nybins = %d" % h1.GetNbinsY()
    #print "xbins = %s (%d bins)" % (",".join(str(x) for x in xbins),len(xbins)-1)
    #print "ybins = %s" % ",".join(str(x) for x in ybins)
    hnew = ROOT.TH2D(hnewname,"",len(xbins)-1,array('d',xbins),len(ybins)-1,array('d',ybins))
    #hnew.SetDirectory(0)

    hnew.GetZaxis().SetTitle(h1.GetZaxis().GetTitle())
    hnew.SetTitle(h1.GetTitle())

    for ix in range(1,hnew.GetNbinsX()+1):
        for iy in range(1,hnew.GetNbinsY()+1):
            val1 = h1.GetBinContent(ix,iy)
            val2 = h2.GetBinContent(ix,iy)
            err1 = h1.GetBinError(ix,iy)
            err2 = h2.GetBinError(ix,iy)
            newval = (w1*val1 + w2*val2)/(w1+w2)
            if uncertaintyRule == "max":
                newerr = max(err1,err2)
            elif uncertaintyRule == "weight":
                newerr = math.sqrt(w1**2 * err1**2 + w2**2 * err2**2)/(w1+w2)
            elif uncertaintyRule == "diff":
                newerr = 0.5 * abs(val1-val2)
            else:
                print "Error in makeSFweightedAverage(): unknown uncertainty rule %s. Exit" % uncertaintyRule
                quit()
            hnew.SetBinContent(ix,iy,newval)
            hnew.SetBinError(ix,iy,newerr)

    return hnew

def copyHisto(h1, h2):

    if h1.GetDimension() != h2.GetDimension():
        print "Error in copyHisto(): histograms have different dimensions. Dim(h1)=%d  Dim(h2)=%d. Exit" % (h1.GetDimension(),h2.GetDimension())
        quit()

    if h1.GetDimension() == 1:
        for ix in range(h2.GetNbinsX()+2):
                h1.SetBinContent(ix,h2.GetBinContent(ix,iy))
                h1.SetBinError(ix,h2.GetBinError(ix,iy))
    elif h1.GetDimension() == 2:
        for ix in range(h2.GetNbinsX()+2):
            for iy in range(h2.GetNbinsY()+2):
                h1.SetBinContent(ix,iy,h2.GetBinContent(ix,iy))
                h1.SetBinError(ix,iy,h2.GetBinError(ix,iy))
    else:
        print "Error in copyHisto(): function not implemented for dimension > 2. Exit"
        quit()        


def fitTurnOn(hist, key, outname, mc, channel="el", hist_chosenFunc=0, drawFit=True, isIso=False, isTrigger=False, isFullID=False):

    #drawFit = False
    # isIso is mainly for muons, for which ID ad ISO are separate

    isEle = True if channel == "el" else False
    if not isEle: isFullID = False
    #print "isEle",str(isEle)
    #mc = "MC" if isMC else "Data"
    outdir = "{out}{mc}/".format(out=outname,mc=mc)
    createPlotDirAndCopyPhp(outdir)

    canvas = ROOT.TCanvas("canvas_%s_%s" % (mc,key),"",700,700)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.14)
    canvas.SetBottomMargin(0.12)
    canvas.SetRightMargin(0.06)
    canvas.cd()                           

    setTDRStyle()

    hist.SetLineColor(ROOT.kBlack)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1)

    hist.GetXaxis().SetTitle("%s p_{T} [GeV]" % ("electron" if channel == "el" else "muon"))
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetLabelSize(0.04)
    if mc == "SF":
        hist.GetYaxis().SetTitle("Data/MC scale factor")
    else:
        hist.GetYaxis().SetTitle("{mc} efficiency".format(mc=mc))
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetLabelSize(0.04)
    if isTrigger or isFullID: hist.GetYaxis().SetRangeUser(0.98*hist.GetMinimum(), 1.02* hist.GetMaximum())
    else: 
        diff = hist.GetMaximum() - hist.GetMinimum()
        hist.GetYaxis().SetRangeUser(hist.GetMinimum() - diff, diff + hist.GetMaximum())
    #hist.GetYaxis().SetRangeUser(0.3,1.2)
    hist.SetStats(0)
    hist.Draw("EP")

    if isTrigger or isFullID:
        fitopt = "QMFS+"  
        maxFitRange = hist.GetXaxis().GetBinLowEdge(1+hist.GetNbinsX())
        minFitRange = hist.GetXaxis().GetBinLowEdge(1)
    else:
        #print "check"
        fitopt = "RQMFS+"
        maxFitRange = 120 #if isIso else 60
        minFitRange = 20

    ###################
    # fits
    ####################
    tf1_erf = ROOT.TF1("tf1_erf","[0]*TMath::Erf((x-[1])/[2])",minFitRange,maxFitRange) 
    tf1_erf.SetParameter(0,1.0)
    tf1_erf.SetParameter(1,38)
    if not isEle and isIso and key == 6: tf1_erf.SetParameter(1,30)
    tf1_erf.SetParameter(2,3.0)
    tf1_erf.SetLineWidth(2)
    tf1_erf.SetLineColor(ROOT.kOrange+1)

    tf1_ln = ROOT.TF1("tf1_ln","[0]*TMath::Log([2]*x-[1])",minFitRange,maxFitRange)
    tf1_ln.SetParameter(0,1)
    tf1_ln.SetParameter(1,30)
    tf1_ln.SetParameter(2,1)
    tf1_ln.SetLineWidth(2)
    tf1_ln.SetLineColor(ROOT.kBlue)

    # tf1_ln2 = ROOT.TF1("tf1_ln2","[0]*TMath::Log([2]*x-[1]) + [3] + [4]*x",minFitRange,maxFitRange)
    # tf1_ln2.SetParameter(0,1)
    # tf1_ln2.SetParameter(1,30)
    # tf1_ln2.SetParameter(2,1)
    # tf1_ln2.SetParameter(3,0)    
    # tf1_ln2.SetParameter(4,0); tf1_ln2.SetParLimits(4,0, 1e12) 
    # tf1_ln2.SetLineWidth(2)
    # tf1_ln2.SetLineColor(ROOT.kOrange+1)

    tf1_pol1 = ROOT.TF1("tf1_pol1","pol1",minFitRange,maxFitRange)
    tf1_pol1.SetLineWidth(2)
    tf1_pol1.SetLineColor(ROOT.kBlue)

    tf1_pol2 = ROOT.TF1("tf1_pol2","pol2",minFitRange,maxFitRange)
    tf1_pol2.SetLineWidth(2)
    tf1_pol2.SetLineColor(ROOT.kGreen+2)

    tf1_pol3 = ROOT.TF1("tf1_pol3","pol3",minFitRange,maxFitRange)
    tf1_pol3.SetLineWidth(2)
    tf1_pol3.SetLineColor(ROOT.kCyan+1)

    # tf1_sqrt = ROOT.TF1("tf1_sqrt","[0]*TMath::Sqrt([2]*x-[1])",minFitRange,maxFitRange)
    # tf1_sqrt.SetParameter(0,1)
    # tf1_sqrt.SetParameter(1,30)
    # tf1_sqrt.SetParameter(2,1)
    # tf1_sqrt.SetLineWidth(2)
    # tf1_sqrt.SetLineColor(ROOT.kOrange+2)

    # tf1_exp = ROOT.TF1("tf1_exp","[0] - [3]*TMath::Exp([1]*x-[2])",minFitRange,maxFitRange)
    # tf1_exp.SetParameter(0,1)
    # tf1_exp.SetParameter(1,-5);    tf1_exp.SetParLimits(1,-1e12,0.0) # negative
    # tf1_exp.SetParameter(2,30)
    # tf1_exp.SetParameter(3,1); #    tf1_exp.SetParLimits(3,0,10) # positive
    # tf1_exp.SetLineWidth(2)
    # tf1_exp.SetLineColor(ROOT.kOrange+1)

    tf1_erf2 = ROOT.TF1("tf1_erf2","[0]*TMath::Erf((x-[1])/[2]) + [4] + [3]*x",minFitRange,maxFitRange) 
    tf1_erf2.SetParameter(0,1.0); #tf1_erf2.SetParLimits(0,0.1,1e12)
    tf1_erf2.SetParameter(1,32)
    tf1_erf2.SetParameter(2,3.0)
    if isTrigger:
        if isEle:
            if mc =="Data":
                if key == 29:
                    tf1_erf2.SetParameter(1,31)
            elif mc == "MC":
                if any(key == x for x in [20,29]):
                    tf1_erf2.SetParameter(1,31)                
                elif key == 12:
                    tf1_erf2.SetParameter(1,32)
                elif key == 17:
                    tf1_erf2.SetParameter(1,31.2)
        else:
            tf1_erf2.SetParameter(1,30)
            if mc == "Data":
                if key == 27:
                    tf1_erf2.SetParameter(1,32)
                # elif key == 13:
                #     tf1_erf2.SetParameter(1,29.5)
                elif any(key == x for x in [11,13,14,26]):
                    tf1_erf2.SetParameter(1,29.5)
                elif any(key == x for x in [10,23]):
                    tf1_erf2.SetParameter(1,30.5)
                    tf1_erf2.SetParameter(2,4.0)
            if mc == "MC":
                if any(key == x for x in [0,8,26,27,37]):
                    tf1_erf2.SetParameter(1,29.5)
                elif key == 33:
                    tf1_erf2.SetParameter(1,27)
            
    tf1_erf2.SetParameter(3,0.0); #tf1_erf2.SetParLimits(3,0,1e12)
    tf1_erf2.SetParameter(4,0)
    tf1_erf2.SetLineWidth(2)
    tf1_erf2.SetLineColor(ROOT.kRed+1)

    # fit and draw (if required)
    if isTrigger or isFullID:
        hist.Fit(tf1_erf,fitopt)        
        # hist.Fit(tf1_ln2,fitopt)        
        hist.Fit(tf1_ln,fitopt)
        if isFullID: 
            tf1_pol3.SetLineColor(ROOT.kRed+1)
        else:
            hist.Fit(tf1_erf2,fitopt)        
        # hist.Fit(tf1_sqrt,fitopt)        
        # hist.Fit(tf1_exp,fitopt)        
        hist.Fit(tf1_pol2,fitopt)        
        hist.Fit(tf1_pol3,fitopt)        
    else:
        if isIso: hist.Fit(tf1_erf,fitopt)        
        else:
            tf1_pol3.SetLineColor(ROOT.kRed+1)
            hist.Fit(tf1_pol1,fitopt)        
            hist.Fit(tf1_pol2,fitopt)        
            hist.Fit(tf1_pol3,fitopt)        

    leg = ROOT.TLegend(0.5, 0.2, 0.9, 0.45 if isEle else 0.3 if isIso else 0.4)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)

    legEntry = {}
    legEntry[tf1_erf2.GetName()] = "Erf[x] + ax + b"
    legEntry[tf1_erf.GetName()]  = "Erf[x]"
    legEntry[tf1_pol2.GetName()] = "pol2"
    legEntry[tf1_pol3.GetName()] = "pol3"
    legEntry[tf1_ln.GetName()]   = "a ln(bx + c)"
    legEntry[tf1_pol1.GetName()] = "pol1"

    if isTrigger or isFullID:
        if not isFullID: leg.AddEntry(tf1_erf2, "Erf[x] + ax + b", 'LF')
        leg.AddEntry(tf1_ln,  "a ln(bx + c)", "LF")
        #leg.AddEntry(tf1_ln2, "a ln(bx + c) + dx + e", "LF")
        #leg.AddEntry(tf1_sqrt,"a #sqrt{bx +c} +cx", "LF")
        leg.AddEntry(tf1_erf, "Erf[x]", 'LF')
        #leg.AddEntry(tf1_exp, "a - exp(-bx -c) +dx", "LF")
        leg.AddEntry(tf1_pol2,"pol2", "LF")
        leg.AddEntry(tf1_pol3,"pol3", "LF")
    else:
        if isIso: leg.AddEntry(tf1_erf, "Erf[x]", 'LF')
        else:
            leg.AddEntry(tf1_pol1,"pol1", "LF")
            leg.AddEntry(tf1_pol2,"pol2", "LF")
            leg.AddEntry(tf1_pol3,"pol3", "LF")
    leg.Draw('same')
    ###################
    # fits
    ####################

    canvas.RedrawAxis("sameaxis")

    setTDRStyle()
    ROOT.gStyle.SetOptTitle(1)  # use histogram title with binning as canvas title

    # for ext in ["pdf","png"]:
    #     if mc == "SF":
    #         canvas.SaveAs("{out}ScaleFactorVsPt_{mc}_{ch}_eta{b}.{ext}".format(out=outdir,mc=mc,ch=channel,b=key,ext=ext))            
    #     else:
    #         canvas.SaveAs("{out}effVsPt_{mc}_{ch}_eta{b}.{ext}".format(out=outdir,mc=mc,ch=channel,b=key,ext=ext))            

    if isTrigger or isFullID:
        fit_pol2 = hist.GetFunction(tf1_pol2.GetName())
        fit_pol3 = hist.GetFunction(tf1_pol3.GetName())
        fit_erf =  hist.GetFunction(tf1_erf.GetName())
        if not isFullID: fit_erf2 =  hist.GetFunction(tf1_erf2.GetName())
        fit_ln =   hist.GetFunction(tf1_ln.GetName())
        #fit_ln2 = hist.GetFunction(tf1_ln2.GetName())
        #fit_sqrt = hist.GetFunction(tf1_sqrt.GetName())
        #fit_exp = hist.GetFunction(tf1_exp.GetName())
    else:
        if isIso: fit_erf = hist.GetFunction(tf1_erf.GetName())
        else:
            fit_pol1 = hist.GetFunction(tf1_pol1.GetName())
            fit_pol2 = hist.GetFunction(tf1_pol2.GetName())
            fit_pol3 = hist.GetFunction(tf1_pol3.GetName())

    functions = {}
    if isTrigger or isFullID:
        functions[tf1_pol2.GetName()] = fit_pol2
        functions[tf1_pol3.GetName()] = fit_pol3
        functions[tf1_erf.GetName()] = fit_erf
        functions[tf1_ln.GetName()] = fit_ln
        #functions[tf1_ln2.GetName()] = fit_ln2
        #functions[tf1_sqrt.GetName()] = fit_sqrt
        if not isFullID: functions[tf1_erf2.GetName()] = fit_erf2
        #functions[tf1_exp.GetName()] = fit_exp
    else:
        if isIso: functions[tf1_erf.GetName()] = fit_erf
        else:
            functions[tf1_pol1.GetName()] = fit_pol1
            functions[tf1_pol2.GetName()] = fit_pol2
            functions[tf1_pol3.GetName()] = fit_pol3

    lat = ROOT.TLatex()
    line = ""
    lat.SetNDC();
    lat.SetTextSize(0.045);
    lat.SetTextFont(42);
    lat.SetTextColor(1);
    xmin = 0.20 
    yhi = 0.85

    if isEle==False and isIso:
        if hist_chosenFunc: hist_chosenFunc.Fill(tf1_erf.GetName(),1)
        retFunc = fit_erf
    else:
        chi2 = 1000000.0
        funcMinChi2 = 0
        for name,f in functions.iteritems():        
            #print "Name: %s func %s" % (name, f) 
            if f.GetNDF() == 0: continue
            if name == tf1_pol3.GetName(): continue
            if f.GetChisquare() < chi2: 
                chi2 = f.GetChisquare()
                funcMinChi2 = f
        #print "Function %s had the best Chi2/Ndof: %.3f/%d among non-pol3" % (funcMinChi2.GetName(),funcMinChi2.GetChisquare(),funcMinChi2.GetNDF())
        #print "pol3 Chi2/Ndof: %.3f/%d" % (fit_pol3.GetChisquare(),fit_pol3.GetNDF())

        if funcMinChi2 == 0:
            print "="*20
            print "Warning: no function had more than 0 degrees of freedom. Returning pol2 function"
            print "="*20
            if hist_chosenFunc: hist_chosenFunc.Fill(fit_pol2.GetName(),1) 
            #return fit_pol3
            line = "Best fit: pol2 (forced)"
            retFunc = fit_pol2

        else:
            nChi2Sigma = abs(funcMinChi2.GetChisquare()-funcMinChi2.GetNDF())/math.sqrt(2.0*funcMinChi2.GetNDF())  # Chi2 variance is 2*Ndof
            nChi2Sigma_pol3 = abs(fit_pol3.GetChisquare()-fit_pol3.GetNDF())/math.sqrt(2.0*fit_pol3.GetNDF()) if fit_pol3.GetNDF() else 999

            # pol3 will generally fit very well also in case of weird points
            # for good looking points, pol3 might be better because it can change curvature, while other functions cannot (which would be more physical)
            # allow non-pol3 fit to have Chi2 within 2 standard deviation from the expected one
            # in this case choose that value, otherwise use the one closer to expected Chisquare
            if nChi2Sigma < 3:  
                if hist_chosenFunc: hist_chosenFunc.Fill(funcMinChi2.GetName(),1)
                retFunc = funcMinChi2
            elif nChi2Sigma_pol3 < nChi2Sigma:
                if hist_chosenFunc: hist_chosenFunc.Fill(fit_pol3.GetName(),1)
                retFunc = fit_pol3
            else:
                if hist_chosenFunc: hist_chosenFunc.Fill(funcMinChi2.GetName(),1)
                retFunc = funcMinChi2                
            line = "Best fit: " + legEntry[retFunc.GetName()]

        #return funcMinChi2

            
    lat.DrawLatex(xmin,yhi,line);
    for ext in ["pdf","png"]:
        if mc == "SF":
            canvas.SaveAs("{out}ScaleFactorVsPt_{mc}_{ch}_eta{b}.{ext}".format(out=outdir,mc=mc,ch=channel,b=key,ext=ext))            
        else:
            canvas.SaveAs("{out}effVsPt_{mc}_{ch}_eta{b}.{ext}".format(out=outdir,mc=mc,ch=channel,b=key,ext=ext))                            
            
    return retFunc
        
        # get pol2 parameters (parameter number 0,1,2 correspond to constant term, x, x^2 respectively)
        #return fit.GetParameter(0),fit.GetParError(0),fit.GetParameter(1),fit.GetParError(1),fit.GetParameter(2),fit.GetParError(2)    
        #return fit_pol2
        
if __name__ == "__main__":
            
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-i','--input', dest='inputfile', default='', type='string', help='input root file with TH2. Not needed when using --make-weighted-average')
    parser.add_option('-o','--outdir', dest='outdir', default='', type='string', help='output directory to save things')
    parser.add_option('-n','--outfilename', dest='outfilename', default='', type='string', help='Name of output file to save fit results')
    parser.add_option('-c','--channel', dest='channel', default='', type='string', help='name of the channel (mu or el)')
    parser.add_option('-e','--era',     dest='era',     default='', type='string', help='For muons: select data era GtoH or BtoF as -e GH or -e BF')
    parser.add_option('-v','--var',     dest='variable',default='', type='string', help='For muons: select variable: ISO or ID')
    parser.add_option('-w','--width-pt',     dest='widthPt',default='0.2', type='float', help='Pt bin width for the smoothed histogram')
    parser.add_option('-t','--trigger', dest='isTriggerScaleFactor',action="store_true", default=False, help='Says if using trigger scale factors (electron and muon share the same root file content)')
    parser.add_option(     '--fullID', dest='isFullIDScaleFactor',action="store_true", default=False, help='For electrons: says if using fullID scale factor')
    parser.add_option(     '--save-TF1', dest='saveTF1',action="store_true", default=False, help='Save TF1 as well, not just TH2 with many bins (note that they are not saved when making averages between eras')
    parser.add_option(     '--make-weighted-average', dest='isWeightedAverage',action="store_true", default=False, help='To be used if you are averaging the scale factors (must use other options to pass the inputs. For example, muons have two sets of ID and ISO scale factors for different eras')
    parser.add_option(    '--files-average', dest='filesAverage',default='', type='string', help='Comma separated list of two files')
    parser.add_option(    '--weights-average', dest='weightsAverage',default='', type='string', help='Comma separated list of two weights (can be luminosity for era BF and GH, they are 19.91/fb and 16.30/fb respectively)')
    parser.add_option(    '--average-uncertainty-mode', dest='averageUncertaintyMode',default='max', type='string', help='Can be max, weight, diff (see function makeSFweightedAverage)')
    parser.add_option(    '--hists-average', dest='histsAverage',default='', type='string', help='Comma separated list of histograms to average (pass names, one for each average, it is assumed the names in each file are the same)')
    (options, args) = parser.parse_args()


    ROOT.TH1.SetDefaultSumw2()

    if options.isTriggerScaleFactor:
        if options.era or options.variable:
            print "Error: option -t is incompatible with -e and -v. Exit"
            quit()
        if options.isFullIDScaleFactor:
            print "Error: option -t is incompatible with --fullID. Exit"
            quit()
    
    channel = options.channel
    if channel not in ["el","mu"]:
        print "Error: unknown channel %s (select 'el' or 'mu')" % channel
        quit()
    isEle = True if channel == "el" else False
    lepton = "electron" if channel == "el" else "muon"

    if not isEle and not options.isTriggerScaleFactor:
        print "Warning: you didn't use option -t, so I assume these are not trigger scale factors"
        if not options.isWeightedAverage:
            if options.era not in ["BF","GH"]:
                print "Error: you should specify a data range for muons using option -e BF|GH. Exit"                                
                quit()                                                          
        elif options.era:
            print "Error: option --make-weighted-average is incompatible with -e. Exit"                                
            quit()                                                          

        if options.variable not in ["ID","ISO",]:
            print "Error: you should specify a variable with option -v ID|ISO. Exit"                                
            quit()                                    

    if not isEle and options.isFullIDScaleFactor:
        print "Error: option --fullID is only for electrons. Exit"
        quit()



    if options.outdir:
        outname = options.outdir
        addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
        createPlotDirAndCopyPhp(outname)
    else:
        print "Error: you should specify an output folder using option -o <name>. Exit"
        quit()

    if not options.outfilename:
        print "Error: you should specify an output file name using option -n <name>. Exit"
        quit()
    outfilename = ".".join(options.outfilename.split('.')[:-1]) + "_{lep}".format(lep="electrons" if channel=="el" else "muons")
    if options.era: outfilename = outfilename + "_" + options.era
    if options.isWeightedAverage: outfilename = outfilename + "_full2016"
    if options.variable: outfilename = outfilename + "_" + options.variable
    if options.isTriggerScaleFactor: outfilename = outfilename + "_trigger"
    if options.isFullIDScaleFactor: outfilename = outfilename + "_fullID"
    outfilename += ".root"

    #########################################    
    #########################################    

    if options.isWeightedAverage:

        hists = options.histsAverage.split(",")
        f1,f2 = options.filesAverage.split(",")
        w1,w2 = list(float(x) for x in options.weightsAverage.split(","))
        hlist = {}
        for h in hists:
            hScaleFactorAverage = makeSFweightedAverage(f1,f2,w1,w2,str(h),"new_"+str(h),uncertaintyRule=options.averageUncertaintyMode)
            hScaleFactorAverage.SetDirectory(0)
            hlist[str(h)] = hScaleFactorAverage
        tf = ROOT.TFile.Open(outname+outfilename,'recreate')
        for k in hlist:
            hlist[k].Write(k)            
            zaxisName = hlist[k].GetZaxis().GetTitle()
            if options.variable == "ID" and "scaleFactor" in k:
                zaxisName += "::0.98,1.01"
            drawCorrelationPlot(hlist[k],"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),zaxisName,
                                "combinedFull2016_{n}".format(n=k),"ForceTitle",outname,1,1,False,False,False,1,palette=55)
        tf.Close()
        print ""
        print "Created file %s" % (outname+outfilename)
        print ""
        quit()
    #########################################    
    #########################################    

    hmc = 0
    hdata = 0
    hsf = 0
    if options.inputfile:
        tf = ROOT.TFile.Open(options.inputfile)        
        if options.isTriggerScaleFactor or options.isFullIDScaleFactor:
            hmc =   tf.Get("EGamma_EffMC2D")
            hdata = tf.Get("EGamma_EffData2D")
            hsf = tf.Get("EGamma_SF2D")
            if (hsf == 0):
                print "Error: could not retrieve hsf from input file %s. Exit" % options.inputfile
                quit()
        elif not isEle:            
            hmc   = tf.Get("eff%s_mc_%s" % (options.variable,options.era))
            hdata = tf.Get("eff%s_data_%s" % (options.variable,options.era))            
        else:
            print "Error: you are doing electrons, but you didn't specify option -t for trigger or --fullID for full ID. This setup is currently not implemented. Exit"
            quit()

        if (hmc == 0 or hdata == 0):
            print "Error: could not retrieve hdata or hmc from input file %s. Exit" % options.inputfile
            quit()
        else:
            hmc.SetDirectory(0)
            hdata.SetDirectory(0)
            if isEle or options.isTriggerScaleFactor: hsf.SetDirectory(0)
        tf.Close()
    else:
        print "Error: you should specify an input file using option -i <name>. Exit"
        quit()
        

    etabins = hdata.GetXaxis().GetXbins()
    ptbins = hdata.GetYaxis().GetXbins()
    #etaBinHisto = ROOT.TH1F("etaBinEdges","The x axis of this histogram has the eta binning",len(etabins)-1,array('d',etabins))
        
    # for muons must create original scale factor as well, unless it was trigger (which also have the SF in the input file)
    if isEle or options.isTriggerScaleFactor:
        pass
    else:
        hsf = ROOT.TH2D("hsf","",
                        len(etabins)-1,array('d',etabins),
                        len(ptbins)-1,array('d',ptbins)
                        )
        copyHisto(hsf,hdata)  
        hsf.Divide(hmc)
    

    hist_chosenFunc = ROOT.TH1D("chosenFitFunc","Best fit function for each eta bin",5,0,5)
    hist_chosenFunc.GetXaxis().SetBinLabel(1,"tf1_erf")
    hist_chosenFunc.GetXaxis().SetBinLabel(2,"tf1_erf2")
    if isEle or options.isTriggerScaleFactor:        
        hist_chosenFunc.GetXaxis().SetBinLabel(3,"tf1_ln")
    else:
        hist_chosenFunc.GetXaxis().SetBinLabel(3,"tf1_pol1")
    hist_chosenFunc.GetXaxis().SetBinLabel(4,"tf1_pol2")
    hist_chosenFunc.GetXaxis().SetBinLabel(5,"tf1_pol3")
    ######################
    # to make ratio
    ######################
    ratioData =  ROOT.TH2D("dataEfficiencyRatio","Original/smooth Data efficiency ratio",
                           len(etabins)-1,array('d',etabins),
                           len(ptbins)-1,array('d',ptbins)
                           )
    ratioMC =  ROOT.TH2D("mcEfficiencyRatio","Original/smooth MC efficiency ratio",
                           len(etabins)-1,array('d',etabins),
                           len(ptbins)-1,array('d',ptbins)
                           )
    ratioSF =  ROOT.TH2D("scaleFactorRatio","Original/smooth scale factor ratio",
                         len(etabins)-1,array('d',etabins),
                         len(ptbins)-1,array('d',ptbins)
                         )
    copyHisto(ratioData,hdata)
    copyHisto(ratioMC,hmc)
    copyHisto(ratioSF,hsf)

    ######################
    # These histograms will contain the parameters of the function used to smooth the efficiencies
    ######################
    # currently not used anymore

    # hdataSmoothEff = ROOT.TH2D("hdataSmoothEff","Data efficiency: fit parameters --> a_{0} + a_{1}*x + a_{2}*x^{2}",
    #                            len(etabins)-1,array('d',etabins),
    #                            3,0.5,3.5
    #                            )
    # hmcSmoothEff = ROOT.TH2D("hmcSmoothEff","MC efficiency: fit parameters --> a_{0} + a_{1}*x + a_{2}*x^{2}",
    #                          len(etabins)-1,array('d',etabins),
    #                          3,0.5,3.5
    #                          )

    #############
    # these will be used to check the smoothed efficiency
    ###############
    hdataSmoothCheck = ROOT.TH2D("hdataSmoothCheck","Data smoothed efficiency",
                                 len(etabins)-1,array('d',etabins),
                                 int(math.ceil(hdata.GetYaxis().GetBinLowEdge(1+hdata.GetNbinsY()) - hdata.GetYaxis().GetBinLowEdge(1))/options.widthPt), # bins of 0.2 GeV
                                 hdata.GetYaxis().GetBinLowEdge(1),hdata.GetYaxis().GetBinLowEdge(1+hdata.GetNbinsY())
                                 )
    hmcSmoothCheck = ROOT.TH2D("hmcSmoothCheck","MC smoothed efficiency",
                               len(etabins)-1,array('d',etabins),
                               int(math.ceil(hmc.GetYaxis().GetBinLowEdge(1+hmc.GetNbinsY()) - hmc.GetYaxis().GetBinLowEdge(1))/options.widthPt),
                               hmc.GetYaxis().GetBinLowEdge(1),hmc.GetYaxis().GetBinLowEdge(1+hmc.GetNbinsY())
                               )
    hsfSmoothCheck = ROOT.TH2D("hsfSmoothCheck","Data/MC smoothed scale factor",
                               len(etabins)-1,array('d',etabins),
                               int(math.ceil(hdata.GetYaxis().GetBinLowEdge(1+hdata.GetNbinsY()) - hdata.GetYaxis().GetBinLowEdge(1))/options.widthPt),
                               hdata.GetYaxis().GetBinLowEdge(1),hdata.GetYaxis().GetBinLowEdge(1+hdata.GetNbinsY())
                               )
    hdataSmoothCheck_origBinPt = ROOT.TH2D("hdataSmoothCheck_origBinPt","Data smoothed efficiency",
                                           len(etabins)-1,array('d',etabins),
                                           len(ptbins)-1,array('d',ptbins)
                                           )
    hmcSmoothCheck_origBinPt = ROOT.TH2D("hmcSmoothCheck_origBinPt","MC smoothed efficiency",
                                         len(etabins)-1,array('d',etabins),
                                         len(ptbins)-1,array('d',ptbins)
                                         )
    hsfSmoothCheck_origBinPt = ROOT.TH2D("hsfSmoothCheck_origBinPt","Data/MC smoothed scale factor",
                                         len(etabins)-1,array('d',etabins),
                                         len(ptbins)-1,array('d',ptbins)
                                         )

    # hmc and hdata have eta on X and pt on Y
    # we select slices at constant eta and fit along pt with some function
    # let's select an error function 

    hmcpt = {}
    for x in range(1,hmc.GetNbinsX()+1):
        bin = x-1
        hmcpt[bin] = ROOT.TH1D("hmcpt_{b}".format(b=str(bin)),
                               "MC: %.4g <= #eta < %.4g" % (hmc.GetXaxis().GetBinLowEdge(x), hmc.GetXaxis().GetBinLowEdge(x+1)),
                               len(ptbins)-1,array('d',ptbins)
                               )
        for y in range(1,hmc.GetNbinsY()+1):
            hmcpt[bin].SetBinContent(y,hmc.GetBinContent(x,y))         
            #hmcpt[bin].SetBinError(y,hmc.GetBinError(x,y))  # no error on MC, have to assign a random value to fit
            hmcpt[bin].SetBinError(y,0.5*hdata.GetBinError(x,y))  # use half the corresponding uncertainty on data (any value which is smaller would do)

    hdatapt = {}
    for x in range(1,hdata.GetNbinsX()+1):
        bin = x-1
        hdatapt[bin] = ROOT.TH1D("hdatapt_{b}".format(b=str(bin)),
                                 "Data: %.4g <= #eta < %.4g" % (hdata.GetXaxis().GetBinLowEdge(x), hdata.GetXaxis().GetBinLowEdge(x+1)),
                                 len(ptbins)-1,array('d',ptbins)
                                 )
        for y in range(1,hdata.GetNbinsY()+1):
            hdatapt[bin].SetBinContent(y,hdata.GetBinContent(x,y))         
            hdatapt[bin].SetBinError(y,hdata.GetBinError(x,y))
            #hdatapt[bin].SetBinError(y,0.01)

    hsfpt = {}
    for x in range(1,hsf.GetNbinsX()+1):
        bin = x-1
        hsfpt[bin] = ROOT.TH1D("hsfpt_{b}".format(b=str(bin)),
                               "Data/MC: %.4g <= #eta < %.4g" % (hsf.GetXaxis().GetBinLowEdge(x), hsf.GetXaxis().GetBinLowEdge(x+1)),
                               len(ptbins)-1,array('d',ptbins)
                               )
        for y in range(1,hsf.GetNbinsY()+1):
            hsfpt[bin].SetBinContent(y,hsf.GetBinContent(x,y))         
            hsfpt[bin].SetBinError(y,hsf.GetBinError(x,y)) 

    ###########################
    # first MC
    ###########################
    bestFit_MC = {}
    for key in hmcpt:

        #fitpol2 = fitTurnOn(hmcpt[key],key,outname, "MC",channel=channel,hist_chosenFunc=hist_chosenFunc,isTrigger=options.isTriggerScaleFactor)
        # for pol2 only
        # for ipar in range(3):
        #     hmcSmoothEff.SetBinContent(key+1,fitpol2.GetParameter(ipar))
        #     hmcSmoothEff.SetBinError(key+1,fitpol2.GetParError(ipar))
        # a0,a1,a2 = fitpol2.GetParameter(0),fitpol2.GetParameter(1),fitpol2.GetParameter(2)
        bestFitFunc = fitTurnOn(hmcpt[key],key,outname, "MC",channel=channel,hist_chosenFunc=hist_chosenFunc, 
                                isIso=True if options.variable=="ISO" else False,
                                isTrigger=options.isTriggerScaleFactor,
                                isFullID=options.isFullIDScaleFactor)
        bestFit_MC["smoothFunc_MC_ieta%d" % key] = bestFitFunc
        for ipt in range(1,hmcSmoothCheck.GetNbinsY()+1):
            ptval = hmcSmoothCheck.GetYaxis().GetBinCenter(ipt)
            #hmcSmoothCheck.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
            hmcSmoothCheck.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))
        for ipt in range(1,hmcSmoothCheck_origBinPt.GetNbinsY()+1):
            ptval = hmcSmoothCheck_origBinPt.GetYaxis().GetBinCenter(ipt)
            #hmcSmoothCheck_origBinPt.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
            hmcSmoothCheck_origBinPt.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))

    ###########################
    # now data
    ###########################
    bestFit_Data = {}
    for key in hdatapt:

        # fitpol2 = fitTurnOn(hdatapt[key],key,outname,"Data",channel=channel,hist_chosenFunc=hist_chosenFunc)
        # for ipar in range(3):
        #     hdataSmoothEff.SetBinContent(key+1,fitpol2.GetParameter(ipar))
        #     hdataSmoothEff.SetBinError(key+1,fitpol2.GetParError(ipar))
        # a0,a1,a2 = fitpol2.GetParameter(0),fitpol2.GetParameter(1),fitpol2.GetParameter(2)
        # for ipt in range(1,hdataSmoothCheck.GetNbinsY()+1):
        #     ptval = hdataSmoothCheck.GetYaxis().GetBinCenter(ipt)
        #     hdataSmoothCheck.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
        # for ipt in range(1,hdataSmoothCheck_origBinPt.GetNbinsY()+1):
        #     ptval = hdataSmoothCheck_origBinPt.GetYaxis().GetBinCenter(ipt)
        #     hdataSmoothCheck_origBinPt.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)

        bestFitFunc = fitTurnOn(hdatapt[key],key,outname, "Data",channel=channel,hist_chosenFunc=hist_chosenFunc, 
                                isIso=True if options.variable=="ISO" else False,
                                isTrigger=options.isTriggerScaleFactor,
                                isFullID=options.isFullIDScaleFactor)
        bestFit_Data["smoothFunc_Data_ieta%d" % key] = bestFitFunc
        for ipt in range(1,hdataSmoothCheck.GetNbinsY()+1):
            ptval = hdataSmoothCheck.GetYaxis().GetBinCenter(ipt)
            #hdataSmoothCheck.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
            hdataSmoothCheck.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))
        for ipt in range(1,hdataSmoothCheck_origBinPt.GetNbinsY()+1):
            ptval = hdataSmoothCheck_origBinPt.GetYaxis().GetBinCenter(ipt)
            #hdataSmoothCheck_origBinPt.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
            hdataSmoothCheck_origBinPt.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))

    ###########################
    # now data/MC
    ###########################
    bestFit_SF = {}
    for key in hsfpt:

        # fitpol2 = fitTurnOn(hsfpt[key],key,outname,"SF",channel=channel,hist_chosenFunc=hist_chosenFunc)
        # for ipar in range(3):
        #     hsfSmoothEff.SetBinContent(key+1,fitpol2.GetParameter(ipar))
        #     hsfSmoothEff.SetBinError(key+1,fitpol2.GetParError(ipar))
        # a0,a1,a2 = fitpol2.GetParameter(0),fitpol2.GetParameter(1),fitpol2.GetParameter(2)
        # for ipt in range(1,hsfSmoothCheck.GetNbinsY()+1):
        #     ptval = hsfSmoothCheck.GetYaxis().GetBinCenter(ipt)
        #     hsfSmoothCheck.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
        # for ipt in range(1,hsfSmoothCheck_origBinPt.GetNbinsY()+1):
        #     ptval = hsfSmoothCheck_origBinPt.GetYaxis().GetBinCenter(ipt)
        #     hsfSmoothCheck_origBinPt.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)

        # do not fill control histogram when fitting scale factors (these fits are not really used, sf are actually obtained dividing the fitted efficiencies)
        bestFitFunc = fitTurnOn(hsfpt[key],key,outname, "SF",channel=channel,hist_chosenFunc=0, 
                                isIso=True if options.variable=="ISO" else False,
                                isTrigger=options.isTriggerScaleFactor,
                                isFullID=options.isFullIDScaleFactor)
        bestFit_SF["smoothFunc_SF_ieta%d" % key] = bestFitFunc
        for ipt in range(1,hsfSmoothCheck.GetNbinsY()+1):
            ptval = hsfSmoothCheck.GetYaxis().GetBinCenter(ipt)
            #hsfSmoothCheck.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
            hsfSmoothCheck.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))
        for ipt in range(1,hsfSmoothCheck_origBinPt.GetNbinsY()+1):
            ptval = hsfSmoothCheck_origBinPt.GetYaxis().GetBinCenter(ipt)
            #hsfSmoothCheck_origBinPt.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
            hsfSmoothCheck_origBinPt.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))

    #################################
    # start to make plots
    #################################
    zaxisRange = "0.4,1.1" if isEle else "0.96,1.01"
    if not isEle and options.isTriggerScaleFactor: zaxisRange = "0.8,1.01"
    if options.variable == "ISO":
        zaxisRange = "0.85,1.01"

    # for muons, plot also official scale factors
    if not isEle and not options.isTriggerScaleFactor:        
        tf = ROOT.TFile.Open("official_muon_ScaleFactors_2016/Run{era}_SF_{var}.root".format(era=options.era,var=options.variable))        
        if options.variable == "ISO":
            hsf_official = tf.Get("NUM_TightRelIso_DEN_MediumID_eta_pt")
        else:
            hsf_official = tf.Get("NUM_MediumID_DEN_genTracks_eta_pt")
        if (hsf_official == 0):
            print "Error: could not retrieve hsf_official from input file %s. Exit" % options.inputfile
            quit()
        hsf_official.SetDirectory(0)
        tf.Close()
        drawCorrelationPlot(hsf_official,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Official data/MC scale factor::%s" % zaxisRange,
                            "input_official_scaleFactor,","",outname,1,1,False,False,False,1,palette=55)

    # plot original histograms

    drawCorrelationPlot(hmc,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC efficiency::%s" % zaxisRange,
                        "inputEfficiency_MC","",outname,1,1,False,False,False,1,palette=55)
    drawCorrelationPlot(hdata,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data efficiency::%s" % zaxisRange,
                        "inputEfficiency_Data","",outname,1,1,False,False,False,1,palette=55)
    drawCorrelationPlot(hsf,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC scale factor::%s" % zaxisRange,
                        "inputScaleFactor","",outname,1,1,False,False,False,1,palette=55)

    # get a smoothed version of those input histograms
    drawCorrelationPlot(hmc,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC efficiency::%s" % zaxisRange,
                        "inputEfficiency_MC_smooth","",outname,1,1,True,False,False,1,palette=55)
    drawCorrelationPlot(hdata,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data efficiency::%s" % zaxisRange,
                        "inputEfficiency_Data_smooth","",outname,1,1,True,False,False,1,palette=55)
    drawCorrelationPlot(hsf,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC scale factor::%s" % zaxisRange,
                        "inputScaleFactor_smooth","",outname,1,1,True,False,False,1,palette=55)


    # now the new ones

    # make a sanity check plot: fill eta-pt with smoothed efficiency
    drawCorrelationPlot(hmcSmoothCheck,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC smoothed efficiency::%s" % zaxisRange,
                        "smoothEfficiency_MC","ForceTitle",outname,1,1,False,False,False,1,palette=55)
    drawCorrelationPlot(hdataSmoothCheck,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data smoothed efficiency::%s" % zaxisRange,
                        "smoothEfficiency_Data","ForceTitle",outname,1,1,False,False,False,1,palette=55)
    drawCorrelationPlot(hsfSmoothCheck,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC smoothed scale factor::%s" % zaxisRange,
                        "smoothScaleFactorDirectly","ForceTitle",outname,1,1,False,False,False,1,palette=55)

    # scale factor: data/MC
    scaleFactor = ROOT.TH2D("scaleFactor","Scale factor",
                            len(etabins)-1,array('d',etabins),
                            int(math.ceil((hdata.GetYaxis().GetBinLowEdge(1+hdata.GetNbinsY()) - hdata.GetYaxis().GetBinLowEdge(1))/options.widthPt)),
                            hdata.GetYaxis().GetBinLowEdge(1),hdata.GetYaxis().GetBinLowEdge(1+hdata.GetNbinsY())
                            )

    copyHisto(scaleFactor, hdataSmoothCheck)
    scaleFactor.Divide(hmcSmoothCheck)
    scaleFactor.SetMinimum(scaleFactor.GetBinContent(scaleFactor.GetMinimumBin()))
    scaleFactor.SetMaximum(scaleFactor.GetBinContent(scaleFactor.GetMaximumBin()))
    drawCorrelationPlot(scaleFactor,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC scale factor::%s" % zaxisRange,
                        "smoothScaleFactor","ForceTitle",outname,1,1,False,False,False,1,palette=55)

    #################################
    # plot also with oiginal binning
    ################################
    
    # divide before drawing the denominator, whose axis settings are modified by drawCorrelationPlot and seems to affect the ratio as well if divided afterwards
    ratioData.Divide(hdataSmoothCheck_origBinPt)
    ratioMC.Divide(hmcSmoothCheck_origBinPt)

    drawCorrelationPlot(hmcSmoothCheck_origBinPt,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC smoothed efficiency::%s" % zaxisRange,
                        "smoothEfficiency_MC_origBinPt","ForceTitle",outname,1,1,False,False,False,1,palette=55)
    drawCorrelationPlot(hdataSmoothCheck_origBinPt,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data smoothed efficiency::%s" % zaxisRange,
                        "smoothEfficiency_Data_origBinPt","ForceTitle",outname,1,1,False,False,False,1,palette=55)
    drawCorrelationPlot(hsfSmoothCheck_origBinPt,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC smoothed scale factor::%s" % zaxisRange,
                        "smoothScaleFactorDirectly_origBinPt","ForceTitle",outname,1,1,False,False,False,1,palette=55)

    # scale factor: data/MC
    scaleFactor_origBinPt = ROOT.TH2D("scaleFactor_origBinPt","Scale factor",
                                      len(etabins)-1,array('d',etabins),
                                      len(ptbins)-1,array('d',ptbins)
                                      )
    copyHisto(scaleFactor_origBinPt, hdataSmoothCheck_origBinPt)
    scaleFactor_origBinPt.Divide(hmcSmoothCheck_origBinPt)
    scaleFactor_origBinPt.SetMinimum(scaleFactor_origBinPt.GetBinContent(scaleFactor_origBinPt.GetMinimumBin()))
    scaleFactor_origBinPt.SetMaximum(scaleFactor_origBinPt.GetBinContent(scaleFactor_origBinPt.GetMaximumBin()))

    # to make ratio, divide before passing to function, to avoid changes in the histogram
    ratioSF.Divide(scaleFactor_origBinPt)

    #scaleFactor_origBinPt.GetZaxis().SetTitle("Data/MC scale factor")
    drawCorrelationPlot(scaleFactor_origBinPt,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC scale factor::%s" % zaxisRange,
                        "smoothScaleFactor_origBinPt","ForceTitle",outname,1,1,False,False,False,1,palette=55)


    ######################
    # finally SF(smooth)/SF(original)
    ######################
    drawCorrelationPlot(ratioData,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data efficiency ratio (original/smooth)::0.98,1.02",
                        "dataEfficiencyRatio","ForceTitle",outname,1,1,False,False,False,1,palette=55)
    drawCorrelationPlot(ratioMC,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC efficiency ratio (original/smooth)::0.98,1.02",
                        "mcEfficiencyRatio","ForceTitle",outname,1,1,False,False,False,1,palette=55)
    drawCorrelationPlot(ratioSF,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"scale factor ratio (original/smooth)::0.98,1.02",
                        "scaleFactorRatio","ForceTitle",outname,1,1,False,False,False,1,palette=55)


    ######################
    # See the difference between smoothing Data and MC efficiency and taking the ratio or smoothing directly the efifciency ratio
    ######################    
    ratioSF_smoothNumDen_smoothRatio = ROOT.TH2D("ratioSF_smoothNumDen_smoothRatio","SF ratio: smooth eff or ratio directly",
                                                 len(etabins)-1,array('d',etabins),
                                                 len(ptbins)-1,array('d',ptbins)
                                                 )

    copyHisto(ratioSF_smoothNumDen_smoothRatio,scaleFactor_origBinPt)
    ratioSF_smoothNumDen_smoothRatio.Divide(hsfSmoothCheck_origBinPt)
    drawCorrelationPlot(ratioSF_smoothNumDen_smoothRatio,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),
                        "SF ratio: smooth eff or ratio directly::0.98,1.02",
                        "ratioSF_smoothNumDen_smoothRatio","ForceTitle",outname,1,1,False,False,False,1,palette=55)

    c = ROOT.TCanvas("c","",700,700)
    c.SetTickx(1)
    c.SetTicky(1)
    c.cd()
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.06)
    c.cd()
    hist_chosenFunc.Scale(1./hist_chosenFunc.Integral())
    hist_chosenFunc.GetXaxis().SetTitleOffset(1.2)
    hist_chosenFunc.GetXaxis().SetTitleSize(0.05)
    hist_chosenFunc.GetXaxis().SetLabelSize(0.06)
    #hist_chosenFunc.GetXaxis().LabelsOption("v")
    hist_chosenFunc.GetYaxis().SetTitle("Fraction of events")
    hist_chosenFunc.GetYaxis().SetTitleOffset(1.15)
    hist_chosenFunc.GetYaxis().SetTitleSize(0.05)
    hist_chosenFunc.GetYaxis().SetLabelSize(0.04)
    hist_chosenFunc.Draw("HE")
    for ext in ["png","pdf"]:
        c.SaveAs("{out}bestFitFunction{ch}.{ext}".format(out=outname,ch=channel,ext=ext))

    # before saving things, assign an uncertainty from the original input (will need to devise a better way to estimate them)
    # following histograms have same eta-pt binning
    for ix in range(1,hdataSmoothCheck.GetNbinsX()+1):
        for iy in range(1,hdataSmoothCheck.GetNbinsY()+1):
            ieta = hdata.GetXaxis().FindFixBin(hdataSmoothCheck.GetXaxis().GetBinCenter(ix))
            ipt = hdata.GetYaxis().FindFixBin(hdataSmoothCheck.GetYaxis().GetBinCenter(ix))
            hdataSmoothCheck.SetBinError(ix,iy,hdata.GetBinError(ieta,ipt))
            hmcSmoothCheck.SetBinError(ix,iy,0.5*hdata.GetBinError(ieta,ipt))
            scaleFactor.SetBinError(ix,iy,hsf.GetBinError(ieta,ipt))

    ###########################
    # Now save things
    ###########################
    tf = ROOT.TFile.Open(outname+outfilename,'recreate')
    # hdataSmoothEff.Write()    
    # hmcSmoothEff.Write()    
    hdataSmoothCheck.Write()
    hmcSmoothCheck.Write()
    scaleFactor.Write()
    hsf.Write("scaleFactorOriginal")
    if options.saveTF1:
        for key in bestFit_MC:
            bestFit_MC[key].Write(key)
        for key in bestFit_Data:
            bestFit_Data[key].Write(key)
    tf.Close()
    print ""
    print "Created file %s" % (outname+outfilename)
    print ""

                               
         
