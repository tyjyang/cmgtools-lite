#!/usr/bin/env python3

import os, re, array, math
import time
import argparse

## safe batch mode                                 
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfolder", type=str, nargs=1)
    #parser.add_argument("outputfolder",   type=str, nargs=1)
    parser.add_argument("-w", "--working-points", dest="workingPoints", type=str, default="alwaystrue,onemuon,muonID,trigAndMatch,pfRelIso04", help="Comma separated list of working points")
    parser.add_argument("-p", "--process", default=None, required=True, type=str, help="Process to pick histogram")
    parser.add_argument("--dz", default=0.1, type=float, help="dz used for plots, in cm")
    parser.add_argument("-r", "--eff-range", dest="effRange", default=(0.5,1.05), type=float, nargs=2, help="y axis range for efficiency plot")
    parser.add_argument("-e", "--eta-ranges", dest="etaRange", default=[], type=float, nargs=2, action="append", metavar=('min','max'), help="Z axis ranges (etamin, etamax) to select, when available (remember that upper bin edge belongs to next bin in ROOT)")
    parser.add_argument("-v", "--variable", default="wpt", choices=["wpt", "mupt"], help="Variable to make efficiency as a function of it (needs to appear in histogram name inside root file)")
    parser.add_argument("--hname", default="dzGenRecoVtx", help="Root of histogram name inside root file")
    parser.add_argument("--postfix", type=str, default=None, help="Postfix for output folder")
    parser.add_argument("--print-uncertainty", dest="printUncertainty", action="store_true", help="Print also efficiency uncertainty in summary plot with selection efficiency")
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    inputfolder = args.inputfolder[0] # "plots/testNanoAOD/vertexStudy/Wboson/"

    workingPoints = [str(x) for x in args.workingPoints.split(',')]# ["alwaystrue", "onemuon", "muonID", "trigAndMatch", "pfRelIso04"]
    colors = [ROOT.kBlack, ROOT.kCyan+1, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kMagenta, ROOT.kAzure+2]
    markers = [ROOT.kFullCircle, ROOT.kFullCross, ROOT.kOpenSquare, ROOT.kOpenTriangleUp, ROOT.kOpenCircle, ROOT.kOpenSquareDiagonal, ROOT.kFullCircle, ROOT.kFullSquare]

    #range for efficiency plot
    effLow = args.effRange[0]
    effHigh = args.effRange[1]

    process = args.process
    absDzVal = 0.1 # 1 mm        
    # Wpt from 0 to 100 with 2 GeV width, dz from -1.0 to 1.0 with 0.01 width
    # preFSR lepton pT from 26 to 100 with 1 GeV width, dz from -1.0 to 1.0 with 0.01 width
    var = args.variable
    genWPtLow = 0.0
    genWPtHigh = 100.0
    rebinWPt = 1 # default binning is 2 GeV, so rebinWpt = 2 makes it 4 GeV
    genLepPtLow = 26.0
    genLepPtHigh = 100.0
    rebinLepPt = 2 # default binning is 1 GeV, so rebinLepPt = 2 makes it 2 GeV    
    hname = f"{args.hname}_{var}_{process}"

    outdir = inputfolder + f"/compareEfficiency_{process}"
    if args.postfix:
        outdir += f"_{args.postfix}"
    outdir += "/"
    
    #sortkeys = inputfiles.keys()
    #sortkeys = sorted(sortkeys, key = lambda x: int(x.split("::")[0]))
    sortkeys = workingPoints
    #print sortkeys
    #quit()

    hists = {}

    heffSel_passDz = ROOT.TH1D("heffSel_passDz","",len(workingPoints)-1,0.5,0.5+len(workingPoints)-1)
    heffSel_failDz = ROOT.TH1D("heffSel_failDz","",len(workingPoints)-1,0.5,0.5+len(workingPoints)-1)
    for b in range(1,1+heffSel_passDz.GetNbinsX()):
        heffSel_passDz.GetXaxis().SetBinLabel(b,workingPoints[b])
        heffSel_failDz.GetXaxis().SetBinLabel(b,workingPoints[b])
    
    # store integrals with no cut to compute selection efficiency as ratio of yields in each cut step
    backupIntInDz_noCuts = 0.0
    backupIntOutDz_noCuts = 0.0

    hnumtmp = ROOT.TH1D("hnumtmp","",1,0,1)
    hdentmp = ROOT.TH1D("hdentmp","",1,0,1)

    for ik,key in enumerate(sortkeys):

        print("="*30)
        print("Cut: " + str(key))
        print("-"*30)
        print
        inputfile = inputfolder + key + "/plots_vertexStudy.root"
        tf = ROOT.TFile.Open(inputfile,"READ")        
        if not tf or tf.IsZombie():
            print(f"Error opening file {inputfile}")
            quit()
        htmp = tf.Get(hname)
        print(f"file name: {inputfile}")
        checkNullObj(htmp, hname)
        htmp.SetDirectory(0)
        tf.Close()

        if htmp.GetDimension() == 3:
            h2 = htmp.Project3D("yxe") # yxe is to make TH2 with y axis versus x axis, this keeps under/overflow
            h2.SetName("histogram2D")
            if len(args.etaRange):
                h2.Reset("ICESM")
                for pair in args.etaRange:
                    binLow = htmp.GetZaxis().FindFixBin(pair[0]) 
                    binHigh = htmp.GetZaxis().FindFixBin(pair[1])
                    print(f"Projecting from range [{htmp.GetZaxis().GetBinLowEdge(binLow)}, {htmp.GetZaxis().GetBinUpEdge(binHigh)}]")
                    htmp.GetZaxis().SetRange(binLow, binHigh)
                    htmpRange = htmp.Project3D("yxe")
                    h2.Add(htmpRange)
        else:
            h2 = htmp
            
        # 0.0001 is an offset to get bin to the right of left of the edge
        dzBinLow = h2.GetYaxis().FindFixBin(-1.0 * absDzVal + 0.0001) # get bin with -absDzVal as left edge
        dzBinHigh = h2.GetYaxis().FindFixBin(absDzVal - 0.0001) # get bin with absDzVal as right edge
        dzBinsTot = h2.GetNbinsY()
        wptBinsTot = h2.GetNbinsX()
        intTot =  h2.Integral(0,1+wptBinsTot,0,1+dzBinsTot)# include under/over-flow bins
        intInDz = h2.Integral(0,1+wptBinsTot,dzBinLow,dzBinHigh)# include under/over-flow bins
        wptlow = h2.GetXaxis().GetBinLowEdge(1)
        wpthigh = h2.GetXaxis().GetBinLowEdge(wptBinsTot+1)

        if ik == 0:
            pass
        elif ik == 1:
            backupIntInDz_noCuts = intInDz
            backupIntOutDz_noCuts = intTot - intInDz
            heffSel_passDz.SetBinContent(ik,1.0)            
            heffSel_failDz.SetBinContent(ik,1.0)            
        else:
            grAsErr = ROOT.TGraphAsymmErrors()
            hnumtmp.SetBinContent(1,intInDz)
            hnumtmp.SetBinError(1,ROOT.TMath.Sqrt(intInDz))
            hdentmp.SetBinContent(1,backupIntInDz_noCuts)
            hdentmp.SetBinError(1,ROOT.TMath.Sqrt(backupIntInDz_noCuts))
            grAsErr.Divide(hnumtmp,hdentmp,"cl=0.683 b(1,1) mode")
            err = max(grAsErr.GetErrorYlow(0),grAsErr.GetErrorYhigh(0))
            heffSel_passDz.SetBinContent(ik,intInDz/backupIntInDz_noCuts)
            heffSel_passDz.SetBinError(ik,err)
            #
            hnumtmp.SetBinContent(1,intInDz)
            hnumtmp.SetBinError(1,ROOT.TMath.Sqrt(intInDz))
            hdentmp.SetBinContent(1,backupIntInDz_noCuts)
            hdentmp.SetBinError(1,ROOT.TMath.Sqrt(backupIntInDz_noCuts))
            grAsErr.Divide(hnumtmp,hdentmp,"cl=0.683 b(1,1) mode")
            err = max(grAsErr.GetErrorYlow(0),grAsErr.GetErrorYhigh(0))
            heffSel_failDz.SetBinContent(ik,(intTot-intInDz)/backupIntOutDz_noCuts)
            heffSel_failDz.SetBinError(ik,err)

        print("Inclusive efficiency: %.1f" % (100.*intInDz/intTot))


        # if ik == 5:
        #     rebinLepPt = 3 * rebinLepPt
        #     rebinWPt = 3 * rebinWPt
        #wptBins = [4.0 * i for i in range(0,26)]    
        wptBins = [genWPtLow + (rebinWPt*2.0) * i for i in range(0,1+int((genWPtHigh-genWPtLow+0.001)/(2.0*rebinWPt)))]
        if var == "mupt":
            wptBins = [genLepPtLow + (rebinLepPt*1.0) * i for i in range(0,1+int((genLepPtHigh-genLepPtLow+0.001)/rebinLepPt))]
        
        hists[key] = ROOT.TH1D("eff_"+key,"",len(wptBins)-1,array("d",wptBins))
        #gr[key] = 

        for iwpt in range(len(wptBins)-1):

            wptBinLow = h2.GetXaxis().FindFixBin(wptBins[iwpt] + 0.0001)
            wptBinHigh = h2.GetXaxis().FindFixBin(wptBins[iwpt+1] - 0.0001)
            intTot_inWptRange = h2.Integral(wptBinLow,wptBinHigh,0,1+dzBinsTot)
            intInDz_inWptRange = h2.Integral(wptBinLow,wptBinHigh,dzBinLow,dzBinHigh)
            #print "%.1f   %.1f" % (intInDz_inWptRange, intTot_inWptRange)
            eff = 0.0
            if intTot_inWptRange > 0.0:            
                eff = intInDz_inWptRange/intTot_inWptRange
            else:
                eff = 1.05
            hists[key].SetBinContent(iwpt+1, eff)            
            print("bin %d) Wpt in [%.0f, %.0f] GeV: eff = %.1f" % (iwpt+1, wptBins[iwpt], wptBins[iwpt+1],100.*eff))
        #
        #overflow
        intTot_inWptRange = h2.Integral(wptBinsTot+1,wptBinsTot+1,0,1+dzBinsTot)
        intInDz_inWptRange = h2.Integral(wptBinsTot+1,wptBinsTot+1,dzBinLow,dzBinHigh)    
        if intTot_inWptRange > 0.0:
            eff = 100. * intInDz_inWptRange/intTot_inWptRange
            print("Overflow %s > %.0f GeV: eff = %.1f" % (var,wptBins[-1],eff))
        else:
            print("Empty overflow bin")
        print("\n\n")
    
    #sortkeys = sortkeys[:5]

    adjustSettings_CMS_lumi()    
    createPlotDirAndCopyPhp(outdir)
    canvas = ROOT.TCanvas("canvas","",900,900)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.cd()
    if len(sortkeys) > len(colors):
        print("Warning: need to extend color array")
        quit()
    legendEntries = []
    for i,k in enumerate(sortkeys):
        legendEntries.append(k)
        hists[k].SetStats(0)
        hists[k].SetLineColor(colors[i])
        hists[k].SetMarkerColor(colors[i])
        hists[k].SetLineWidth(2)
        hists[k].SetMarkerStyle(markers[i])        
        print("%d) %s " % (i,k))
        if i:
            hists[k].Draw("LPSAME")
        else:
            hists[k].Draw("LP")
            hists[k].GetXaxis().SetTitle("W p_{T} [GeV]")
            if var == "mupt":
                hists[k].GetXaxis().SetTitle("preFSR muon p_{T} [GeV]")
            hists[k].GetXaxis().SetTitleOffset(1.2)
            hists[k].GetXaxis().SetTitleSize(0.05)
            hists[k].GetXaxis().SetLabelSize(0.04)
            hists[k].GetYaxis().SetTitle("efficiency: vertex dz(gen,reco) < %g mm" % (10.*absDzVal))
            hists[k].GetYaxis().SetTitleOffset(1.15)
            hists[k].GetYaxis().SetTitleSize(0.05)
            hists[k].GetYaxis().SetLabelSize(0.04)
            hists[k].GetYaxis().SetRangeUser(effLow,effHigh)
    leg = ROOT.TLegend(0.15,0.15,0.9,0.35+(0.05 if len(sortkeys) > 6 else 0.0))
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetFillColorAlpha(0,0.6)
    leg.SetShadowColor(0)
    leg.SetBorderSize(0)
    leg.SetNColumns(2)
    for il,le in enumerate(legendEntries):
        leg.AddEntry(hists[sortkeys[il]],"{n}) {s}".format(n=il,s=le),"LP")
    leg.Draw("same")

    canvas.RedrawAxis("sameaxis")
    setTDRStyle()
    #
    for ext in ["png","pdf"]:
        canvas.SaveAs(outdir + "/vertexEfficiency_{v}.{ext}".format(v=var,ext=ext))

    csel = ROOT.TCanvas("csel","",1000,900)
    csel.SetTickx(1)
    csel.SetTicky(1)
    csel.SetGrid(1)
    csel.cd()
    csel.SetLeftMargin(0.12)
    csel.SetRightMargin(0.04)
    csel.cd()
    heffSel_passDz.GetXaxis().SetTitle("Selection step")
    heffSel_passDz.GetYaxis().SetTitle("Selection efficiency")
    heffSel_passDz.SetStats(0)
    heffSel_passDz.SetLineColor(ROOT.kGreen+2)
    heffSel_passDz.SetFillColor(ROOT.kGreen+2)
    heffSel_passDz.SetFillStyle(3002)
    heffSel_passDz.SetMarkerColor(ROOT.kGreen+2)
    heffSel_passDz.SetLineWidth(2)
    heffSel_passDz.GetXaxis().SetTitleOffset(1.2)
    heffSel_passDz.GetXaxis().SetTitleSize(0.05)
    heffSel_passDz.GetXaxis().SetLabelSize(0.04)
    heffSel_passDz.GetYaxis().SetTitleOffset(1.15)
    heffSel_passDz.GetYaxis().SetTitleSize(0.05)
    heffSel_passDz.GetYaxis().SetLabelSize(0.04)
    heffSel_passDz.GetYaxis().SetRangeUser(0,1.2)
    heffSel_passDz.Draw("HIST")    
    heffSel_failDz.SetLineColor(ROOT.kRed+2)
    heffSel_failDz.SetMarkerColor(ROOT.kRed+2)
    heffSel_failDz.SetLineWidth(2)
    heffSel_failDz.Draw("HISTSAME")

    leg = ROOT.TLegend(0.25,0.84,0.9,0.9)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetFillColorAlpha(0,0.6)
    leg.SetShadowColor(0)
    leg.SetBorderSize(0)
    leg.SetNColumns(2)
    leg.AddEntry(heffSel_passDz,"dz < %g mm" % (10.*absDzVal),"LF")
    leg.AddEntry(heffSel_failDz,"dz > %g mm" % (10.*absDzVal),"LF")
    leg.Draw("same")

    csel.RedrawAxis("sameaxis")
    setTDRStyle()

    lat = ROOT.TLatex()
    #lat.SetNDC();
    lat.SetTextFont(42)        
    lat.SetTextSize(0.03)
    lat.SetTextAlign(12)
    lat.SetTextColor(ROOT.kBlack)
    for ib in range(1,1+heffSel_passDz.GetNbinsX()):
        if ib == 1:
            lat.DrawLatex(-0.2+heffSel_passDz.GetXaxis().GetBinCenter(ib),
                          0.75,
                          "#varepsilon(i)/#varepsilon(i-1)")            
        else:
            a  = heffSel_passDz.GetBinContent(ib)
            da = heffSel_passDz.GetBinError(ib)
            b  = heffSel_passDz.GetBinContent(ib-1)
            db = heffSel_passDz.GetBinError(ib-1)
            ratio = a/b 
            lat.DrawLatex(-0.2+heffSel_passDz.GetXaxis().GetBinCenter(ib),
                          0.75,
                          "{:.1f}%".format(100.*ratio))
            if args.printUncertainty:
                err = ratio * ROOT.TMath.Sqrt((da*da)/(a*a) + (db*db)/(b*b))
                lat.DrawLatex(-0.2+heffSel_passDz.GetXaxis().GetBinCenter(ib),
                              0.7,
                              "#pm{:.1f}%".format(100.*err))
    lat.SetTextColor(ROOT.kRed+2)
    for ib in range(1,1+heffSel_failDz.GetNbinsX()):
        if ib == 1:
            lat.DrawLatex(-0.2+heffSel_passDz.GetXaxis().GetBinCenter(ib),
                          0.6,
                          "#varepsilon(i)/#varepsilon(i-1)")            
        else:
            a  = heffSel_failDz.GetBinContent(ib)
            da = heffSel_failDz.GetBinError(ib)
            b  = heffSel_failDz.GetBinContent(ib-1)
            db = heffSel_failDz.GetBinError(ib-1)
            ratio = a/b 
            lat.DrawLatex(-0.2+heffSel_failDz.GetXaxis().GetBinCenter(ib),
                          0.6,
                          "{:.1f}%".format(100.*ratio))
            if args.printUncertainty:
                err = ratio * ROOT.TMath.Sqrt((da*da)/(a*a) + (db*db)/(b*b))
                lat.DrawLatex(-0.2+heffSel_failDz.GetXaxis().GetBinCenter(ib),
                             0.55,
                             "#pm{:.1f}%".format(100.*err))
    #
    for ext in ["png","pdf"]:
        csel.SaveAs(outdir + "/selectionEfficiency_{v}.{ext}".format(v=var,ext=ext))

        
        
