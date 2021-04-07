#!/usr/bin/env python3

# compare data histograms per era from these files
# plots/testNanoAOD/testMuonPrefire/Zmumu_plus_EventsPerRun_ReReco//plots_test.root
# plots/testNanoAOD/testMuonPrefire/Zmumu_plus_EventsPerRun_UltraLegacy//plots_test.root

import re
import os, os.path
import logging
import argparse
import shutil

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *
sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

sys.path.append(os.getcwd() + "/lumiStuff/")
from runPerEra import runsForEra as rfe

logging.basicConfig(level=logging.INFO)

def getFromFile(fname, mydict, plots=[], procs=[]):
    f = ROOT.TFile.Open(fname)
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    for pl in plots:
        for p in procs:
            name = f"{pl}_{p}"
            mydict[pl] = f.Get(name)
            if mydict[pl] == None:
                logging.info(" Cannot find {name} in file {fname}")
                quit()
            mydict[pl].SetDirectory(0)
    f.Close()
    return mydict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfileUL", type=str, nargs=1)
    parser.add_argument("rootfileRR", type=str, nargs=1)
    parser.add_argument("outdir",   type=str, nargs=1)
    parser.add_argument("-n", '--norm-lumi', dest='normLumi',  action='store_true', help='Normalize event yields by run luminosity')
    args = parser.parse_args()

    outdir = args.outdir[0] + "/"
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
        htmlpath = "./templates/index.php"
        shutil.copy(htmlpath, outdir)

        
    histoUL = {}
    histoRR = {}

    histoUL = getFromFile(args.rootfileUL[0], histoUL, plots=["runNumber"], procs=["data_UL"])
    histoRR = getFromFile(args.rootfileRR[0], histoRR, plots=["runNumber"], procs=["data_RR"])
    k = "runNumber"   # only one key for now, ugly but fine

    intLumi = None
    lumiPerRun = {}

    tmpIntLumi = 0.0
    if args.normLumi:
        
        intLumi = histoUL[k].Clone("intLumi")
        intLumi.Reset("ICESM")
        print(f"Integrated luminosity = {intLumi.Integral()}")
        with open("./lumiStuff/runlumi.csv") as f:
            for line in f:
                if line.startswith('#'):
                    continue
                tokens = line.strip().split(',')
                run,lumi = tokens[0],tokens[1]
                run = int(run)
                lumi = float(lumi)
                lumiPerRun[run] = lumi

        runs = list(lumiPerRun.keys())
        for ibin in range(1, 1 + intLumi.GetNbinsX()):
            binCenter = int(intLumi.GetBinCenter(ibin)) # it is a run on the axis
            #print(f"Run {binCenter}")
            if binCenter in runs:
                lumi = float(lumiPerRun[binCenter])
                histoUL[k].SetBinContent(ibin, histoUL[k].GetBinContent(ibin) / lumi)
                histoRR[k].SetBinContent(ibin, histoRR[k].GetBinContent(ibin) / lumi)
                #print(f"Divide by {lumi}")
            else:
                lumi = 0.0
            tmpIntLumi += lumi
            intLumi.SetBinContent(ibin, tmpIntLumi)
            #print(f"Lumi = {lumi} --> integrated luminosity = {tmpIntLumi}")        
        intLumi.SetStats(0)

    averageRunB = 0.0
    countB = 0
    for ibin in range(1, 1 + histoUL[k].GetNbinsX()):
        binCenter = int(histoUL[k].GetBinCenter(ibin))
        if binCenter >= rfe["B"][0] and binCenter <= rfe["B"][-1]:
            if histoUL[k].GetBinContent(ibin) > 0:
                averageRunB += histoUL[k].GetBinContent(ibin)
                countB += 1
        elif binCenter > rfe["B"][-1]:
            break
    averageRunB = averageRunB / countB
        
    adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas","",3000,1000)
    # for k in list(histoUL.keys()):        
    #     hists = [histoUL[k], histoRR[k]]
    #     legEntries = ["Data UltraLegacy", "Data ReReco"]
    #     # note: for 2 histograms h[0] is the numerator, while for more it becomes denominator
    #     drawNTH1(hists, legEntries,
    #              hists[0].GetXaxis().GetTitle(), "Events (normalized to int. lumi)",
    #              k,
    #              outdir,
    #              draw_both0_noLog1_onlyLog2=1,
    #              leftMargin=0.06,
    #              labelRatioTmp="UL/RR::0.95,1.05",
    #              legendCoords="0.12,0.52,0.84,0.9;2",
    #              passCanvas=canvas,
    #              drawLumiLatex=True,
    #              drawVertLines=runRanges
    #     )
    hists = [histoUL[k], histoRR[k]]                                                                             
    legEntries = ["Data UltraLegacy", "Data ReReco"] 

    labelXtmp = hists[0].GetXaxis().GetTitle()
    labelYtmp = "Events / fb^{-1}" if args.normLumi else "Events"
    draw_both0_noLog1_onlyLog2=1
    leftMargin = 0.06
    rightMargin = 0.04
    labelRatioTmp = "UL/RR::0.95,1.05"
    legendCoords = "0.06,0.48,0.84,0.9;2"
    yAxisTitleOffset = 0.6
    lowerPanelHeight = 0.3
    canvasName = f"{k}_normLumi" if args.normLumi else k
    
    # run ranges for eras from B to H, prior to actual json cut (array has beginning of each except for B)
    textForLines = sorted(list(rfe.keys())) # this is correct sorting
    textForLines = [x for x in textForLines if x != "F_postVFP"] # too short of a span, remove it
    #print(textForLines)
    runStart = [rfe[k][0] for k in textForLines]
    labelPosition = [int(run+200) for run in runStart]
    
    moreText = ""
    moreTextLatex = ""
    lumi = 35.9
    skipLumi = False
    drawLumiLatex = True
    drawLineLowerPanel = ""
    draw_both0_noLog1_onlyLog2 = 1
    
    xAxisName,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    yAxisName,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    yRatioAxisName,setRatioYAxisRangeFromUser,yminRatio,ymaxRatio = getAxisRangeFromUser(labelRatioTmp)
    
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.cd()

    pad2 = 0
    if lowerPanelHeight: 
        canvas.SetBottomMargin(lowerPanelHeight)
        pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
        pad2.SetTopMargin(1-lowerPanelHeight)
        pad2.SetRightMargin(rightMargin)
        pad2.SetLeftMargin(leftMargin)
        pad2.SetFillColor(0)
        pad2.SetGridy(1)
        pad2.SetFillStyle(0)
    else:
        canvas.SetBottomMargin(0.15)


    h1 = hists[0]
    hnums = [hists[i] for i in range(1,len(hists))]
    frame = h1.Clone("frame")
    frame.GetXaxis().SetLabelSize(0.04)
    frame.SetStats(0)

    h1.SetLineColor(ROOT.kBlack)
    h1.SetMarkerColor(ROOT.kBlack)
    h1.SetMarkerStyle(20)
    #h1.SetMarkerSize(0)

    colors = [ROOT.kRed+2, ROOT.kBlue, ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kAzure+2, ROOT.kPink+7]
    for ic,h in enumerate(hnums):
        # h.SetLineColor(colors[ic])
        # h.SetFillColor(colors[ic])
        # if ic==0: h.SetFillStyle(3004)   
        # if ic==2: h.SetFillStyle(3002)   
        # h.SetFillColor(colors[ic])
        # h.SetMarkerSize(0)
        h.SetLineColor(colors[ic])
        h.SetFillColor(colors[ic])
        h.SetMarkerSize(0)
        if ic==0: 
            h.SetFillStyle(3004)   
        if ic==1: 
            h.SetFillColor(0) 
            h.SetLineWidth(2) 
        if ic==2: 
            h.SetFillStyle(3002)           
        if ic==3:
            h.SetFillColor(0)
            h1.SetMarkerColor(ROOT.kGray+3)
            h1.SetMarkerStyle(25)
            #h1.SetMarkerSize(2)
            
    if ymin == ymax == 0.0:
        ymin = 9999.9
        ymax = -9999.9
        for h in hists:
            if h.GetBinContent(h.GetMaximumBin()) > ymax: ymax = h.GetBinContent(h.GetMaximumBin())
            if h.GetBinContent(h.GetMinimumBin()) < ymin: ymin = h.GetBinContent(h.GetMinimumBin())
        if ymin < 0: ymin = 0
        ymax *= (1.5 if args.normLumi else 1.2)
        
    if lowerPanelHeight:
        h1.GetXaxis().SetLabelSize(0)
        h1.GetXaxis().SetTitle("")  
    else:
        h1.GetXaxis().SetTitle(xAxisName)
        h1.GetXaxis().SetTitleOffset(1.2)
        h1.GetXaxis().SetTitleSize(0.05)
        h1.GetXaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetTitle(yAxisName)
    h1.GetYaxis().SetTitleOffset(yAxisTitleOffset) 
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetRangeUser(ymin, ymax)
    h1.GetYaxis().SetTickSize(0.01)
    if setXAxisRangeFromUser: h1.GetXaxis().SetRangeUser(xmin,xmax)
    h1.Draw("PE")
    # draw them below
    #for h in hnums:
    #    h.Draw("HIST SAME")

    nColumnsLeg = 1
    if ";" in legendCoords: 
        nColumnsLeg = int(legendCoords.split(";")[1])
    legcoords = [float(x) for x in (legendCoords.split(";")[0]).split(',')]
    lx1,lx2,ly1,ly2 = legcoords[0],legcoords[1],legcoords[2],legcoords[3]
    leg = ROOT.TLegend(lx1,ly1,lx2,ly2)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetFillColorAlpha(0,0.6)
    leg.SetShadowColor(0)
    leg.SetBorderSize(0)
    leg.SetNColumns(nColumnsLeg)
    for il,le in enumerate(legEntries):
        leg.AddEntry(hists[il],le,"PE" if il == 0 else "FL")
    leg.Draw("same")
    canvas.RedrawAxis("sameaxis")

    for htmp in hists:
        htmp.SetStats(0)

    # overlay int lumi plot
    rightmax = intLumi.GetBinContent(intLumi.GetMaximumBin())
    print(f"intLumi.Integral() = {intLumi.Integral()}")
    print(f"h1.GetBinContent(h1.GetMaximumBin())/rightmax = {h1.GetBinContent(h1.GetMaximumBin())} / {rightmax}")
    scale = h1.GetBinContent(h1.GetMaximumBin())/rightmax
    intLumi.SetLineColor(ROOT.kAzure)
    intLumi.SetFillColor(ROOT.kCyan)
    intLumi.SetFillColorAlpha(ROOT.kCyan, 0.5)
    intLumi.Scale(scale)
    intLumi.Draw("LF same")
    # draw an axis on the right side (doesn't seem to work, though)
    axis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
                       ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(),
                       0, rightmax, 510, "+L")
    axis.SetLineColor(ROOT.kBlue)
    axis.SetTextColor(ROOT.kBlue)
    axis.SetTitle("Integrated Luminosity (fb^{-1 })")
    axis.SetTitleOffset(0.1)
    axis.SetLabelSize(0.04)
    axis.Draw("same");
    legIntLumi = ROOT.TLegend(0.53,ly1,0.9,ly2)
    legIntLumi.SetFillColor(0)
    legIntLumi.SetFillStyle(0)
    legIntLumi.SetFillColorAlpha(0, 0.6)
    legIntLumi.SetShadowColor(0)
    legIntLumi.SetBorderSize(0)
    legIntLumi.SetNColumns(2)
    legIntLumi.AddEntry(intLumi, "Integrated luminosity (fb^{-1 })", "LF")
    h1.Draw("PE SAME")
    for h in hnums:
        h.Draw("HIST SAME")

    lineAverageB = ROOT.TF1("horiz_line",str(averageRunB),h1.GetXaxis().GetBinLowEdge(1),h1.GetXaxis().GetBinLowEdge(h1.GetNbinsX()+1))
    lineAverageB.SetLineColor(ROOT.kGray+3)
    lineAverageB.SetLineWidth(1)
    lineAverageB.SetLineStyle(2)
    lineAverageB.Draw("Lsame")
    legIntLumi.AddEntry(lineAverageB, "Average Run B", "L")
    legIntLumi.Draw("same")
        
    canvas.RedrawAxis("sameaxis")

    vertline = ROOT.TLine(36,0,36,canvas.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(3)
    bintext = ROOT.TLatex()
    #bintext.SetNDC()
    bintext.SetTextSize(0.05)  # 0.03
    bintext.SetTextFont(42)

    if len(runStart):
        nLines = len(runStart)
        for i in range(nLines):
            vertline.DrawLine(runStart[i],0,runStart[i],ymax)
        if len(textForLines):
            for i in range(len(textForLines)): # we need nLines
                #texoffset = 0.1 * (4 - (i%4))
                #ytext = (1. + texoffset)*ymax/2.  
                ytext = 0.75 * ymax
                bintext.DrawLatex(labelPosition[i], ytext, textForLines[i])

    # redraw legend, or vertical lines appear on top of it
    leg.Draw("same")

    if len(moreText):
        realtext = moreText.split("::")[0]
        x1,y1,x2,y2 = 0.7,0.8,0.9,0.9
        if "::" in moreText:
            x1,y1,x2,y2 = (float(x) for x in (moreText.split("::")[1]).split(","))
        pavetext = ROOT.TPaveText(x1,y1,x2,y2,"NB NDC")
        for tx in realtext.split(";"):
            pavetext.AddText(tx)
        pavetext.SetFillColor(0)
        pavetext.SetFillStyle(0)
        pavetext.SetBorderSize(0)
        pavetext.SetLineColor(0)
        pavetext.Draw("same")

    if len(moreTextLatex):
        realtext = moreTextLatex.split("::")[0]
        x1,y1,ypass,textsize = 0.75,0.8,0.08,0.035
        if "::" in moreTextLatex:
            x1,y1,ypass,textsize = (float(x) for x in (moreTextLatex.split("::")[1]).split(","))            
        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)        
        lat.SetTextSize(textsize)
        for itx,tx in enumerate(realtext.split(";")):
            lat.DrawLatex(x1,y1-itx*ypass,tx)

    setTDRStyle()
    if not skipLumi:
        if not drawLumiLatex:
            if lumi != None: CMS_lumi(canvas,lumi,True,False)
            else:            CMS_lumi(canvas,"",True,False)
        else:
            latCMS = ROOT.TLatex()
            latCMS.SetNDC();
            latCMS.SetTextFont(42)
            latCMS.SetTextSize(0.045)
            latCMS.DrawLatex(0.1, 0.95, '#bf{CMS} #it{Preliminary}')
            if lumi != None: latCMS.DrawLatex(0.85, 0.95, '%s fb^{-1} (13 TeV)' % lumi)
            else:            latCMS.DrawLatex(0.90, 0.95, '(13 TeV)')

    if lowerPanelHeight:
        pad2.Draw()
        pad2.cd()

        frame.Reset("ICES")
        if setRatioYAxisRangeFromUser: frame.GetYaxis().SetRangeUser(yminRatio,ymaxRatio)
        #else:                          
        #frame.GetYaxis().SetRangeUser(0.5,1.5)
        frame.GetYaxis().SetNdivisions(5)
        frame.GetYaxis().SetTitle(yRatioAxisName)
        frame.GetYaxis().SetTitleOffset(yAxisTitleOffset)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().CenterTitle()
        frame.GetXaxis().SetTitle(xAxisName)
        if setXAxisRangeFromUser: frame.GetXaxis().SetRangeUser(xmin,xmax)
        frame.GetXaxis().SetTitleOffset(1.2)
        frame.GetXaxis().SetTitleSize(0.05)

        if len(hists) == 2:
            ratio = h1.Clone("ratio")
            den = hnums[0].Clone("den")
            den_noerr = hnums[0].Clone("den_noerr")
            for iBin in range (1,den_noerr.GetNbinsX()+1):
                den_noerr.SetBinError(iBin,0.)
            den.Divide(den_noerr)
            ratio.Divide(den_noerr)
            den.SetFillColor(ROOT.kGray)
            den.SetFillStyle(1001)
            #den_noerr.SetFillColor(ROOT.kGray)
            frame.Draw()
            frame.SetMarkerSize(0)
            frame.SetMarkerStyle(0) # important to remove dots at y = 1
            den.Draw("E2same")
            ratio.Draw("EPSAME")
        else:
            ratio = h1.Clone("ratio")
            den_noerr = h1.Clone("den_noerr")
            for iBin in range (1,den_noerr.GetNbinsX()+1):
                den_noerr.SetBinError(iBin,0.)

            ratio.Divide(den_noerr)
            ratio.SetFillColor(ROOT.kGray)
            ratio.SetFillStyle(1001)
            #den_noerr.SetFillColor(ROOT.kGray)
            frame.Draw()
            ratio.SetMarkerSize(0)
            ratio.SetMarkerStyle(0) # important to remove dots at y = 1
            ratio.Draw("E2same")

            ratios = []
            for i,h in enumerate(hnums):
                ratios.append(h.Clone("ratio_"+str(i+1)))
                ratios[-1].Divide(den_noerr)
                #ratios[-1].SetLineColor(h.GetLineColor())
                #ratios[-1].SetMarkerSize(0)
                #ratios[-1].SetMarkerStyle(0)
                #ratios[-1].SetFillColor(0)
                if h.GetFillColor():
                    ratios[-1].Draw("E2 SAME")
                else:
                    ratios[-1].Draw("HIST SAME")

        line = ROOT.TF1("horiz_line","1",ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
        line.SetLineColor(ROOT.kBlack)
        line.SetLineWidth(1)
        line.Draw("Lsame")
        
        if drawLineLowerPanel:
            legEntry,yline = drawLineLowerPanel.split('::')
            line2 = ROOT.TF1("horiz_line_2",str(1+float(yline)),ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
            line3 = ROOT.TF1("horiz_line_3",str(1-float(yline)),ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
            line2.SetLineColor(ROOT.kBlue)
            line2.SetLineWidth(1)
            line2.Draw("Lsame")
            line3.SetLineColor(ROOT.kBlue)
            line3.SetLineWidth(1)
            line3.Draw("Lsame")
            x1leg2 = 0.2 if leftMargin > 0.1 else 0.07
            x2leg2 = 0.5 if leftMargin > 0.1 else 0.27
            y1leg2 = 0.25 if leftMargin > 0.1 else 0.3
            y2leg2 = 0.35 if leftMargin > 0.1 else 0.35
            leg2 = ROOT.TLegend(x1leg2, y1leg2, x2leg2, y2leg2)
            leg2.SetFillColor(0)
            leg2.SetFillStyle(0)
            leg2.SetBorderSize(0)
            leg2.AddEntry(line2,legEntry,"L")
            leg2.Draw("same")

        
        pad2.RedrawAxis("sameaxis")


    if draw_both0_noLog1_onlyLog2 != 2:
        canvas.SaveAs(outdir + canvasName + ".png")
        canvas.SaveAs(outdir + canvasName + ".pdf")

    if draw_both0_noLog1_onlyLog2 != 1:        
        if yAxisName == "a.u.": 
            h1.GetYaxis().SetRangeUser(max(0.0001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        else:
            h1.GetYaxis().SetRangeUser(max(0.001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        canvas.SetLogy()
        canvas.SaveAs(outdir + canvasName + "_logY.png")
        canvas.SaveAs(outdir + canvasName + "_logY.pdf")
        canvas.SetLogy(0)
