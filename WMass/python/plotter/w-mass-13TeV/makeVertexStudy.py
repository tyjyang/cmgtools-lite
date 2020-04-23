#!/bin/env python

import ROOT, os, sys, re, array, math
import time

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gROOT.SetBatch(True)

# plots/vertexStudy/ask1orMoreLep
outdir = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/rel_slc7/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/plots/vertexStudy/ask1orMoreLep/compareEfficiency"

inputfolder = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/rel_slc7/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/plots/vertexStudy/ask1orMoreLep/"
inputfiles = {"1::genEtaPt"                : "W_94X_noRecoCuts_genEtaPt26to56",
              "2::fullSel_NoMtNoIsoNoID"   : "W_94X_allRecoCuts_noMtNoIsoNoID",
              "3::fullSel_NoMtNoIso"       : "W_94X_allRecoCuts_noMtNoIso",
              "4::fullSel_NoMt"            : "W_94X_allRecoCuts_noMt",
              "5::fullSel"                 : "W_94X_allRecoCuts",
              "6::fullSel_NoTrigger"       : "W_94X_allRecoCuts_noTrigger",
              "7::1lepInAccept"            : "W_94X_1lepInAccept",
          }


absDzVal = 0.1 # 1 mm        
hname = "dzVertex_gen_primary__Wpt_Wnopt" # Wpt from 0 to 100 with 2 GeV width, dz from -1.0 to 1.0 with 0.1 width

if __name__ == "__main__":

    sortkeys = inputfiles.keys()
    sortkeys = sorted(sortkeys, key = lambda x: int(x.split("::")[0]))
    #print sortkeys
    #quit()

    hists = {}

    for key in sortkeys:

        inputfile = inputfolder + inputfiles[key] + "/plots_forTest.root"
        tf = ROOT.TFile.Open(inputfile)        
        h2 =   tf.Get(hname)
        h2.SetDirectory(0)
        tf.Close()

        # 0.001 is an offset to get bin to the right of left of the edge
        dzBinLow = h2.GetYaxis().FindFixBin(-1.0 * absDzVal + 0.001) # get bin with -absDzVal as left edge
        dzBinHigh = h2.GetYaxis().FindFixBin(absDzVal - 0.001) # get bin with absDzVal as right edge
        dzBinsTot = h2.GetNbinsY()
        wptBinsTot = h2.GetNbinsX()
        intTot =  h2.Integral(0,1+wptBinsTot,0,1+dzBinsTot)# include under/over-flow bins
        intInDz = h2.Integral(0,1+wptBinsTot,dzBinLow,dzBinHigh)# include under/over-flow bins
        wptlow = h2.GetXaxis().GetBinLowEdge(1)
        wpthigh = h2.GetXaxis().GetBinLowEdge(wptBinsTot+1)

        print "Inclusive efficiency: %.1f" % (100.*intInDz/intTot)

        wptBins = [4.0 * i for i in range(0,26)]
        
        hists[key] = ROOT.TH1D("eff_"+key.split("::")[0],"",len(wptBins)-1,array("d",wptBins))
        #gr[key] = 

        for iwpt in range(len(wptBins)-1):

            wptBinLow = h2.GetXaxis().FindFixBin(wptBins[iwpt] + 0.01)
            wptBinHigh = h2.GetXaxis().FindFixBin(wptBins[iwpt+1] - 0.01)
            intTot_inWptRange = h2.Integral(wptBinLow,wptBinHigh,0,1+dzBinsTot)
            intInDz_inWptRange = h2.Integral(wptBinLow,wptBinHigh,dzBinLow,dzBinHigh)
            #print "%.1f   %.1f" % (intInDz_inWptRange, intTot_inWptRange)
            eff = intInDz_inWptRange/intTot_inWptRange
            hists[key].SetBinContent(iwpt+1, eff)            
            print "bin %d) Wpt in [%.0f, %.0f] GeV: eff = %.1f" % (iwpt+1, wptBins[iwpt], wptBins[iwpt+1],100.*eff)
        #
        #overflow
        intTot_inWptRange = h2.Integral(wptBinsTot+1,wptBinsTot+1,0,1+dzBinsTot)
        intInDz_inWptRange = h2.Integral(wptBinsTot+1,wptBinsTot+1,dzBinLow,dzBinHigh)    
        eff = 100. * intInDz_inWptRange/intTot_inWptRange
        print "Overflow Wpt > %.0f GeV: eff = %.1f" % (wptBins[-1],eff)
    
    adjustSettings_CMS_lumi()    
    createPlotDirAndCopyPhp(outdir)
    canvas = ROOT.TCanvas("canvas","",900,900)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.cd()
    colors = [ROOT.kBlack, ROOT.kCyan, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kMagenta]
    markers = [ROOT.kFullCircle, ROOT.kOpenCross, ROOT.kFullSquare, ROOT.kOpenTriangleUp, ROOT.kOpenCircle, ROOT.kOpenSquareDiagonal, ROOT.kFullCircle]
    legendEntries = []
    for i,k in enumerate(sortkeys):
        legendEntries.append(k.split("::")[1])
        hists[k].SetStats(0)
        hists[k].SetLineColor(colors[i])
        hists[k].SetMarkerColor(colors[i])
        hists[k].SetLineWidth(2)
        hists[k].SetMarkerStyle(markers[i])
        if i:
            hists[k].Draw("LPSAME")
        else:
            hists[k].Draw("LP")
            hists[k].GetXaxis().SetTitle("W p_{T} [GeV]")
            hists[k].GetXaxis().SetTitleOffset(1.2)
            hists[k].GetXaxis().SetTitleSize(0.05)
            hists[k].GetXaxis().SetLabelSize(0.04)
            hists[k].GetYaxis().SetTitle("efficiency: vertex dz(gen,reco) < 1 mm")
            hists[k].GetYaxis().SetTitleOffset(1.15)
            hists[k].GetYaxis().SetTitleSize(0.05)
            hists[k].GetYaxis().SetLabelSize(0.04)
            hists[k].GetYaxis().SetRangeUser(0.7,1)
    leg = ROOT.TLegend(0.15,0.15,0.9,0.35)
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
        canvas.SaveAs(outdir + "/vertexEfficiency.{ext}".format(ext=ext))
