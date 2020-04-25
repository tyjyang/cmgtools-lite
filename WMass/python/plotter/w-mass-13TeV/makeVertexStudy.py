#!/bin/env python

import os, re, array, math
import time

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

# plots/vertexStudy/ask1orMoreLep
outdir = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/rel_slc7/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/plots/vertexStudy/newNtuplesNoSkim_v3/compareEfficiency"

inputfolder = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/rel_slc7/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/plots/vertexStudy/newNtuplesNoSkim_v3/"
# inputfiles = {"1::genEtaPt"                : "W_94X_noRecoCuts_genEtaPt26to56",
#               "2::fullSel_NoMtNoIsoNoID"   : "W_94X_allRecoCuts_noMtNoIsoNoID",
#               "3::fullSel_NoMtNoIso"       : "W_94X_allRecoCuts_noMtNoIso",
#               "4::fullSel_NoMt"            : "W_94X_allRecoCuts_noMt",
#               "5::fullSel"                 : "W_94X_allRecoCuts",
#               "6::fullSel_NoTrigger"       : "W_94X_allRecoCuts_noTrigger",
#               "7::1lepInAccept"            : "W_94X_1lepInAccept",
#           }

workingPoints = ["alwaystrue", "genEtaPt", "vertexPresel", "muonInAccept", "muMediumId", "muTightIso", "mtl1pf40", "trigger"]

absDzVal = 0.1 # 1 mm        
var = "Wpt"
#var = "dressedLepPt"
genLepPtLow = 26.0
genLepPtHigh = 56.0
hname = "dzVertex_gen_primary__{v}_Wnopt".format(v=var)
# Wpt from 0 to 100 with 2 GeV width, dz from -1.0 to 1.0 with 0.1 width
# dressed lepton pT from 26 to 56 with 2 GeV width, dz from -1.0 to 1.0 with 0.1 width

if __name__ == "__main__":

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

    for ik,key in enumerate(sortkeys):

        print "="*30
        print "Cut: " + str(key)
        print "-"*30
        print
        inputfile = inputfolder + key + "/plots_forTest.root"
        tf = ROOT.TFile.Open(inputfile)        
        h2 =   tf.Get(hname)
        h2.SetDirectory(0)
        tf.Close()

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
            heffSel_passDz.SetBinContent(ik,intInDz/backupIntInDz_noCuts)
            heffSel_failDz.SetBinContent(ik,(intTot-intInDz)/backupIntOutDz_noCuts)

        print "Inclusive efficiency: %.1f" % (100.*intInDz/intTot)

        wptBins = [4.0 * i for i in range(0,26)]
        if var == "dressedLepPt":
            wptBins = [genLepPtLow + 1.0 * i for i in range(0,1+int(genLepPtHigh-genLepPtLow+0.001))]
        
        hists[key] = ROOT.TH1D("eff_"+key,"",len(wptBins)-1,array("d",wptBins))
        #gr[key] = 

        for iwpt in range(len(wptBins)-1):

            wptBinLow = h2.GetXaxis().FindFixBin(wptBins[iwpt] + 0.0001)
            wptBinHigh = h2.GetXaxis().FindFixBin(wptBins[iwpt+1] - 0.0001)
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
        if intTot_inWptRange > 0.0:
            eff = 100. * intInDz_inWptRange/intTot_inWptRange
            print "Overflow %s > %.0f GeV: eff = %.1f" % (var,wptBins[-1],eff)
        else:
            print "Empty overflow bin"
        print "\n"*2
    
    adjustSettings_CMS_lumi()    
    createPlotDirAndCopyPhp(outdir)
    canvas = ROOT.TCanvas("canvas","",900,900)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.cd()
    colors = [ROOT.kBlack, ROOT.kCyan+1, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kMagenta, ROOT.kAzure+2]
    markers = [ROOT.kFullCircle, ROOT.kFullCross, ROOT.kOpenSquare, ROOT.kOpenTriangleUp, ROOT.kOpenCircle, ROOT.kOpenSquareDiagonal, ROOT.kFullCircle, ROOT.kFullSquare]
    if len(sortkeys) > len(colors):
        print "Warning: need to extend color array"
        quit()
    legendEntries = []
    for i,k in enumerate(sortkeys):
        legendEntries.append(k)
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
            if var == "dressedLepPt":
                hists[k].GetXaxis().SetTitle("dressed lepton p_{T} [GeV]")
            hists[k].GetXaxis().SetTitleOffset(1.2)
            hists[k].GetXaxis().SetTitleSize(0.05)
            hists[k].GetXaxis().SetLabelSize(0.04)
            hists[k].GetYaxis().SetTitle("efficiency: vertex dz(gen,reco) < %g mm" % (10.*absDzVal))
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
                          0.7,
                          "#varepsilon(i)/#varepsilon(i-1)")            
        else:
            lat.DrawLatex(-0.2+heffSel_passDz.GetXaxis().GetBinCenter(ib),
                          0.7,
                          "{:.1f}%".format(100.*heffSel_passDz.GetBinContent(ib)/heffSel_passDz.GetBinContent(ib-1)))
    lat.SetTextColor(ROOT.kRed+2)
    for ib in range(1,1+heffSel_failDz.GetNbinsX()):
        if ib == 1:
            lat.DrawLatex(-0.2+heffSel_passDz.GetXaxis().GetBinCenter(ib),
                          0.6,
                          "#varepsilon(i)/#varepsilon(i-1)")            
        else:
            lat.DrawLatex(-0.2+heffSel_failDz.GetXaxis().GetBinCenter(ib),
                          0.6,
                          "{:.1f}%".format(100.*heffSel_failDz.GetBinContent(ib)/heffSel_failDz.GetBinContent(ib-1)))
    #
    for ext in ["png","pdf"]:
        csel.SaveAs(outdir + "/selectionEfficiency.{ext}".format(v=var,ext=ext))

        
