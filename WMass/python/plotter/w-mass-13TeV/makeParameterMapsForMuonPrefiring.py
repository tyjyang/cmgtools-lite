#!/usr/bin/env python3

# transform TF1 parametrization of muon prefiring probability provided by Jan into histograms containing those parameters

import os, re
import argparse
from array import array

import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

import utilities
utilities = utilities.util()

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(1)
ROOT.TH1.SetDefaultSumw2()

muonprefiringfile = "testMuonSF/L1MuonPrefiringParametriations.root"
output_muonprefiringfile = "testMuonSF/L1MuonPrefiringParametriations_histograms.root"
etabins = ["0.0", "0.2", "0.3", "0.55", "0.83", "1.24", "1.4", "1.6", "1.8", "2.1", "2.25", "2.4"]
etabinsfloat = [round(float(x), 2) for x in etabins]
plotFolder = "plots/testNanoAOD/testMuonPrefire/newParametricPrefiring_etaPt/"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    args = parser.parse_args()

    createPlotDirAndCopyPhp(plotFolder)
    
    #print(etabinsfloat)
    infile = safeOpenFile(muonprefiringfile, mode="READ")
    eras = ["preVFP", "postVFP", "BG", "H"]
    # check that the pseudodata histogram does not exist, unless one wants to update it
    htemplate = ROOT.TH2D("htemplate", "", len(etabinsfloat)-1, array('d', etabinsfloat), 3, 0.5, 3.5)
    outfile = safeOpenFile(output_muonprefiringfile, mode="RECREATE")

    hplot = ROOT.TH2D("hplot", "", len(etabinsfloat)-1, array('d', etabinsfloat), 100, 10, 60)
    adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas", "", 800, 800) 

    fhotspot = safeGetObject(infile, "L1prefiring_muonparam_HotSpot_2016", quitOnFail=True, silent=True, detach=False)
    hhotspot = ROOT.TH2D("L1prefiring_muonparam_2016_hotspot", "1.24 < #eta < 1.6, 2.44346 < #phi < 2.79253", 1, 1.24, 1.6, 3, 0.5, 3.5)
    hplot.Reset("ICESM")
    hplot.SetTitle(f"2016 hotspot, 1.24 < #eta < 1.6, 2.44346 < #phi < 2.79253")
    hhotspot.SetBinContent(1, 1, fhotspot.GetParameter(0))
    hhotspot.SetBinContent(1, 2, fhotspot.GetParameter(1))
    hhotspot.SetBinContent(1, 3, fhotspot.GetParameter(2))
    hhotspot.SetBinError(1, 1, fhotspot.GetParError(0))
    hhotspot.SetBinError(1, 2, fhotspot.GetParError(1))
    hhotspot.SetBinError(1, 3, fhotspot.GetParError(2))
    ibin = hplot.GetXaxis().FindFixBin(1.3)
    for ipt in range(1, hplot.GetNbinsY()+1):
        val = fhotspot.Eval(hplot.GetYaxis().GetBinCenter(ipt))
        hplot.SetBinContent(ibin,   ipt, val)
        hplot.SetBinContent(ibin+1, ipt, val)
    hhotspot.Write()
    drawCorrelationPlot(hplot,
                        "Muon #eta",
                        "Muon p_{T} (GeV)",
                        f"Prefiring probability::0.0,{hplot.GetBinContent(hplot.GetMaximumBin())}",
                        f"muonPrefiringProbability_2016_hotspot", plotLabel="ForceTitle", outdir=plotFolder,
                        draw_both0_noLog1_onlyLog2=1, nContours=51, palette=87,
                        passCanvas=canvas, skipLumi=True)
    
    
    hPrefVsEta_pt40 = {}
    
    for era in eras:
        hplot.SetTitle(f"2016{era}")
        h = copy.deepcopy(htemplate.Clone(f"L1prefiring_muonparam_2016{era}"))
        hplot.Reset("ICESM")
        for ib,binLow in enumerate(etabins[:-1]):
            ibin = ib+1
            name = f"L1prefiring_muonparam_{binLow}To{etabins[ibin]}_2016{era}"
            # print(f"{binLow}   {name}")
            ftmp = safeGetObject(infile, name, quitOnFail=True, silent=True, detach=False)
            h.SetBinContent(ibin, 1, ftmp.GetParameter(0))
            h.SetBinContent(ibin, 2, ftmp.GetParameter(1))
            h.SetBinContent(ibin, 3, ftmp.GetParameter(2))
            h.SetBinError(ibin, 1, ftmp.GetParError(0))
            h.SetBinError(ibin, 2, ftmp.GetParError(1))
            h.SetBinError(ibin, 3, ftmp.GetParError(2))
            for ipt in range(1, hplot.GetNbinsY()+1):
                hplot.SetBinContent(ibin, ipt, ftmp.Eval(hplot.GetYaxis().GetBinCenter(ipt)))
                hplot.SetBinError(ibin, ipt, ftmp.GetParError(2)) # just use uncertainty on the plateau
        h.Write()

        ptBin = hplot.GetYaxis().FindFixBin(40.1)
        hPrefVsEta_pt40[era] = hplot.ProjectionX(f"muonPrefiringVsEta_pt40GeV_{era}", ptBin, ptBin, "e")
        
        drawCorrelationPlot(hplot,
                            "Muon |#eta|",
                            "Muon p_{T} (GeV)",
                            "Prefiring probability::0.0,0.035",
                            f"muonPrefiringProbability_2016{era}", plotLabel="ForceTitle", outdir=plotFolder,
                            draw_both0_noLog1_onlyLog2=1, nContours=51, palette=87,
                            passCanvas=canvas, skipLumi=True)

    hlist = [hPrefVsEta_pt40[era] for era in ["BG", "postVFP", "H"]]
    adjustSettings_CMS_lumi()
    drawNTH1(hlist, ["preVFP", "postVFP", "H"],
             "Muon |#eta|", "Muon L1 prefiring probability",
             "muonPrefiringProbability_prePostVFP",
             plotFolder,
             draw_both0_noLog1_onlyLog2=1,
             #legendCoords="0.18,0.9,0.82,0.9;3",
             legendCoords="0.18,0.48,0.63,0.9;1;Data 2016",
             lowerPanelHeight=0,
             passCanvas=canvas,
             lumi=36.3,
             onlyLineColor=True,
             drawErrorAll=True,
             #markerStyleFirstHistogram=25,
             useLineFirstHistogram=True,
             #moreTextLatex="25 < p_{T}^{#mu} < 60 GeV::0.18,0.65,0.08,0.05"  # if legend has 3 columns
             moreTextLatex="25 < p_{T}^{#mu} < 60 GeV::0.55,0.85,0.08,0.05"
    )


    print(f"Writing histograms in {output_muonprefiringfile}")
    outfile.Close()
    infile.Close()
