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
    canvas = ROOT.TCanvas("canvas", "", 800, 800) 
    
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
        h.Write()
        drawCorrelationPlot(hplot,
                            "Muon |#eta|",
                            "Muon p_{T} (GeV)",
                            "Prefiring probability::0.0,0.035",
                            f"muonPrefiringProbability_2016{era}", plotLabel="ForceTitle", outdir=plotFolder,
                            draw_both0_noLog1_onlyLog2=1, nContours=51, palette=87,
                            passCanvas=canvas, skipLumi=True)


    print(f"Writing histograms in {output_muonprefiringfile}")
    outfile.Close()
    infile.Close()
