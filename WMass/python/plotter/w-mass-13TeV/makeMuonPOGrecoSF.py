#!/usr/bin/env python3

import json
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

muonfiles = {"preVFP"  : "testMuonSF/Efficiency_muon_generalTracks_Run2016preVFP_UL_trackerMuon.json",
             "postVFP" : "testMuonSF/Efficiency_muon_generalTracks_Run2016postVFP_UL_trackerMuon.json"
}


plotFolder = "plots/testNanoAOD/testMuonPOGrecoSF/sf/"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    args = parser.parse_args()

    createPlotDirAndCopyPhp(plotFolder)

    histo = {}
    canvas = ROOT.TCanvas("canvas","",800,800)
    
    # Opening JSON file
    for k in muonfiles.keys():

        with open(muonfiles[k]) as json_file:
            foo = json.load(json_file)
            #print(foo["NUM_TrackerMuons_DEN_genTracks"]["abseta_pt"])
            data = foo["NUM_TrackerMuons_DEN_genTracks"]["abseta_pt"]

            etabinsfloat = []
            ptbinsfloat = []

            for eta in data.keys():
                #print(eta)
                etaRange = eta.split(":")[1].replace("[","").replace("]","")
                #print(etaRange)
                for x in etaRange.split(","):
                    if round(float(x),2) not in etabinsfloat:
                        etabinsfloat.append(round(float(x),2))  
                etabinsfloat = sorted(etabinsfloat)
                for pt in data[eta].keys():
                    ptRange = pt.split(":")[1].replace("[","").replace("]","")
                    for x in ptRange.split(","):
                        if round(float(x),2) not in ptbinsfloat:
                            ptbinsfloat.append(round(float(x),2))  
                    ptbinsfloat = sorted(ptbinsfloat)

            print(etabinsfloat)
            print(ptbinsfloat)
            
        histo[k] = ROOT.TH2D(f"muonPOGrecoSF_{k}", f" {k} muon POG reco SF", len(etabinsfloat)-1, array('d', etabinsfloat), len(ptbinsfloat)-1, array('d', ptbinsfloat))

        for eta in data.keys():
            etaRange = eta.split(":")[1].replace("[","").replace("]","")
            etalow, etahigh = map(float, etaRange.split(","))
            etacenter = 0.5 * (etahigh + etalow)
            etabin = histo[k].GetXaxis().FindFixBin(etacenter)
            for pt in data[eta].keys():
                ptRange = pt.split(":")[1].replace("[","").replace("]","")
                ptlow, pthigh = map(float, ptRange.split(","))
                ptcenter = 0.5 * (pthigh + ptlow)
                ptbin = histo[k].GetYaxis().FindFixBin(ptcenter)
                histo[k].SetBinContent(etabin, ptbin, data[eta][pt]["value"])
                histo[k].SetBinError(etabin, ptbin, data[eta][pt]["error"])
            
        drawCorrelationPlot(histo[k], "muon #eta", "muon p_{T} (GeV)", f"data/MC scale factor::0.95,1.05",
                            f"muonPOGrecoSF_{k}", plotLabel="ForceTitle", outdir=plotFolder, passCanvas=canvas)
        drawCorrelationPlot(histo[k], "muon #eta", "muon p_{T} (GeV)", f"absolute uncertainty on SF::0.0,1.00",
                            f"absoluteUncertainty_muonPOGrecoSF_{k}", plotLabel="ForceTitle", outdir=plotFolder, passCanvas=canvas,
                            plotError=True, drawOption="COLZ0")
