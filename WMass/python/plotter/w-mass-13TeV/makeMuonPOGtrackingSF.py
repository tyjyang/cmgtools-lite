#!/usr/bin/env python3

# transform the graphs into histograms (uncertainties are symmetric with very good approximation)

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

muonfiles = {"preVFP"  : "testMuonSF/fits_preVFP.root",
             "postVFP" : "testMuonSF/fits_postVFP.root"
}
output_muonfile = "testMuonSF/muonPOGtrackingSF.root"
etabinsfloat = [-2.4, -2.1, -1.6, -1.1, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.1, 1.6, 2.1, 2.4]
plotFolder = "plots/testNanoAOD/testMuonPOGtrackingSF/eff_and_sf/"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    args = parser.parse_args()

    createPlotDirAndCopyPhp(plotFolder)

    histo = {}
    
    for era in muonfiles.keys():
        histo[era] = ROOT.TH1D(f"muonPOGtrackingSF_{era}", f"{era} tracking SF (all tracks) from muon POG", len(etabinsfloat)-1, array('d', etabinsfloat))
        #print(histo[era])
        histo[era].SetDirectory(0)
        infile = safeOpenFile(muonfiles[era], mode="READ")
        graph = safeGetObject(infile, "ratio_eff_eta3_dr030e030_corr", quitOnFail=True, silent=True, detach=False)
        yval = list(graph.GetY())
        #print(yval)
        for i,y in enumerate(yval):
            histo[era].SetBinContent(i+1, y)
            histo[era].SetBinError(i+1, 0.5 * (graph.GetErrorYhigh(i) + graph.GetErrorYlow(i))) # set uncertainty as average between high and low from the graph (they are almost equal)
            
    infile.Close()

    drawNTH1([histo["preVFP"], histo["postVFP"]], ["preVFP", "postVFP"], "muon #eta", "muon POG tracking SF::0.98,1.01", "muonPOGtrackingSF_allTracks", plotFolder,
             labelRatioTmp="pre/post::0.98,1.01", legendCoords="0.15,0.85,0.82,0.9;2", canvasSize="800,800", onlyLineColor=True, drawErrorAll=True, yAxisExtendConstant=1.00)#, setOnlyLineRatio=True)
    
    outfile = safeOpenFile(output_muonfile, mode="RECREATE")
    print()
    print(f"Created file {output_muonfile}")
    print()
    outfile.cd()
    for era in muonfiles.keys():
        histo[era].Write()
    outfile.Close()
