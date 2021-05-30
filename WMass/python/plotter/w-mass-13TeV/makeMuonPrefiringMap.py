#!/usr/bin/env python3

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

sys.path.append(os.getcwd())
from lumiStuff.runPerEra import lumiForEra_UL_new as lpe

def getHistoFromFile(fname, hname=[]):

    ret = []
    f = ROOT.TFile.Open(fname)
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    for h in hname:
        hist = f.Get(h)
        if not hist:
            raise RuntimeError(f"Could not get histogram {h}")
        hist.SetDirectory(0)
        ret.append(hist)
    f.Close()
    return ret
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="output file name with path (can be in different folder than plots)")
    parser.add_argument("outdir",   type=str, nargs=1, help="output folder name for plots")
    
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()
    
    outdir = args.outdir[0] + "/"
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
        htmlpath = "./templates/index.php"
        shutil.copy(htmlpath, outdir)

    hists = getHistoFromFile("./testMuonSF/muonPrefiring_fineEta.root",
                             ["muonPrefiring_beforeRun2016H", "muonPrefiring_Run2016H"])

    hbefore = hists[0]
    hafter = hists[1]

    hPreVFP = copy.deepcopy(hbefore.Clone("muonPrefiring_preVFP"))
    hPostVFP = copy.deepcopy(hbefore.Clone("muonPrefiring_postVFP"))
    hPostVFP.Reset("ICESM")
    hRunH = copy.deepcopy(hafter.Clone("muonPrefiring_runH"))
    
    # for postVFP need to keep the 8 barrel bins (4 for each side) from hbefore, while for the 8 endcap bins need
    # weighted average of hbefore and hafter (weight given by luminosity in (F_postVFP + G) and (H)
    weight_H = lpe["H"]
    weight_FG = lpe["F_postVFP"] + lpe["G"] 
    hPostVFP.Add(hbefore, hafter, weight_FG, weight_H)
    hPostVFP.Scale(1.0 / (weight_FG + weight_H))
    # set again central bins to hbefore, because the effect in EB was not actually fixed in Run H
    # also set histogram for RunH to fixed in endcap
    for i in range(1, 17): # we know there are 16 bins, 8 for each eta side, let's keep it like this
        if i > 4 and i < 13:
            hPostVFP.SetBinContent(i, hbefore.GetBinContent(i))
            hPostVFP.SetBinError(  i, hbefore.GetBinError(i))
            hRunH.SetBinContent(i, hbefore.GetBinContent(i))
            hRunH.SetBinError(  i, hbefore.GetBinError(i))

    adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas", "", 800, 800)

    drawNTH1([hbefore, hafter], ["Before fix", "After fix"],
             "Muon #eta", "L1 prefiring probability",
             "muonPrefiringProbability_prePostRunH",
             outdir,
             draw_both0_noLog1_onlyLog2=1,
             legendCoords="0.18,0.9,0.82,0.9;2",
             lowerPanelHeight=0,
             passCanvas=canvas,
             lumi=35.9,
             onlyLineColor=True,
             drawErrorAll=True
    )

    
    drawNTH1([hPreVFP, hPostVFP, hRunH], ["preVFP", "postVFP", "H"],
             "Muon #eta", "L1 prefiring probability",
             "muonPrefiringProbability_prePostVFP",
             outdir,
             draw_both0_noLog1_onlyLog2=1,
             legendCoords="0.18,0.9,0.82,0.9;3",
             lowerPanelHeight=0,
             passCanvas=canvas,
             lumi=36.3,
             onlyLineColor=True,
             drawErrorAll=True
    )

    fname = args.rootfile[0]
    f = ROOT.TFile.Open(fname, "RECREATE")
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    hPreVFP.Write()
    hPostVFP.Write()
    hRunH.Write()
    f.Close()
    print(f"Histograms saved in file {fname}")
    print()
