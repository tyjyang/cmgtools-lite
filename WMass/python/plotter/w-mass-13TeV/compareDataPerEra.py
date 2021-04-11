#!/usr/bin/env python3

# compare data histograms per era from these files
# plots/testNanoAOD/testMuonPrefire/Zmumu_plus_groupEra_ReReco_vsB/plots_test.root
# plots/testNanoAOD/testMuonPrefire/Zmumu_plus_groupEra_UltraLegacy_vsB/plots_test.root 

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

logging.basicConfig(level=logging.INFO)

def getFromFile(fname, mydict, plots=[], procs=[]):
    f = ROOT.TFile.Open(fname)
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    for pl in plots:
        for p in procs:
            name = f"{pl}_{p}"
            keyname = "B" if p == "data" else p.replace("data_","")
            mydict[(pl,keyname)] = f.Get(name)
            if mydict[(pl,keyname)] == None:
                logging.info(" Cannot find {name} in file {fname}")
                quit()
            mydict[(pl,keyname)].SetDirectory(0)
    f.Close()
    return mydict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfileUL", type=str, nargs=1)
    parser.add_argument("rootfileRR", type=str, nargs=1)
    parser.add_argument("outdir",   type=str, nargs=1)
    parser.add_argument("-p", "--procs", default="data,data_CD,data_EF,data_FG_postVFP,data_H", type=str, nargs=1)    
    parser.add_argument("--plots", default="muon_eta_fine,muon_eta_other_fine,muon_pt,muon_pt_other,muon_iso04,muon_iso04_other,zmass,muon_dz", type=str, nargs=1)    
    args = parser.parse_args()

    outdir = args.outdir[0] + "/"
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
        htmlpath = "./templates/index.php"
        shutil.copy(htmlpath, outdir)

        
    histoUL = {}
    histoRR = {}

    histoUL = getFromFile(args.rootfileUL[0], histoUL, plots=args.plots.split(','), procs=args.procs.split(','))
    histoRR = getFromFile(args.rootfileRR[0], histoRR, plots=args.plots.split(','), procs=args.procs.split(','))

    adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas","",800,800)
    for k in list(histoUL.keys()):
        hists = [histoUL[k], histoRR[k]]
        legEntries = [f"{k[1]} UltraLegacy", f"{k[1]} ReReco"]
        # note: for 2 histograms h[0] is the numerator, while for more it becomes denominator
        drawNTH1(hists, legEntries,
                 hists[0].GetXaxis().GetTitle(), "Events (normalized to Run B)",
                 hists[0].GetName().split("data")[0] + "_" + str(k[1]),
                 outdir,
                 draw_both0_noLog1_onlyLog2=1,
                 labelRatioTmp="UL/RR::0.95,1.05",
                 legendCoords="0.12,0.92,0.84,0.9;2",
                 passCanvas=canvas,
                 skipLumi=True
        )
