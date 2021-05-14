#!/usr/bin/env python3

# plot ratios of histogram variation over nominal
# example
# python w-mass-13TeV/makeSystRatios.py cards/wmass_fixMassWeights/Wmunu_plus_shapes.root plots/testNanoAOD/WmassPlots_jetEta2p4_fixMassWeight/fakeRateRegion_postVFP_plus_systTH3/postprocessing/testHistoFits/ -p "Zmumu,Wmunu,data_fakes" -s ".*muonL1Prefire(1|3|16)Up|.*pdf12|.*muRmuF9PlusUp|.*muRUp|.*effStatTnP(1|112|500)PlusUp|.*lumi.*"

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input root file")
    parser.add_argument("outdir",   type=str, nargs=1, help="Folder for plots")
    parser.add_argument("-s", "--systematics",    type=str, default=".*pdf.*", help="Comma separated list of regular expressions to select systematics to make ratios with nominal")
    parser.add_argument("-p", "--processes",    type=str, default="Wmunu_plus", help="Comma separated list of processes to plot (full name please)")
    parser.add_argument(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (let it be odd, so the central strip is white if not using --abs-value and the range is symmetric)')
    parser.add_argument(     '--palette'  , dest='palette',      default=55, type=int, help='Set palette: use 0 for a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette' , default=False , action='store_true',   help='Inverte color ordering in palette')
    args = parser.parse_args()

    fname = args.rootfile[0]
    outdir = args.outdir[0] + "/"
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
        htmlpath = "./templates/index.php"
        shutil.copy(htmlpath, outdir)

    processes = args.processes.split(',')
    regexp_syst = re.compile(args.systematics.replace(',','|'))

    nominals = {p : None for p in processes}
    print(nominals)
    ratios = {p : [] for p in processes}

    canvas = ROOT.TCanvas("canvas","",900,800) 
    
    f = ROOT.TFile.Open(fname)
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    # get nominals
    for p in processes:
        nominals[p] = f.Get(f"x_{p}")
        if not nominals[p]:
            print(f"Error getting nomnial histogram for process {p}")
            quit()
        else:
            nominals[p].SetDirectory(0)
    # now make a full loop for systematics
    for k in f.GetListOfKeys():
        name = k.GetName()
        if not regexp_syst.match(name): continue
        tokens = name.split("_") # remove "x_" and name of nuisance
        #print(tokens)
        pname = "_".join(tokens[1:-1])
        if pname not in processes: continue
        sname = tokens[-1]
        ratio = f.Get(name)
        ratio.SetDirectory(0)
        #ratios[pname].append(htmp.Divide(nominals[pname]))
        ratio.Divide(nominals[pname])
        ratio.SetTitle(f"syst: {sname}")
        drawCorrelationPlot(ratio, "Muon #eta", "Muon p_{T} (GeV)", "syst / nominal",
                            name, plotLabel="ForceTitle", outdir=outdir,
                            smoothPlot=False, drawProfileX=False, scaleToUnitArea=False, draw_both0_noLog1_onlyLog2=1,
                            palette=args.palette, nContours=args.nContours, invertePalette=args.invertePalette
                            passCanvas=canvas, drawOption="COLZ0")

        
    print()
