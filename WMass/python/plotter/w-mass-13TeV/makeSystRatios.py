#!/usr/bin/env python3

# make SF for data/data and MC/MC from efficiency ratios

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
            print(f"Error getting nomnal histogram for process {p}")
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
                            smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                            draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, drawOption="COLZ0")

        
    print()
