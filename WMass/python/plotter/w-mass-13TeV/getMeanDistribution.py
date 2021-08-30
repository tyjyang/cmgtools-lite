#!/usr/bin/env python3

# python w-mass-13TeV/getMeanDistribution.py plots/testNanoAOD/WmassPlots_distributionsVsIsoMt/fakeRateRegion_postVFP_plus_highIso_highMt/plots_fakerate.root plots/testNanoAOD/WmassPlots_distributionsVsIsoMt/fakeRateRegion_postVFP_plus_highIso_highMt/getMeanDistribution/

import os, re, array, math
import time
import argparse

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
from cropNegativeTemplateBins import cropNegativeContent


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", type=str, nargs=1)
    parser.add_argument("outputfolder",   type=str, nargs=1)
    #parser.add_argument("--hname", default="mt_pt_eta", help="Root of histogram name inside root file")
    parser.add_argument("-x", "--x-axis-name", dest="xAxisName", default="m_{T} (GeV)", help="x axis name")
    parser.add_argument(      "--hname", dest="hName", default="mt_MET", help="histogram prefix to get from file")
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    #outfolder = args.outputfolder[0]
    outfolder = args.outputfolder[0]
    createPlotDirAndCopyPhp(outfolder)

    xAxisName = args.xAxisName
    
    histo = None
    rfile = safeOpenFile(args.inputfile[0])
    datakey = f"{args.hName}_data"
    histo = safeGetObject(rfile, datakey)
    #data1D = histo.ProjectionZ("data_mt_lowIso", 0, -1, 0, -1, "e")
    #mc1D = []
    for key in rfile.GetListOfKeys():
        name = key.GetName()
        if "data" in name or "background" in name: continue
        if args.hName not in name: continue
        obj = key.ReadObj()
        if obj.ClassName() != "TH1D": continue
        obj.SetDirectory(0)
        histo.Add(obj, -1.0)
        #mc1D.append(obj.ProjectionZ(f"{key}_mt_lowIso", 0, -1, 0, -1, "e"))
    rfile.Close()

    histo.SetFillColor(ROOT.kGray+1)
    drawTH1(histo, xAxisName, "Events", outfolder, histo.GetName(), moreTextLatex="data - MC::0.18,0.8,0.08,0.05")
    
