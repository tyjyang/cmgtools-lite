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
from cropNegativeTemplateBins import cropNegativeContent

logging.basicConfig(level=logging.INFO)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input file with histograms")
    parser.add_argument("outdir",   type=str, nargs=1, help="Ouput folder for plots")
    parser.add_argument("-b", "--ptBins", required=True, type=str, help = "Binning for pt in a single region, formatted as in TH1 constructor, e.g. '29,26,55' (uniform binning expected for now)")
    args = parser.parse_args()

    if len(args.ptBins.split(',')) != 3:
        print("Error: the pt binning is expected to be of the form 'nPt,ptLow,ptHigh'. Abort")
        quit()
              
    fname = args.rootfile[0]
    outdir = args.outdir[0]
    if not outdir.endswith('/'):
        outdir += '/'
    createPlotDirAndCopyPhp(outdir)

    ROOT.TH1.SetDefaultSumw2()
    
    f = ROOT.TFile.Open(fname)
    hists = {}
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    for k in f.GetListOfKeys():
        name = k.GetName()  
        if not name.startswith("nominal_"):
            continue
        proc = name.split("__")[1] 
        hists[proc] = f.Get(name)
        hists[proc].SetDirectory(0)
    f.Close()

    # get data-MC for all regions
    hDataSubMC = hists["data"].Clone("dataSubMC")
    for k in hists.keys():
        if k == "data": continue
        hDataSubMC.Add(hists[k],-1.0)

    # now unpack the four regions within the histogram
    regionId = {0: "highIso_lowMt",
                1: "lowIso_lowMt",
                2: "highIso_highMt",
                3: "lowIso_highMt"}

    hFakes = {}
    nEtaBins = hDataSubMC.GetNbinsX()
    etaLow   = hDataSubMC.GetXaxis().GetBinLowEdge(1)
    etaHigh  = hDataSubMC.GetXaxis().GetBinLowEdge(1+nEtaBins)
    tokens = args.ptBins.split(',')
    nPtBins = int(tokens[0])
    ptLow   = float(tokens[1])
    ptHigh  = float(tokens[2])

    xAxisName = "Muon #eta"
    yAxisName = "Muon p_{T} (GeV)"
    
    canvas = ROOT.TCanvas("canvas", "", 800, 700)
    adjustSettings_CMS_lumi()
    
    for k in regionId.keys():
        hFakes[k] = ROOT.TH2D(f"hFakes_{regionId[k]}", regionId[k],
                              nEtaBins, etaLow, etaHigh,
                              nPtBins,  ptLow,  ptHigh)
        for ieta in range(1, 1 + nEtaBins):
            for ipt in range(1, 1 + nPtBins):
                content = hDataSubMC.GetBinContent(ieta, ipt + k * nPtBins)
                error   = hDataSubMC.GetBinError(  ieta, ipt + k * nPtBins)
                hFakes[k].SetBinContent(ieta, ipt, content)
                hFakes[k].SetBinError(  ieta, ipt, error)

        # check there are no negative bins
        wasCropped = cropNegativeContent(hFakes[k], silent=False, cropError=False)
        if wasCropped:
            print(f"Histogram cropped to 0 in {regionId[k]} region")
                
        drawCorrelationPlot(hFakes[k], hists["data"].GetXaxis().GetTitle(), hists["data"].GetYaxis().GetTitle(), "Events",
                            f"yields_dataSubMC_{regionId[k]}", plotLabel="ForceTitle", outdir=outdir,
                            passCanvas=canvas, drawOption="COLZ0")


    templateQCD = copy.deepcopy(hFakes[1].Clone("templateQCD"))
    templateQCD.Divide(hFakes[0])
    templateQCD.Multiply(hFakes[2])
    templateQCD.SetTitle("QCD template")
    drawCorrelationPlot(templateQCD, hists["data"].GetXaxis().GetTitle(), hists["data"].GetYaxis().GetTitle(), "Events",
                        f"yields_templateQCD", plotLabel="ForceTitle", outdir=outdir,
                        passCanvas=canvas, drawOption="COLZ0")
    drawCorrelationPlot(templateQCD,
                        hists["data"].GetXaxis().GetTitle(), hists["data"].GetYaxis().GetTitle(), "Absolute uncertainty",
                        f"absUncertainty_templateQCD", plotLabel="ForceTitle", outdir=outdir,
                        passCanvas=canvas, plotError=True, drawOption="COLZ0")
    drawCorrelationPlot(templateQCD,
                        hists["data"].GetXaxis().GetTitle(), hists["data"].GetYaxis().GetTitle(), "Relative uncertainty",
                        f"relUncertainty_templateQCD", plotLabel="ForceTitle", outdir=outdir,
                        passCanvas=canvas, plotRelativeError=True, drawOption="COLZ0")
    
