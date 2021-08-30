#!/usr/bin/env python3

# python w-mass-13TeV/makeFakeValidation.py plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion/validationFR_lowMt_lowIso/postVFP_plus/plots_fakerate.root plots/testNanoAOD/WmassPlots_jetEta2p4_fixMassWeight_splitW/fakeRateRegion_postVFP_plus_systTH3/postprocessing/templatesQCD/histogramsQCD.root plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion/validationFR_lowMt_highIso/postVFP_plus/plots_fakerate.root plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion/validationFR_lowMt_lowIso/postVFP_plus/plotWithFakes/

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
#from plotUtils.utility import *

sys.path.append(os.getcwd())
from cropNegativeTemplateBins import cropNegativeContent

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input file with histograms")
    parser.add_argument("rootfileFake", type=str, nargs=1, help="Input file with histograms to make fake rate factor in low mT region")
    parser.add_argument("rootfileHighIso", type=str, nargs=1, help="Input file with histogram in region with high iso and low mT")
    parser.add_argument("outdir",   type=str, nargs=1, help="Ouput folder for plots")
    args = parser.parse_args()
              
    fname = args.rootfile[0]
    fnameFake = args.rootfileFake[0]
    outdir = args.outdir[0]
    if not outdir.endswith('/'):
        outdir += '/'
    createPlotDirAndCopyPhp(outdir)

    ROOT.TH1.SetDefaultSumw2()

    xAxisName = "Muon #eta"
    yAxisName = "Muon p_{T} (GeV)"
    
    f = safeOpenFile(fname)
    hdata = safeGetObject(f, "muon_pt_eta_data")
    hmc = safeGetObject(f, "muon_pt_eta_background")
    f.Close()

    ff = safeOpenFile(fnameFake)
    # hnum = safeGetObject(ff, "yields_dataSubMC_lowIso_lowMt")
    # hden = safeGetObject(ff, "yields_dataSubMC_highIso_lowMt")
    # ff.Close()
    # hfakerateFactor = copy.deepcopy(hnum.Clone("hfakerateFactor"))
    # hfakerateFactor.Divide(hden)
    hfakerateFactor = safeGetObject(ff, "fakerateFactor")
    hfakerateFactor.SetTitle(">=1 jet && m_{T} < 40 GeV")

    canvas = ROOT.TCanvas("canvas", "", 800, 700)

    drawCorrelationPlot(hfakerateFactor,
                        xAxisName,
                        yAxisName,
                        f"fakerate factor",
                        hfakerateFactor.GetName(),
                        plotLabel="ForceTitle",
                        outdir=outdir,
                        passCanvas=canvas,
                        drawOption="COLZ0")

    fhiso = safeOpenFile(fname)
    hdataHighIso = safeGetObject(fhiso, "muon_pt_eta_data")
    hmcHighIso = safeGetObject(fhiso, "muon_pt_eta_background")
    fhiso.Close()

    hdataSubMChighIso = copy.deepcopy(hdataHighIso.Clone("hdataSubMChighIso"))
    hdataSubMChighIso.Add(hmcHighIso, -1.0)
    hdataSubMChighIso.SetTitle("High iso && m_{T} < 40 GeV")
    drawCorrelationPlot(hdataSubMChighIso,
                        xAxisName,
                        yAxisName,
                        f"yields: data - MC",
                        hdataSubMChighIso.GetName(),
                        plotLabel="ForceTitle",
                        outdir=outdir,
                        passCanvas=canvas,
                        drawOption="COLZ0")

    hqcd = copy.deepcopy(hdataSubMChighIso.Clone("data_fakes"))
    hqcd.Multiply(hfakerateFactor)
    hqcd.SetTitle("m_{T} < 40 GeV")
    
    ratio = copy.deepcopy(hdata.Clone("ratio_dataSubMCoverQCD"))
    ratio.Add(hmc, -1.0)
    ratio.Divide(hqcd)    
    ratio.SetTitle("m_{T} < 40 GeV")
    
    drawCorrelationPlot(hqcd,
                        xAxisName,
                        yAxisName,
                        f"QCD template",
                        hqcd.GetName(),
                        plotLabel="ForceTitle",
                        outdir=outdir,
                        passCanvas=canvas,
                        drawOption="COLZ0")
    drawCorrelationPlot(ratio,
                        xAxisName,
                        yAxisName,
                        f"ratio: (data - MC) / QCD::0.8,1.2",
                        ratio.GetName(),
                        plotLabel="ForceTitle",
                        outdir=outdir,
                        passCanvas=canvas,
                        drawOption="COLZ0")

    hqcd.SetMarkerSize(0)
    hmc.SetMarkerSize(0)
    hdata_eta = hdata.ProjectionX("data_eta", 0, -1, "e")
    hmc_eta   = hmc.ProjectionX("mc_eta", 0, -1, "e")
    hqcd_eta   = hqcd.ProjectionX("qcd_eta", 0, -1, "e")
    hdata_pt  = hdata.ProjectionY("data_pt", 0, -1, "e")
    hmc_pt    = hmc.ProjectionY("mc_pt", 0, -1, "e")
    hqcd_pt    = hqcd.ProjectionY("qcd_pt", 0, -1, "e")

    stack_eta = ROOT.THStack("stack_eta", "signal and backgrounds")
    hmc_eta.SetFillColor(ROOT.kRed+2)
    hmc_eta.SetLineColor(ROOT.kBlack)
    hmc_eta.SetMarkerSize(0)
    stack_eta.Add(hmc_eta)
    hqcd_eta.SetFillColor(ROOT.kGray)
    hqcd_eta.SetLineColor(ROOT.kBlack)
    hqcd_eta.SetMarkerSize(0)
    stack_eta.Add(hqcd_eta)
    
    stack_pt = ROOT.THStack("stack_pt", "signal and backgrounds")
    hmc_pt.SetFillColor(ROOT.kRed+2)
    hmc_pt.SetLineColor(ROOT.kBlack)
    hmc_pt.SetMarkerSize(0)
    stack_pt.Add(hmc_pt)
    hqcd_pt.SetFillColor(ROOT.kGray)
    hqcd_pt.SetLineColor(ROOT.kBlack)
    hqcd_pt.SetMarkerSize(0)
    stack_pt.Add(hqcd_pt)
    
    legend = ROOT.TLegend(0.2,0.72,0.95,0.92)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetNColumns(2)
    legend.AddEntry(hdata_eta, "Data", "EP")
    legend.AddEntry(hmc_eta, "MC prompt", "F")
    legend.AddEntry(hqcd_eta, "QCD multijet", "F")

    canvas1D = ROOT.TCanvas("canvas1D", "", 800, 900)

    drawTH1dataMCstack(hdata_eta, stack_eta, xAxisName, "Events", "muon_eta_validationRegion",
                       outdir, legend, ratioPadYaxisNameTmp="Data/MC::0.8,1.2", passCanvas=canvas1D,
                       drawLumiLatex=True, xcmsText=0.3, noLegendRatio=True
    )
    drawTH1dataMCstack(hdata_pt, stack_pt, yAxisName, "Events", "muon_pt_validationRegion",
                       outdir, legend, ratioPadYaxisNameTmp="Data/MC::0.8,1.2", passCanvas=canvas1D,
                       drawLumiLatex=True, xcmsText=0.3, noLegendRatio=True
    )

