#!/bin/env python

# script to make PU weights from existing PU profiles per era in 2016 UL for Wmass analysis
import os, re, array, math
import argparse
import copy

## safe batch mode                                 
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str, default="plots/testNanoAOD/check_newPUweights_perRunEra/", help="Folder to store plots and weights")
    parser.add_argument("-e", "--eras", type=str, default="B,C,D,E,F,F_postVFP,G,H", help="Eras to be used in plots")
    parser.add_argument("--normalizeIntegralUpToBin", type=int, default=100, help="Use integral up to this bin to normalize PU profiles")
    parser.add_argument("--maxWeightWarning", type=float, default=5.0, help="At the end issue a warning if weight > this value")
    parser.add_argument("--cropHighWeight", type=float, default=5.0, help="If larger than 0, crop values larger than maxWeightWarning, setting them to this value (i.e. use 1 for no reweighting, or same as maxWeightWarning)")
    args = parser.parse_args()

    outdir = args.outdir + "/"
    createPlotDirAndCopyPhp(outdir)

    dataname = "pileup"
    
    datahists = {}
    datahistsRef = {}

    eras = args.eras.split(',')
    refEra = "B" if "B" in eras else "2016"
    for era in eras:
        # file 1
        if era == "2016":
            #datafile = "pileupStuff/pileupProfileData_2016Legacy_TESTFULL2016_99bins.root"
            datafile = "pileupStuff/pileupProfileData_2016Legacy_TESTFULL2016_99bins_customJsonNanoFilter.root"
        else:
            datafile = f"pileupStuff/pileupProfileData_2016Legacy_Run{era}_04June2021.root"
        tf = ROOT.TFile.Open(datafile)
        datahists[era] = tf.Get(dataname)        
        datahists[era].SetDirectory(0)
        tf.Close()
        # reference
        if era == "2016":
           datafile = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2016-69200ub-99bins.root"
        else:
           datafile = f"pileupStuff/pileupProfileData_2016Legacy_Run{era}.root"
        tf = ROOT.TFile.Open(datafile)
        datahistsRef[era] = tf.Get(dataname)        
        datahistsRef[era].SetDirectory(0)
        tf.Close()

    adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas", "", 900, 900)

    hists = [copy.deepcopy(datahists[era].Clone(f"tmp{era}")) for era in eras]
    legEntries = [era for era in eras]
    drawNTH1(hists,legEntries, "number of true interactions::0,70", "Normalized events", f"PU_compareErasData_noScaling",
             outdir,
             lowerPanelHeight=0.0, legendCoords="0.6,0.94,0.7,0.9;2", passCanvas=canvas,
             skipLumi=True, onlyLineColor=True, drawErrorAll=True)
    for h in hists:
        h.Scale(datahists[refEra].Integral()/h.Integral())
    refLabelRatio = "ref" if refEra == "2016" else refEra
    drawNTH1(hists,legEntries, "number of true interactions::0,70", "Events", f"PU_compareErasData",
             outdir,
             labelRatioTmp=f"X/{refLabelRatio}::0.0,3.0", legendCoords="0.6,0.94,0.7,0.9;2", passCanvas=canvas,
             skipLumi=True, onlyLineColor=True, drawErrorAll=True)
    
    for era in eras:
        print(f"Checking era {era}")
        outsubdir = args.outdir + "/Run{era}/".format(era=era)
        createPlotDirAndCopyPhp(outsubdir)

        datahist = datahists[era]
        datahistRef = datahistsRef[era]
        if not datahist:
            print(f"Error getting data histogram for era {era}")
            quit()
        if datahist.GetNbinsX() != datahistRef.GetNbinsX():
            print("Warning: histograms have different number of bins")
            quit()
        nbins = datahist.GetNbinsX()
        if datahist.GetXaxis().GetBinLowEdge(1) != datahistRef.GetXaxis().GetBinLowEdge(1):
            print("Warning: histograms have different low edge")
            quit()
        if datahist.GetXaxis().GetBinLowEdge(1+nbins) != datahistRef.GetXaxis().GetBinLowEdge(1+nbins):
            print("Warning: histograms have different up edge")
            quit()

        puWeights = []
            
        datahistRef.Scale(datahist.Integral(0, args.normalizeIntegralUpToBin) / datahistRef.Integral(0, args.normalizeIntegralUpToBin))
        ratio = datahist.Clone(f"puWeight_2016_UL_Run{era}")
        ratio.Reset("ICESM")
        for i in range(1, 1 + nbins):
            puratio = 1.0
            if datahistRef.GetBinContent(i) == 0:
                ratio.SetBinContent(i, puratio)
            else:
                puratio = datahist.GetBinContent(i)/datahistRef.GetBinContent(i)
                ratio.SetBinContent(i,puratio)
            puWeights.append(puratio)

        fileWgt = outsubdir + "/puWeights_Run{era}.txt".format(era=era)
        outf = open(fileWgt, "w")
        outf.write("")
        outf.write("Printing PU weights for Run %s %s\n" % (era, f"(cropped to {args.cropHighWeight} if > {args.maxWeightWarning})" if args.cropHighWeight > 0 else ""))
        outf.write("\n\n")
        for i in range(0,nbins):
            if puWeights[i] > args.maxWeightWarning:
                outf.write("nPU = %d -> weight = %.3f %s\n" % (i+1,puWeights[i], f"cropping to {args.cropHighWeight}" if args.cropHighWeight > 0 else ""))
                if args.cropHighWeight > 0:
                    puWeights[i] = args.cropHighWeight
        outf.write("\n\n")
        outf.write("[{wgt}]\n".format(wgt=", ".join(str(x) for x in puWeights)))
        outf.close()
        print(f"PU weights for {era} saved in {fileWgt}")

        hists = [datahist, datahistRef]
        legEntries = [f"Data {era}", "Reference" if era == "2016" else "Old reference"]
        drawNTH1(hists, legEntries, "number of true interactions", "Normalized events", f"PU_dataVsRef_Run{era}",
                 outsubdir,
                 labelRatioTmp="Data/Ref::0.0,2.0", legendCoords="0.6,0.9,0.76,0.9", passCanvas=canvas,
                 skipLumi=True, drawErrorAll=True)
