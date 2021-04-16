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
    parser.add_argument("-o", "--outdir", type=str, default="plots/testNanoAOD/PU_weights_perRunEra/", help="Folder to store plots and weights")
    parser.add_argument("-e", "--eras", type=str, default="B,C,D,E,F,F_postVFP,G,H", help="Eras to be used in plots")
    parser.add_argument("--normalizeIntegralUpToBin", type=int, default=100, help="Use integral up to this bin to normalize PU profiles")
    parser.add_argument("--maxWeightWarning", type=float, default=5.0, help="At the end issue a warning if weight > this value")
    parser.add_argument("--cropHighWeight", type=float, default=5.0, help="If larger than 0, crop values larger than maxWeightWarning, setting them to this value (i.e. use 1 for no reweighting, or same as maxWeightWarning)")
    args = parser.parse_args()

    outdir = args.outdir + "/"
    createPlotDirAndCopyPhp(outdir)

    mcfile   = "pileupStuff/MyMCPileupHistogram_2016Legacy_noGenWeights_preAndPostVFP.root"    
    mcname   = "Pileup_nTrueInt_Wmunu_preVFP,Pileup_nTrueInt_Wmunu_postVFP"
    dataname = "pileup"
    
    datahists = {}

    eras = args.eras.split(',')
    for era in eras:
        datafile = f"pileupStuff/pileupProfileData_2016Legacy_Run{era}.root"
        tf = ROOT.TFile.Open(datafile)
        datahists[era] = tf.Get(dataname)        
        datahists[era].SetDirectory(0)
        tf.Close()

    mchist = None
    tf = ROOT.TFile.Open(mcfile)
    for i,mcn in enumerate(mcname.split(',')):
        tmp = tf.Get(mcn)
        if i == 0:
            mchist = copy.deepcopy(tmp.Clone("mchist"))
            mchist.SetDirectory(0)
        else:
            mchist.Add(tmp)
    tf.Close()

    if not mchist:
        print("Error getting MC histograms")
        quit()

    adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas", "", 900, 900)

    hists = [copy.deepcopy(datahists[era].Clone(f"tmp{era}")) for era in eras]
    legEntries = [era for era in eras]
    drawNTH1(hists,legEntries, "number of true interactions::0,70", "Events", f"PU_compareErasData_noScaling",
             outdir,
             lowerPanelHeight=0.0, legendCoords="0.6,0.94,0.7,0.9;2", passCanvas=canvas,
             skipLumi=True, onlyLineColor=True, drawErrorAll=True)
    for h in hists:
        h.Scale(datahists["B"].Integral()/h.Integral())
    drawNTH1(hists,legEntries, "number of true interactions::0,70", "Events", f"PU_compareErasData",
             outdir,
             labelRatioTmp="X/B::0.0,3.0", legendCoords="0.6,0.94,0.7,0.9;2", passCanvas=canvas,
             skipLumi=True, onlyLineColor=True, drawErrorAll=True)
    
    for era in eras:
        print(f"Checking era {era}")
        outsubdir = args.outdir + "/Run{era}/".format(era=era)
        createPlotDirAndCopyPhp(outsubdir)

        datahist = datahists[era]
        if not datahist:
            print(f"Error getting data histogram for era {era}")
            quit()
        if datahist.GetNbinsX() != mchist.GetNbinsX():
            print("Warning: histograms have different number of bins")
            quit()
        nbins = datahist.GetNbinsX()
        if datahist.GetXaxis().GetBinLowEdge(1) != mchist.GetXaxis().GetBinLowEdge(1):
            print("Warning: histograms have different low edge")
            quit()
        if datahist.GetXaxis().GetBinLowEdge(1+nbins) != mchist.GetXaxis().GetBinLowEdge(1+nbins):
            print("Warning: histograms have different up edge")
            quit()

        puWeights = []
            
        mchist.Scale(datahist.Integral(0, args.normalizeIntegralUpToBin) / mchist.Integral(0, args.normalizeIntegralUpToBin))
        ratio = datahist.Clone(f"puWeight_2016_UL_Run{era}")
        ratio.Reset("ICESM")
        for i in range(1, 1 + nbins):
            puratio = 1.0
            if mchist.GetBinContent(i) == 0:
                ratio.SetBinContent(i, puratio)
            else:
                puratio = datahist.GetBinContent(i)/mchist.GetBinContent(i)
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

        hists = [datahist, mchist]
        legEntries = [f"Data {era}", "MC"]
        drawNTH1(hists, legEntries, "number of true interactions", "Events", f"PU_dataMC_Run{era}",
                 outsubdir,
                 labelRatioTmp="Data/MC::0.0,2.0", legendCoords="0.6,0.9,0.76,0.9", passCanvas=canvas,
                 skipLumi=True, drawErrorAll=True)
