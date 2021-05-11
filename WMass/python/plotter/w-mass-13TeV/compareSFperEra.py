#!/usr/bin/env python3

# compare SF per era (mainly in preVFP) with the inclusive one. Also check the luminosity weighted average of the SF vs the inclusive

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

sys.path.append(os.getcwd() + "/w-mass-13TeV/")
from makeEffRatioPrePostVFP import getAntiisoEfficiency
sys.path.append(os.getcwd())
from plotUtils.utility import *
from lumiStuff.runPerEra import lumiForEra_UL_new as lpe

logging.basicConfig(level=logging.INFO)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1)
    parser.add_argument("outdir",   type=str, nargs=1)
    parser.add_argument("-e", "--era",    type=str, default="BtoF", choices=["BtoF"], help="Comma separated list of eras for SF in histogram name; default: %(default)s")
    parser.add_argument("-n", "--sfnames", type=str, default="trigger,idip,iso,isonotrig", help="Comma separated list of efficiency names inside root file, which will be used (trigger uses both plus and minus automatically); default: %(default)s, (antiiso,antiisonotrig have to be made, they are not in the file)")
    parser.add_argument("--sub-era", dest="subEra", type=str, default="", help="If given, comma-separated list of eras to use, others in main era will be ignored (including the inclusive one)")
    args = parser.parse_args()

    eraVFP = args.era
    fname = args.rootfile[0]

    outdir = args.outdir[0]
    if not outdir.endswith('/'):
        outdir = outdir + '/'
    createPlotDirAndCopyPhp(outdir)

    ROOT.TH1.SetDefaultSumw2()

        
    productsToMake = {"isoTrigPlus"       : ["iso",           "triggerplus",  "idip"], # "tracking"],
                      "isoTrigMinus"      : ["iso",           "triggerminus", "idip"], # "tracking"],
                      "isoNotrig"         : ["isonotrig",                     "idip"], # "tracking"],
                      "noisoTrigPlus"     : [                 "triggerplus",  "idip"], # "tracking"],
                      "noisoTrigMinus"    : [                 "triggerminus", "idip"], # "tracking"],
                      "noisoNotrig"       : [                                 "idip"], # "tracking"],
                      "antiisoTrigPlus"   : ["antiiso",       "triggerplus",  "idip"], # "tracking"],
                      "antiisoTrigMinus"  : ["antiiso",       "triggerminus", "idip"], # "tracking"],
                      "antiisoNotrig"     : ["antiisonotrig",                 "idip"], # "tracking"],
    }

    hists = {}
    histsMC = {}
    #eras = ["B", "C", "D", "E", "F", "BtoF", "G", "H"] if eraVFP == "BtoF" else ["G", "H", "GtoH"]
    eras = ["B", "C", "D", "E", "F", "BtoF"] if eraVFP == "BtoF" else ["G", "H", "GtoH"]
    if len(args.subEra):
        eras = [str(x) for x in args.subEra.split(',')]
    for era in eras:
        hists[era] = {}
        histsMC[era] = {}
        
    lumi = 19.5 if eraVFP == "BtoF" else 16.8
        
    effKeys = []
    names = args.sfnames.split(',')
    f = ROOT.TFile.Open(fname)
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    for n in names:
        charges = ["plus", "minus"] if n == "trigger" else ["both"]
        for ch in charges:
            for era in eras:
                actualEra = era if era != "F" else "F_preVFP"
                #print(hname)
                hkey = n + (ch if n == "trigger" else "")
                # first data
                hname = f"effData_{n}_{actualEra}_{ch}"
                hists[era][hkey] = f.Get(hname)
                if hists[era][hkey] is None:
                    raise RuntimeError(f"Error when getting histogram {hname}")
                hists[era][hkey].SetDirectory(0)
                # now MC
                hnameMC = f"effMC_{n}_{actualEra}_{ch}"
                histsMC[era][hkey] = f.Get(hnameMC)
                if histsMC[era][hkey] is None:
                    raise RuntimeError(f"Error when getting histogram {hnameMC}")
                histsMC[era][hkey].SetDirectory(0)
                if hkey not in effKeys:
                    effKeys.append(hkey)
    f.Close()

    addAntiIso = False
    addAntiIsonotrig = False
    if "antiiso" not in effKeys:
        addAntiIso = True
    if "antiisonotrig" not in effKeys:
        addAntiIsonotrig = True
        
    for era in eras:
        if addAntiIso:
            hists[era]["antiiso"] = getAntiisoEfficiency(hists[era]["iso"], hists[era]["iso"].GetName().replace("_iso_","_antiiso_"))
            histsMC[era]["antiiso"] = getAntiisoEfficiency(histsMC[era]["iso"], histsMC[era]["iso"].GetName().replace("_iso_","_antiiso_"))
            effKeys.append("antiiso")
        
        if addAntiIsonotrig:
            hists[era]["antiisonotrig"] = getAntiisoEfficiency(hists[era]["isonotrig"], hists[era]["isonotrig"].GetName().replace("_isonotrig_","_antiisonotrig_"))
            histsMC[era]["antiisonotrig"] = getAntiisoEfficiency(histsMC[era]["isonotrig"], histsMC[era]["isonotrig"].GetName().replace("_isonotrig_","_antiisonotrig_"))
            effKeys.append("antiisonotrig")

    print(effKeys)

        
    prodHists = {}
    for era in eras:
        prodHists[era] = {}
    prodHistsMC = {}
    for era in eras:
        prodHistsMC[era] = {}

    for key in productsToMake.keys():
        for era in eras:
            for i,name in enumerate(productsToMake[key]): 
                if i == 0:
                    stringProduct = name
                    prodHists[era][key] = copy.deepcopy(hists[era][name].Clone(key))
                    prodHistsMC[era][key] = copy.deepcopy(histsMC[era][name].Clone(key))
                else:
                    stringProduct = stringProduct + "*" + name
                    if not prodHists[era][key].Multiply(hists[era][name]):
                        print(f"ERROR in multiplication for prodHists[{era}][{key}] with {name}")
                        quit()
                    if not prodHistsMC[era][key].Multiply(histsMC[era][name]):
                        print(f"ERROR in multiplication for prodHistsMC[{era}][{key}] with {name}")
                        quit()
            prodHists[era][key].SetTitle(f"{stringProduct}")            
            prodHistsMC[era][key].SetTitle(f"{stringProduct}")            


    # proceeding to plot these histograms
    canvas = ROOT.TCanvas("canvas","",800,700)
    adjustSettings_CMS_lumi()
    canvas1D = ROOT.TCanvas("canvas1D","",700,700)
    adjustSettings_CMS_lumi()
    
    for n in effKeys:

        outdir_n = outdir + n + "/"
        createPlotDirAndCopyPhp(outdir_n)

        if eraVFP in eras:
            lumiAverage = copy.deepcopy(hists[eraVFP][n].Clone(f"lumiAveEffData_{n}_{eraVFP}"))
            lumiAverage.SetTitle("lumi-weighted average efficiency")
            lumiAverage.Reset("ICESM")

            for era in eras:
                if era == eraVFP or (eraVFP == "BtoF" and any(era == x for x in ["G", "H"])):
                    continue
                lumiAverage.Add(hists[era][n], lpe[era])
            lumiAverage.Scale(1./lumi)
            drawCorrelationPlot(lumiAverage, "muon #eta", "muon p_{T} (GeV)", f"{n} data efficiency",
                                lumiAverage.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                                passCanvas=canvas)

            pull = ROOT.TH1D(f"pull_fullEffData_lumiAverageVsInclusive_{n}_{eraVFP}","(average - inclusive) / uncertainty",100,-5,5)
            for ix in range(1, 1 + lumiAverage.GetNbinsX()):
                for iy in range(1, 1 + lumiAverage.GetNbinsY()):
                    diff = lumiAverage.GetBinContent(ix, iy) - hists[eraVFP][n].GetBinContent(ix, iy)
                    err = math.pow(lumiAverage.GetBinError(ix, iy), 2) + math.pow(hists[eraVFP][n].GetBinError(ix, iy), 2)
                    pull.Fill(diff / math.sqrt(err))
            drawTH1(pull, "pull", "Events", outdir_n, prefix=pull.GetName(),passCanvas=canvas1D)

            ratioToInclusive = copy.deepcopy(lumiAverage.Clone(f"ratio_{lumiAverage.GetName()}"))
            ratioToInclusive.Divide(hists[eraVFP][n])
            ratioToInclusive.SetTitle("lumi-weighted average / inclusive")
            drawCorrelationPlot(ratioToInclusive, "muon #eta", "muon p_{T} (GeV)", f"{n} data efficiency ratio",
                                ratioToInclusive.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                                passCanvas=canvas)

        for era in eras:
            if era == eraVFP or era == "B":
                continue
            #if "trigger" not in n:
            #    continue
            ratioToB = copy.deepcopy(hists["B"][n].Clone(f"ratio_Bover{era}_effData_{n}"))
            ratioToB.Divide(hists[era][n])
            ratioToB.SetTitle(f"{hists[era][n].GetTitle()}: B/{era}")
            drawCorrelationPlot(ratioToB, "muon #eta", "muon p_{T} (GeV)", f"B/{era} data efficiency ratio",
                                ratioToB.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                                passCanvas=canvas)
        
    # now for product, repeating the same code
    # also open root file to save some histograms
    fout = ROOT.TFile.Open(outdir + "productEffAndSFperEra.root", "RECREATE")        
    if not fout or not fout.IsOpen():
        raise RuntimeError(f"Error when opening file {fout.name}")
    fout.cd()

    for n in prodHists[eras[0]].keys():

        nMC = n.replace("fullEffData_","fullEffMC_") # I think this is no longer needed for the products
        outdir_n = outdir + f"products/{n}/"
        createPlotDirAndCopyPhp(outdir_n)

        if eraVFP in eras:

            lumiAverage = copy.deepcopy(prodHists[eraVFP][n].Clone(f"lumiAveEffData_{n}_{eraVFP}"))
            lumiAverage.SetTitle("lumi-weighted average efficiency")
            lumiAverage.Reset("ICESM")

            lumiAverageMC = copy.deepcopy(prodHistsMC[eraVFP][nMC].Clone(f"lumiAveEffMC_{n}_{eraVFP}"))
            lumiAverageMC.SetTitle("lumi-weighted average efficiency")
            lumiAverageMC.Reset("ICESM")

            for era in eras:
                if era == eraVFP or (eraVFP == "BtoF" and any(era == x for x in ["G", "H"])):
                    continue
                lumiAverage.Add(prodHists[era][n], lpe[era])
                lumiAverageMC.Add(prodHistsMC[era][nMC], lpe[era])
            lumiAverage.Scale(1./lumi)
            lumiAverageMC.Scale(1./lumi)
            drawCorrelationPlot(lumiAverage, "muon #eta", "muon p_{T} (GeV)", f"{n} data efficiency (product)",
                                lumiAverage.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                                passCanvas=canvas)
            drawCorrelationPlot(lumiAverageMC, "muon #eta", "muon p_{T} (GeV)", f"{n} MC efficiency (product)",
                                lumiAverageMC.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                                passCanvas=canvas)

            # data
            pull = ROOT.TH1D(f"pull_fullEffData_lumiAverageVsInclusive_{n}_{eraVFP}","(average - inclusive) / uncertainty",100,-5,5)
            for ix in range(1, 1 + lumiAverage.GetNbinsX()):
                for iy in range(1, 1 + lumiAverage.GetNbinsY()):
                    diff = lumiAverage.GetBinContent(ix, iy) - prodHists[eraVFP][n].GetBinContent(ix, iy)
                    err = math.pow(lumiAverage.GetBinError(ix, iy), 2) + math.pow(prodHists[eraVFP][n].GetBinError(ix, iy), 2)
                    pull.Fill(diff / math.sqrt(err))
            drawTH1(pull, "pull", "Events", outdir_n, prefix=pull.GetName(),passCanvas=canvas1D)

            ratioToInclusive = copy.deepcopy(lumiAverage.Clone(f"ratio_{lumiAverage.GetName()}"))
            ratioToInclusive.Divide(prodHists[eraVFP][n])
            ratioToInclusive.SetTitle("lumi-weighted average / inclusive")
            drawCorrelationPlot(ratioToInclusive, "muon #eta", "muon p_{T} (GeV)", f"{n} data efficiency ratio",
                                ratioToInclusive.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                                passCanvas=canvas)
            #MC
            pullMC = ROOT.TH1D(f"pull_fullEffMC_lumiAverageVsInclusive_{n}_{eraVFP}","(average - inclusive) / uncertainty",100,-5,5)
            for ix in range(1, 1 + lumiAverageMC.GetNbinsX()):
                for iy in range(1, 1 + lumiAverageMC.GetNbinsY()):
                    diff = lumiAverageMC.GetBinContent(ix, iy) - prodHistsMC[eraVFP][nMC].GetBinContent(ix, iy)
                    err = math.pow(lumiAverageMC.GetBinError(ix, iy), 2) + math.pow(prodHistsMC[eraVFP][nMC].GetBinError(ix, iy), 2)
                    pullMC.Fill(diff / (math.sqrt(err) if err != 0.0 else 1.0))
            drawTH1(pullMC, "pull", "Events", outdir_n, prefix=pullMC.GetName(),passCanvas=canvas1D)

            ratioToInclusiveMC = copy.deepcopy(lumiAverageMC.Clone(f"ratio_{lumiAverageMC.GetName()}"))
            ratioToInclusiveMC.Divide(prodHistsMC[eraVFP][nMC])
            ratioToInclusiveMC.SetTitle("lumi-weighted average / inclusive")
            drawCorrelationPlot(ratioToInclusiveMC, "muon #eta", "muon p_{T} (GeV)", f"{n} MC efficiency ratio",
                                ratioToInclusiveMC.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                                passCanvas=canvas)

        # prepare to make lumi-weighted average of SF
        lumiAverageDataMCSF = copy.deepcopy(prodHists[eraVFP][n].Clone(f"lumiAveScaleFactor_{n}_{eraVFP}"))
        lumiAverageDataMCSF.SetTitle("lumi-weighted average data/MC SF")
        lumiAverageDataMCSF.Reset("ICESM")

        inclusiveDataMCSF = None # to keep this one for ratio with luminosity-weighted average
        
        for era in eras:

            sfRatio = copy.deepcopy(prodHists[era][n].Clone(f"fullSF2D_{n}_{era}"))
            sfRatio.Divide(prodHistsMC[era][nMC])
            sfRatio.Write(sfRatio.GetName())
            # use it to compute weighted average
            if era == eraVFP or (eraVFP == "BtoF" and any(era == x for x in ["G", "H"])):
                inclusiveDataMCSF = copy.deepcopy(sfRatio.Clone(f"inclusiveDataMCSF_{eraVFP}"))
            else:
                lumiAverageDataMCSF.Add(sfRatio, lpe[era])

            drawCorrelationPlot(sfRatio, "muon #eta", "muon p_{T} (GeV)", f"Global data/MC scale factor ({era})",
                                sfRatio.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                                passCanvas=canvas)


            
            if era == eraVFP or era == "B":
                continue
            ratioToB = copy.deepcopy(prodHists["B"][n].Clone(f"ratio_Bover{era}_fullEffData_{n}"))
            ratioToB.Divide(prodHists[era][n])
            ratioToB.SetTitle(f"{prodHists[era][n].GetTitle()}: B/{era}")
            ratioToB.Write(ratioToB.GetName())
            drawCorrelationPlot(ratioToB, "muon #eta", "muon p_{T} (GeV)", f"B/{era} full data efficiency ratio",
                                ratioToB.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                                passCanvas=canvas)

            ratioToBinMC = copy.deepcopy(prodHistsMC["B"][nMC].Clone(f"ratio_Bover{era}_fullEffMC_{n}"))
            ratioToBinMC.Divide(prodHistsMC[era][nMC])
            ratioToBinMC.SetTitle(f"{prodHists[era][nMC].GetTitle()}: B/{era}")
            ratioToBinMC.Write(ratioToBinMC.GetName())
            drawCorrelationPlot(ratioToBinMC, "muon #eta", "muon p_{T} (GeV)", f"B/{era} full MC efficiency ratio",
                                ratioToBinMC.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                                passCanvas=canvas)
            
        # loop on eras is over, normalize average SF and plot it, also with pull and ratio to inclusive
        lumiAverageDataMCSF.Scale(1./lumi)
        lumiAverageDataMCSF.Write(lumiAverageDataMCSF.GetName())
        drawCorrelationPlot(lumiAverageDataMCSF, "muon #eta", "muon p_{T} (GeV)", f"{n} data/MC SF (product)",
                            lumiAverageDataMCSF.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                            passCanvas=canvas)

        pull = ROOT.TH1D(f"pull_fullDataMCSF_lumiAverageVsInclusive_{n}_{eraVFP}","(average - inclusive) / uncertainty",100,-5,5)
        for ix in range(1, 1 + lumiAverageDataMCSF.GetNbinsX()):
            for iy in range(1, 1 + lumiAverageDataMCSF.GetNbinsY()):
                diff = lumiAverageDataMCSF.GetBinContent(ix, iy) - inclusiveDataMCSF.GetBinContent(ix, iy)
                err = math.pow(lumiAverageDataMCSF.GetBinError(ix, iy), 2) + math.pow(inclusiveDataMCSF.GetBinError(ix, iy), 2)
                pull.Fill(diff / math.sqrt(err))
        drawTH1(pull, "pull", "Events", outdir_n, prefix=pull.GetName(),passCanvas=canvas1D)

        ratioToInclusive = copy.deepcopy(lumiAverageDataMCSF.Clone(f"ratio_{lumiAverageDataMCSF.GetName()}"))
        ratioToInclusive.Divide(inclusiveDataMCSF)
        ratioToInclusive.SetTitle("lumi-weighted average / inclusive")
        drawCorrelationPlot(ratioToInclusive, "muon #eta", "muon p_{T} (GeV)", f"{n} data/MC SF ratio",
                            ratioToInclusive.GetName(), plotLabel="ForceTitle", outdir=outdir_n,
                            passCanvas=canvas)


            
    fout.Close()
