#!/usr/bin/env python3

# plot all SF from the input root file, and make some products for faster usage

# python w-mass-13TeV/plotSF.py testMuonSF/2021-05-31_allSFs_nodz_dxybs.root plots/testNanoAOD/testSF/SFeta0p1_31May2021_nodz_dxybs/globalAndPerEra/ -e 'BtoF,GtoH,B,C,D,E,F,G,H' -n 'trigger,idip,iso,antiiso,isonotrig,antiisonotrig'

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1)
    parser.add_argument("outdir",   type=str, nargs=1)
    parser.add_argument("-e", "--era",    type=str, default="BtoF,GtoH,B,C,D,E,F,G,H", help="Comma separated list of eras for SF in histogram name; default: %(default)s")
    parser.add_argument("-n", "--sfnames", type=str, default="trigger,idip,iso,antiiso,isonotrig,antiisonotrig,tracking,reco", help="Comma separated list of SF names inside root file, which will be plotted (trigger uses both plus and minus automatically); default: %(default)s")
    parser.add_argument("--sf-version", dest="sfversions", type=str, default="nominal,dataAltSig", help="SF versions to plot and to use for the products, usually one would use just nominal and dataAltSig to define the systematic variation; default: %(default)s")
    parser.add_argument(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51')
    parser.add_argument(     '--palette'  , dest='palette',      default=87, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette' , default=False , action='store_true',   help='Inverte color ordering in palette')
    args = parser.parse_args()

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

    #productsToMake = {"trackingReco"      : ["reco", "tracking"], # "tracking"],
    #}
    
    eras = args.era.split(',')
    fname = args.rootfile[0]
    outdirOriginal = args.outdir[0] # to keep it below
    if not outdirOriginal.endswith('/'):
        outdirOriginal = outdirOriginal + '/'
    productSubfolder = "productSF/"
    
    # make folder structure
    foldersToCreate = []
    foldersToCreate.append(outdirOriginal + productSubfolder)
    for era in eras:
        foldersToCreate.append(outdirOriginal + era + "/")
        foldersToCreate.append(outdirOriginal + productSubfolder + era + "/")

    for folder in foldersToCreate:
        createPlotDirAndCopyPhp(folder)
        
    print()
        
    hists = {}
    for era in eras:
        hists[era] = {}
    histsPrePost = {"Data" : {},
                    "MC"   : {}
    }

    sf_version = [str(x) for x in args.sfversions.split(",")]
    
    names = args.sfnames.split(',')
    f = ROOT.TFile.Open(fname)
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    for n in names:
        charges = ["plus", "minus"] if n == "trigger" else ["both"]
        for ch in charges:
            for era in eras:
                tmpEra = "F_preVFP" if era == "F" else era
                # arrays for bad bins in reco sf
                badbinsNomi = []
                badbinsAlt = []
                for v in sf_version:
                    hname = f"SF2D_{v}_{n}_{tmpEra}_{ch}"
                    #if v == "nominal":
                    #    hkey = n
                    #else:
                    #    hkey = f"{v}_{n}"
                    hkey = f"{v}_{n}" 
                    hkey += (ch if n == "trigger" else "")
                    #print(f"{hkey} -> {hname}")
                    hists[era][hkey] = f.Get(hname)
                    if hists[era][hkey] is None:
                        raise RuntimeError(f"Error when getting histogram {hname}")
                    # print("type(%s) = %s" % (hname,type(hists[era][hkey])))
                    hists[era][hkey].SetDirectory(0)
                    # some hardcoded stuff, setting histograms based on tnp plots, if bad we set sf to 0 to flag them out when comparing
                    # there are 48 bins in eta, from 0 to 47
                    # for preVFP see for example https://mdunser.web.cern.ch/mdunser/private/w-mass-13TeV/tnp/results_nodz_dxybs_mcTruth_2021-10-06/efficiencies_BtoF/mu_reco_both/plots/
                    # or https://mdunser.web.cern.ch/mdunser/private/w-mass-13TeV/tnp/results_nodz_dxybs_mcTruth_2021-10-08/efficiencies_BtoF/mu_reco_both/plots/
                    # the latter has gen-reco match with DR=0.1 instead of 0.3
                    if n == "reco":
                        badbins = []
                        if "2021-10-07" in args.rootfile[0]:
                            if era == "BtoF":
                                badbinsNomi = [7, 8, 9, 10, 11, 19, 23, 28, 34, 36, 37, 38, 39, 40]
                                badbinsAlt = [3, 12, 21, 29]
                                if v == "nominal":
                                    badbins = badbinsNomi
                                elif v == "dataAltSig":
                                    badbins = badbinsAlt
                            elif era == "GtoH":
                                badbinsNomi = [7, 8, 9, 10, 11, 13, 19, 30, 34, 35, 36, 37, 38, 39, 40, 41, 42]
                                badbinsAlt = [3, 29, 42, 45]
                                if v == "nominal":
                                    badbins = badbinsNomi
                                elif v == "dataAltSig":
                                    badbins = badbinsAlt
                        elif "2021-10-09" in args.rootfile[0] or "2021-10-12" in args.rootfile[0]:
                            # not implemented yet
                            if era == "BtoF":
                                badbinsNomi = [8, 10, 11, 36, 37, 38, 39]
                                badbinsAlt = [1, 5, 26, 43]
                                if v == "nominal":
                                    badbins = badbinsNomi
                                elif v == "dataAltSig":
                                    badbins = badbinsAlt
                            elif era == "GtoH":
                                badbinsNomi = [7, 8, 9, 10, 11, 13, 34, 36, 37, 38, 39, 40]
                                badbinsAlt = [17, 41, 43] # 43 is a bad fit but it yields numerically the correct efficiency, and the SF would be 1, so let's add it here
                                if v == "nominal":
                                    badbins = badbinsNomi
                                elif v == "dataAltSig":
                                    badbins = badbinsAlt
                        for ib in range(0, hists[era][hkey].GetNbinsX()):
                            if ib in badbins:
                                hists[era][hkey].SetBinContent(ib+1, 1, 1)
                                hists[era][hkey].SetBinError(ib+1, 1, 0)
                            elif hists[era][hkey].GetBinContent(ib+1, 1) > 1.0:
                                hists[era][hkey].SetBinContent(ib+1, 1, 1.0)
                if n == "reco":
                    #print(badbinsNomi)
                    #print(badbinsAlt)
                    recoRegKey = f"regularized_reco"
                    recoRegularized = f"SF2D_regularized_reco_{tmpEra}_{ch}"
                    hists[era][recoRegKey] = copy.deepcopy(hists[era]["nominal_reco"].Clone(recoRegularized))
                    hists[era][recoRegKey].SetDirectory(0)
                    hists[era][recoRegKey].SetTitle(recoRegularized)
                    for ib in range(0, hists[era][recoRegKey].GetNbinsX()):
                        if ib in badbinsNomi:
                            if ib not in badbinsAlt:
                                hists[era][recoRegKey].SetBinContent(ib+1, 1, hists[era]["dataAltSig_reco"].GetBinContent(ib+1, 1))
                                hists[era][recoRegKey].SetBinError(  ib+1, 1, hists[era]["dataAltSig_reco"].GetBinError(  ib+1, 1))
                            else:
                                hists[era][recoRegKey].SetBinContent(ib+1, 1, 1)
                                hists[era][recoRegKey].SetBinError(  ib+1, 1, 0)
                        else:
                            hists[era][recoRegKey].SetBinContent(ib+1, 1, hists[era]["nominal_reco"].GetBinContent(ib+1, 1))
                            hists[era][recoRegKey].SetBinError(  ib+1, 1, hists[era]["nominal_reco"].GetBinError(  ib+1, 1))
    f.Close()
    
    prodHists = {}
    for era in eras:
        prodHists[era] = {}

    namesForCheck = [x for x in names if x != "trigger"]
    if "trigger" in names:
        namesForCheck += ["triggerplus", "triggerminus"]
        
    for key in list(productsToMake.keys()):
        if not all(x in namesForCheck for x in productsToMake[key]):
            print()
            missingFactors = ",".join([x for x in productsToMake[key] if x not in namesForCheck])
            print(f"Warning: skipping product {key} because these factors are missing: {missingFactors}")
            print()
            continue
        for era in eras:
            for sfv in sf_version:
                prodname = f"fullSF2D_{sfv}_{key}_{era}"
                for i,basename in enumerate(productsToMake[key]): 
                    name = f"{sfv}_{basename}"
                    if i == 0:
                        stringProduct = basename
                        prodHists[era][prodname] = copy.deepcopy(hists[era][name].Clone(prodname))
                    else:
                        stringProduct = stringProduct + "*" + basename
                        if not prodHists[era][prodname].Multiply(hists[era][name]):
                            print(f"ERROR in multiplication for prodHists[{era}][{prodname}] with {name}")
                            quit()
                prodHists[era][prodname].SetTitle(f"{stringProduct}")            

    fullout = outdirOriginal + "scaleFactorProduct.root" 
    f = ROOT.TFile.Open(fullout, "RECREATE")
    # put new ones in this files, without copying original histograms
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fullout}")
    f.cd()
    for era in eras:
        # save reco and tracking as they are, if present in the list of histograms
        for tmp in ["reco", "tracking", "alttrack", "altre"]:
            if tmp not in names:
                continue
            for sfv in sf_version:
                hists[era][f"{sfv}_{tmp}"].Write()
        if "regularized_reco" in hists[era].keys():
            hists[era]["regularized_reco"].Write()
        # now the product of scale factors
        for n in list(prodHists[era].keys()):
            prodHists[era][n].Write(n)
    f.Close()

                            
    # proceeding to plot these histograms
    canvas = ROOT.TCanvas("canvas","",800,800)

    minmax = {"triggerplus"  : "0.65,1.15",
              "triggerminus" : "0.65,1.15",
              "idip"         : "0.95,1.01",
              "iso"          : "0.975,1.025",
              "antiiso"      : "0.7,1.25",
              "isonotrig"    : "0.97,1.03",
              "antiisonotrig": "0.7,1.25",
              "tracking"     : "0.99,1.01",
              "alttrack"     : "0.99,1.01",
              "reco"         : "0.9,1.02",
              "altre"        : "0.95,1.05"
    }
    
    for era in eras:

        outdir = outdirOriginal + era + "/"

        for n in list(hists[era].keys()):
            ntmp = str(n.split("_")[1])
            if ntmp in minmax.keys():
                zrange = f"::{minmax[ntmp]}"
            else:
                zrange = ""
                        
            drawCorrelationPlot(hists[era][n], "muon #eta", "muon p_{T} (GeV)", f"{n} data/MC scale factor{zrange}",
                                f"muonSF_{n}", plotLabel="ForceTitle", outdir=outdir,
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                                nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
            # abs. uncertainty (only on nominal, it should be the same for all histograms, hoping the SF and efficiencies were sane in this configuration)
            if "nominal" in n:
                drawCorrelationPlot(hists[era][n], "muon #eta", "muon p_{T} (GeV)", f"abs. uncertainty on {n} SF",
                                    f"absUnc_muonSF_{n}", plotLabel="ForceTitle", outdir=outdir+"absoluteUncertainty/",
                                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                    draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotError=True,
                                    nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
            ## rel. uncertainty
            #drawCorrelationPlot(hists[era][n], "muon #eta", "muon p_{T} (GeV)", f"rel. uncertainty on {n} SF",
            #                    f"relUnc_muonSF_{n}", plotLabel="ForceTitle", outdir=outdir+"relativeUncertainty/",
            #                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
            #                    draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotRelativeError=True)

        outdir = outdirOriginal + productSubfolder + era + "/"

        for n in list(prodHists[era].keys()):
            drawCorrelationPlot(prodHists[era][n], "muon #eta", "muon p_{T} (GeV)", "data/MC scale factor product",
                                f"{n}", plotLabel="ForceTitle", outdir=outdir,
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                                nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
            # plot absolute error
            drawCorrelationPlot(prodHists[era][n], "muon #eta", "muon p_{T} (GeV)", "abs. uncertainty on SF product",
                                f"absUnc_{n}", plotLabel="ForceTitle", outdir=outdir+"absoluteUncertainty/",
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotError=True,
                                nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
            ## plot relative error
            #drawCorrelationPlot(prodHists[era][n], "muon #eta", "muon p_{T} (GeV)", "rel. uncertainty on SF product",
            #                    f"relUnc_{n}", plotLabel="ForceTitle", outdir=outdir+"relativeUncertainty/",
            #                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
            #                    draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotRelativeError=True)
