#!/usr/bin/env python3

# plot all SF from the input root file, and make some products for faster usage
# can also plot those for data/data and MC/MC if they are in the file (can be added with w-mass-13TeV/makeEffRatioPrePostVFP.py, or using option --make-preOverPost which calls that other script)

# python w-mass-13TeV/plotSF.py testMuonSF/2021-05-31_allSFs_nodz_dxybs.root plots/testNanoAOD/testSF/SFeta0p1_31May2021_nodz_dxybs/globalAndPerEra/ -e 'BtoF,GtoH,B,C,D,E,F,G,H' -n 'trigger,idip,iso,antiiso,isonotrig,antiisonotrig' --makePreOverPost

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

    sf_version = ["nominal", "dataAltSig", "MCAltSig", "dataMCAltSig"]
    
    names = args.sfnames.split(',')
    f = ROOT.TFile.Open(fname)
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    for n in names:
        charges = ["plus", "minus"] if n == "trigger" else ["both"]
        for ch in charges:
            for era in eras:
                for v in sf_version:
                    tmpEra = "F_preVFP" if era == "F" else era
                    hname = f"SF2D_{v}_{n}_{tmpEra}_{ch}"
                    if v == "nominal":
                        hkey = n
                    else:
                        hkey = f"{v}_{n}"
                    hkey += (ch if n == "trigger" else "")
                    #print(f"{hkey} -> {hname}")
                    hists[era][hkey] = f.Get(hname)
                    if hists[era][hkey] is None:
                        raise RuntimeError(f"Error when getting histogram {hname}")
                    # print("type(%s) = %s" % (hname,type(hists[era][hkey])))
                    hists[era][hkey].SetDirectory(0)
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
            prodname = f"fullSF2D_{key}_{era}"
            for i,name in enumerate(productsToMake[key]): 
                if i == 0:
                    stringProduct = name
                    prodHists[era][prodname] = copy.deepcopy(hists[era][name].Clone(prodname))
                else:
                    stringProduct = stringProduct + "*" + name
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
        for n in list(prodHists[era].keys()):
            prodHists[era][n].Write(n)
    f.Close()

                            
    # proceeding to plot these histograms
    canvas = ROOT.TCanvas("canvas","",800,800)

    minmax = {"triggerplus"  : "0.65,1.15",
              "triggerminus" : "0.65,1.15",
              "idip"         : "0.98",
              "iso"          : "0.975,1.025"
    }
    
    for era in eras:

        outdir = outdirOriginal + era + "/"

        for n in list(hists[era].keys()):
            drawCorrelationPlot(hists[era][n], "muon #eta", "muon p_{T} (GeV)", f"{n} data/MC scale factor",
                                f"muonSF_{n}", plotLabel="ForceTitle", outdir=outdir,
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                                nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
            # abs. uncertainty
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
