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
    parser.add_argument("-n", "--sfnames", type=str, default="trigger,idip,iso,antiiso,isonotrig,antiisonotrig", help="Comma separated list of SF names inside root file, which will be plotted (trigger uses both plus and minus automatically); default: %(default)s")
    parser.add_argument('--makePreOverPost', action="store_true", help="Make data/data and MC/MC scale factors (in subfolder 'effRatio_preOverPost/')")
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
    prePostSubfolder = "effRatio_preOverPost/" 
    productSubfolder = "productSF/"
    
    # make folder structure
    foldersToCreate = []
    foldersToCreate.append(outdirOriginal + prePostSubfolder)
    foldersToCreate.append(outdirOriginal + productSubfolder)
    foldersToCreate.append(outdirOriginal + productSubfolder + prePostSubfolder)
    for era in eras:
        foldersToCreate.append(outdirOriginal + era + "/")
        foldersToCreate.append(outdirOriginal + productSubfolder + era + "/")

    for folder in foldersToCreate:
        createPlotDirAndCopyPhp(folder)
        
    if args.makePreOverPost:
        command = f"python w-mass-13TeV/makeEffRatioPrePostVFP.py {args.rootfile[0]} {outdirOriginal}{prePostSubfolder}"
        print("Running following command to make preVFP/postVFP scale factors")
        print(command)
        ret = os.system(command)
        if ret:
            print("Problem with command. Abort")
            quit()
        newfname = os.path.basename(args.rootfile[0]).replace(".root","_NEW.root")
        fname = outdirOriginal + prePostSubfolder + newfname
        os.system(f"mv {fname} {outdirOriginal}")
        fname = outdirOriginal + newfname
        print(f"Updated file with data/data and MC/MC scale factors moved in {fname}")
        print()
        
    hists = {}
    for era in eras:
        hists[era] = {}
    histsPrePost = {"Data" : {},
                    "MC"   : {}
    }

    names = args.sfnames.split(',')
    f = ROOT.TFile.Open(fname)
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    for n in names:
        charges = ["plus", "minus"] if n == "trigger" else ["both"]
        for ch in charges:
            for era in eras:
                tmpEra = "F_preVFP" if era == "F" else era
                hname = f"SF2D_{n}_{tmpEra}_{ch}"
                hkey = n + (ch if n == "trigger" else "")
                hists[era][hkey] = f.Get(hname)
                if hists[era][hkey] is None:
                    raise RuntimeError(f"Error when getting histogram {hname}")
                # print("type(%s) = %s" % (hname,type(hists[era][hkey])))
                hists[era][hkey].SetDirectory(0)
            if args.makePreOverPost:
                for dataType in ["Data", "MC"]:
                    hname = f"SF2D_{dataType}_preOverPost_{n}_{ch}"
                    hkey = n + (ch if n == "trigger" else "")
                    histsPrePost[dataType][hkey] = f.Get(hname)
                    if histsPrePost[dataType][hkey] is None:
                        raise RuntimeError(f"Error when getting histogram {hname}")
                    # print("type(%s) = %s" % (hname,type(hists[hkey])))
                    # print(hname)
                    histsPrePost[dataType][hkey].SetDirectory(0)
    f.Close()

    prodHists = {}
    for era in eras:
        prodHists[era] = {}
    prodHistsPrePost = {"Data" : {},
                        "MC"   : {}
    }
    
    for key in list(productsToMake.keys()):
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

        if args.makePreOverPost:
            for dataType in list(prodHistsPrePost.keys()):
                prodnamePrePost = f"fullSF2D_{dataType}_preOverPost_{key}"
                for i,name in enumerate(productsToMake[key]): 
                    if i == 0:
                        stringProduct = name
                        prodHistsPrePost[dataType][prodnamePrePost] = copy.deepcopy(histsPrePost[dataType][name].Clone(prodnamePrePost))
                    else:
                        stringProduct = stringProduct + "*" + name
                        if not prodHistsPrePost[dataType][prodnamePrePost].Multiply(histsPrePost[dataType][name]):
                            print(f"ERROR in multiplication for prodHistsPrePost[{prodnamePrePost}] with {name}")
                            quit()
                prodHistsPrePost[dataType][prodnamePrePost].SetTitle(f"{stringProduct}")

    fullout = outdirOriginal + "scaleFactorProduct.root" 
    f = ROOT.TFile.Open(fullout, "RECREATE")
    # put new ones in this files, without copying original histograms
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fullout}")
    f.cd()
    for era in eras:
        for n in list(prodHists[era].keys()):
            prodHists[era][n].Write(n)
    if args.makePreOverPost:
        for dataType in list(prodHistsPrePost.keys()):
            for n in list(prodHistsPrePost[dataType].keys()):
                prodHistsPrePost[dataType][n].Write(n)
    f.Close()

                            
    # proceeding to plot these histograms
    canvas = ROOT.TCanvas("canvas","",800,800)

    for era in eras:

        outdir = outdirOriginal + era + "/"

        for n in list(hists[era].keys()):
            drawCorrelationPlot(hists[era][n], "muon #eta", "muon p_{T} (GeV)", f"{n} data/MC scale factor",
                                f"muonSF_{n}", plotLabel="ForceTitle", outdir=outdir,
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas)
            # abs. uncertainty
            drawCorrelationPlot(hists[era][n], "muon #eta", "muon p_{T} (GeV)", f"abs. uncertainty on {n} SF",
                                f"absUnc_muonSF_{n}", plotLabel="ForceTitle", outdir=outdir+"absoluteUncertainty/",
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotError=True)
            # rel. uncertainty
            drawCorrelationPlot(hists[era][n], "muon #eta", "muon p_{T} (GeV)", f"rel. uncertainty on {n} SF",
                                f"relUnc_muonSF_{n}", plotLabel="ForceTitle", outdir=outdir+"relativeUncertainty/",
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotRelativeError=True)

        outdir = outdirOriginal + productSubfolder + era + "/"

        for n in list(prodHists[era].keys()):
            drawCorrelationPlot(prodHists[era][n], "muon #eta", "muon p_{T} (GeV)", "data/MC scale factor product",
                                f"{n}", plotLabel="ForceTitle", outdir=outdir,
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas)
            # plot absolute error
            drawCorrelationPlot(prodHists[era][n], "muon #eta", "muon p_{T} (GeV)", "abs. uncertainty on SF product",
                                f"absUnc_{n}", plotLabel="ForceTitle", outdir=outdir+"absoluteUncertainty/",
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotError=True)
            # plot relative error
            drawCorrelationPlot(prodHists[era][n], "muon #eta", "muon p_{T} (GeV)", "rel. uncertainty on SF product",
                                f"relUnc_{n}", plotLabel="ForceTitle", outdir=outdir+"relativeUncertainty/",
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotRelativeError=True)

    if args.makePreOverPost:

        outdir = outdirOriginal + productSubfolder + prePostSubfolder

        for dataType in list(prodHistsPrePost.keys()):

            for n in list(prodHistsPrePost[dataType].keys()):
                drawCorrelationPlot(prodHistsPrePost[dataType][n],
                                    "muon #eta", "muon p_{T} (GeV)", "pre/post scale factor product",
                                    f"{n}", plotLabel="ForceTitle", outdir=outdir,
                                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                    draw_both0_noLog1_onlyLog2=1, passCanvas=canvas)
                # plot absolute error
                drawCorrelationPlot(prodHistsPrePost[dataType][n],
                                    "muon #eta", "muon p_{T} (GeV)", "abs. uncertainty on SF product",
                                    f"absUnc_{n}", plotLabel="ForceTitle", outdir=outdir+"absoluteUncertainty/",
                                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                    draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotError=True)
                # plot relative error
                drawCorrelationPlot(prodHistsPrePost[dataType][n],
                                    "muon #eta", "muon p_{T} (GeV)", "rel. uncertainty on SF product",
                                    f"relUnc_{n}", plotLabel="ForceTitle", outdir=outdir+"relativeUncertainty/",
                                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                    draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotRelativeError=True)
