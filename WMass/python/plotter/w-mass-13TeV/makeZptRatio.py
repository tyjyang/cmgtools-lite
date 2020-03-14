#!/bin/env python

import array
import json
import ROOT
import re 
import os
import sys

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gROOT.SetBatch(True)

use_binApproxSMP17010 = False
inputfile = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/plots/Wlike/TREE_4_WLIKE_MU/Zratio_OldAndnewMC//plots_zmm.root"
outdir = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/plots/Wlike/TREE_4_WLIKE_MU/Zratio_OldAndnewMC/inverseRatio_MGaMCAtNLO_PowhegMiNNLO" # no ending "/", it is added below
outfiletag = "ZptRatio_MGaMCAtNLO_PowhegMiNNLO_dressed" 
var = "genZptDressed" 

if use_binApproxSMP17010: 
    var += "_binApproxSMP17010"
    outdir = outdir + "_binApproxSMP17010"

rescaleBy_withCut_over_noCut = True
rescaleFile = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/plots/Wlike/TREE_4_WLIKE_MU/testZ_newMC_checkCuts_onGenZpt/plots_zmm.root"
# these have 2 GeV binning from 0 to 100
rescaleHistoNum = var + "_newZ_someRecoCuts"
rescaleHistoDen = var + "_newZ"
rescaleHistoPhotosNum = var + "_newZ_photos_someRecoCuts"
rescaleHistoPhotosDen = var + "_newZ_photos"
if rescaleBy_withCut_over_noCut:
    outfiletag += "_noCuts"
else:
    outfiletag += "_someCuts"


outdir += "/"
createPlotDirAndCopyPhp(outdir)
outfile = outdir + outfiletag + ".root"

        
f = ROOT.TFile.Open(inputfile)
oldh = f.Get(var + "_oldZ")
newh = f.Get(var + "_newZ")
newhPhotos = f.Get(var + "_newZ_photos")
oldh.SetDirectory(0)
newh.SetDirectory(0)
newhPhotos.SetDirectory(0)
f.Close()

ratio_old_new = oldh.Clone("pythia8")
ratio_old_newPhotos = oldh.Clone("pythia8_photos")

ratio_old_new.Divide(newh)
ratio_old_newPhotos.Divide(newhPhotos)

hrescale = None
hrescalePhotos = None
if rescaleBy_withCut_over_noCut:
    f = ROOT.TFile.Open(rescaleFile)
    if not f:
        print "Error opening file %s" % f
        quit()
    num = f.Get(rescaleHistoNum)
    den = f.Get(rescaleHistoDen)
    num.SetDirectory(0)
    den.SetDirectory(0)
    numphotos = f.Get(rescaleHistoPhotosNum)
    denphotos = f.Get(rescaleHistoPhotosDen)
    numphotos.SetDirectory(0)
    denphotos.SetDirectory(0)
    f.Close()

    hrescale = num.Clone("hrescale_withCuts_over_noCuts")
    hrescale.Scale(den.Integral()/hrescale.Integral())
    hrescale.Divide(den)
    ratio_old_new.Multiply(hrescale)

    hrescalePhotos = numphotos.Clone("hrescalePhotos_withCuts_over_noCuts")
    hrescalePhotos.Scale(denphotos.Integral()/hrescalePhotos.Integral())
    hrescalePhotos.Divide(denphotos)
    ratio_old_newPhotos.Multiply(hrescalePhotos)

    histsRescale = [hrescale, hrescalePhotos]
    legsRescale = ["pythia8", "pythia8_photos"]    

    
hists = [ratio_old_new, ratio_old_newPhotos]
legs = ["MadGraph5_aMC@NLO / powhegMiNNLO_pythia8", "MadGraph5_aMC@NLO / powhegMiNNLO_pythia8_photos"]
ytitle = "ratio"
#legs = ["PYTHIA8", "PYTHIA8 + PHOTOS"]
#ytitle = "MadGraph5_aMC@NLO / POWHEG-MiNNLO"

adjustSettings_CMS_lumi()
canvas = ROOT.TCanvas("canvas","",1000,900)

drawNTH1(hists,legs,"Z p_{T} (dressed) [GeV]", ytitle, outfiletag, outdir, 
         draw_both0_noLog1_onlyLog2=1,
         #labelRatioTmp="nomi/photos::0.98,1.02",
         legendCoords="0.15,0.97,0.74,0.9",
         lowerPanelHeight=0.0,
         passCanvas=canvas)

if rescaleBy_withCut_over_noCut:
    drawNTH1(histsRescale,legsRescale,"Z p_{T} (dressed) [GeV]", "scaling: with/no cuts", "scalingFactor_PowhegMiNNLO_cutsOverNoCuts", outdir, 
             draw_both0_noLog1_onlyLog2=1,
             #labelRatioTmp="nomi/photos::0.98,1.02",
             legendCoords="0.15,0.97,0.74,0.9",
             lowerPanelHeight=0.0,
             passCanvas=canvas)


of = ROOT.TFile.Open(outfile,"recreate")
of.cd()
for h in hists:
    h.Write()
of.Close()
print "-"*30
print "I wrote histograms in " + outfile
print "-"*30
print ""
