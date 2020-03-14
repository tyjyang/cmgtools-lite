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


inputfile = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/plots/Wlike/TREE_4_WLIKE_MU/Zratio_OldAndnewMC//plots_zmm.root"
outdir = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/plots/Wlike/TREE_4_WLIKE_MU/Zratio_OldAndnewMC/inverseRatio/"
outfiletag = "ZptRatio_newlOverOldDY_dressed_someCuts"

createPlotDirAndCopyPhp(outdir)

outfile = outdir + outfiletag + ".root"

f = ROOT.TFile.Open(inputfile)
oldh = f.Get("genZptDressed_oldZ")
newh = f.Get("genZptDressed_newZ")
newhPhotos = f.Get("genZptDressed_newZ_photos")
oldh.SetDirectory(0)
newh.SetDirectory(0)
newhPhotos.SetDirectory(0)
f.Close()


ratio_old_new = oldh.Clone("ZptRatio_dressed_oldOverNew")
ratio_old_newPhotos = oldh.Clone("ZptRatio_dressed_oldOverNewPhotos")

ratio_old_new.Divide(newh)
ratio_old_newPhotos.Divide(newhPhotos)

hists = [ratio_old_new, ratio_old_newPhotos]
legs = ["MadGraph_aMC@NLO / powhegMiNNLO_pythia8", "MadGraph_aMC@NLO / powhegMiNNLO_pythia8_photos"]

adjustSettings_CMS_lumi()
canvas = ROOT.TCanvas("canvas","",1000,900)

drawNTH1(hists,legs,"Z p_{T} (dressed) [GeV]", "ratio", outfiletag, outdir, 
         draw_both0_noLog1_onlyLog2=1,
         #labelRatioTmp="nomi/photos::0.98,1.02",
         legendCoords="0.15,0.95,0.7,0.9",
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
