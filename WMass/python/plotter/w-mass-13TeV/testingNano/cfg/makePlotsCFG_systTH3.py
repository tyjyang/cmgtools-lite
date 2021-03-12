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

def wptBinsScales(i):
    # 5% quantiles (to be redone on the new MC)
    wptbins = [0.0, 1.971, 2.949, 3.838, 4.733, 5.674, 6.684, 7.781, 8.979, 10.303, 11.777, 13.435, 15.332, 17.525, 20.115, 23.245, 27.173, 32.414, 40.151, 53.858, 13000.0]
    if len(wptbins)<2*i:
        print('You are asking too much from the wpt binning for decorrelation of scales')
    ptlo = wptbins[2*(i-1)]
    pthi = wptbins[2*i]
    return [ptlo, pthi]

###################################
# SOME BASIC CONFIGS #
######################
baseHistName = "muon_eta_pt"

axisNames = "XTitle='Muon #eta', YTitle='Muon p_{T} (GeV)'"
etaptBins = "48,-2.4,2.4,29,26,55" 
####################################

#nominal
print("{n}: Muon_pt[goodMuonsCharge][0]\:Muon_eta[goodMuonsCharge][0]: {b}; {axis} \n".format(n=baseHistName,b=etaptBins,axis=axisNames))

# pdf
syst_key = "pdf"
syst_expr = "indices(LHEPdfWeight)\:scalarToRVec(Muon_pt[goodMuonsCharge][0],LHEPdfWeight)\:scalarToRVec(Muon_eta[goodMuonsCharge][0],LHEPdfWeight)"
syst_binning = "102,0.5,102.5"
syst_axisname = "ZTitle='PDF index'"
syst_weight  = "LHEPdfWeight"
process_regexpr = "W.*|Z.*"

print("{n}_{sk}: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(n=baseHistName, sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr))

# qcd scales (not Vpt binned, that one comes just afterward)
syst_key = "qcdScale"
syst_expr = "indices(LHEScaleWeight)\:scalarToRVec(Muon_pt[goodMuonsCharge][0],LHEScaleWeight)\:scalarToRVec(Muon_eta[goodMuonsCharge][0],LHEScaleWeight)"
syst_binning = "18,-0.5,17.5"
syst_axisname = "ZTitle='QCD scale index'"
syst_weight  = "LHEScaleWeight"
process_regexpr = "W.*|Z.*" # actually one or the other based on Wlike or Wmass

print("{n}_{sk}: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(n=baseHistName, sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr))

# qcd scales (Vpt binned)
NVTPBINS = 10
for ipt in range(1,1+NVTPBINS):
    syst_key = "qcdScaleVptBin%d" % ipt
    ptcut = wptBinsScales(ipt)
    syst_weight  = "qcdScaleWeight_VptBinned(LHEScaleWeight\,ptVgen\,{ptlo}\,{pthi})".format(ptlo=ptcut[0],pthi=ptcut[1])
    process_regexpr = "W.*|Z.*" # actually one or the other based on Wlike or Wmass

    print("{n}_{sk}: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(n=baseHistName, sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr))

