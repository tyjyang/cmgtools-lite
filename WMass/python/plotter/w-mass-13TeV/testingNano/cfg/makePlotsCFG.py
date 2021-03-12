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

def getQcdScaleWgt():
    qcdScaleWgt = {}
    for i,varPair in enumerate([("05","05"), ("05","1"), ("05", "2"), ("1", "05"), \
                            ("1","1"), ("1","2"), ("2", "05"), ("2", "1"), ("2", "2")]):
        qcdScaleWgt["scaleWeightMuR%sMuF%s" % varPair] = "LHEScaleWeight[%d]" % (i*2)
    return qcdScaleWgt
###################################
# SOME BASIC CONFIGS #
######################
baseHistName = "muon_eta_pt_charge"


commonPart = "Muon_charge[goodMuonsCharge][0]\:Muon_pt[goodMuonsCharge][0]\:Muon_eta[goodMuonsCharge][0]: 48,-2.4,2.4,29,26,55,2,-2,2; XTitle='Muon #eta', YTitle='Muon p_{T} (GeV)', ZTitle='Muon charge'"
####################################


#nominal
print("%s: %s \n" % (baseHistName,commonPart))

# pdf
for i in range(1,103):
    systname = "pdf{ipdf}".format(ipdf=i) if i <= 100 else "alphaSUp" if i == 101 else "alphaSDown"
    print("{bn}_{sys}: {main}, AddWeight='LHEPdfWeight[{ipdf}]', ProcessRegexp='W.*|Z.*'\n".format(bn=baseHistName,sys=systname,ipdf=i,main=commonPart))

# QCD scales
qcdScaleWgt = getQcdScaleWgt()
qcdScaleToSimplerName = {
    "muRUp"    : qcdScaleWgt["scaleWeightMuR2MuF1"],
    "muRDown"    : qcdScaleWgt["scaleWeightMuR05MuF1"],
    "muFUp"    : qcdScaleWgt["scaleWeightMuR1MuF2"],
    "muFDown"    : qcdScaleWgt["scaleWeightMuR1MuF05"],
    "muRmuFUp" : qcdScaleWgt["scaleWeightMuR2MuF2"],
    "muRmuFDown" : qcdScaleWgt["scaleWeightMuR05MuF05"]
}

NVPTBINS=10 # usually it would be 10, use less for tests, e.g. 2
for qcd in list(qcdScaleToSimplerName.keys()):
    scalevar = qcdScaleToSimplerName[qcd]
    for ipt in range(1,1+NVPTBINS):
        ptcut = wptBinsScales(ipt)
        #wgtstr = 'pow({sc}\,(Vpt_preFSR>={ptlo}&&Vpt_preFSR<{pthi}))'.format(sc=scalevar,ptlo=ptcut[0],pthi=ptcut[1])
        wgtstr = 'qcdScaleWeight_VptBinned({sc}\,Vpt_preFSR\,{ptlo}\,{pthi})'.format(sc=scalevar,ptlo=ptcut[0],pthi=ptcut[1])
        tokens = qcd.split("Up") if "Up" in qcd else qcd.split("Down")
        systname = tokens[0] + str(ipt)
        systname += ("Up" if "Up" in qcd else "Down")
        print("{bn}_{sys}: {main}, AddWeight='{wgt}', ProcessRegexp='W.*|Z.*'\n".format(bn=baseHistName,sys=systname,ipdf=1,main=commonPart,wgt=wgtstr))
    
