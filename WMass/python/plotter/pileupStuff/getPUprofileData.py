#!/usr/bin/env python3

# run brilcalc to get pu profiles per era (for the inclusive one you can just run 1 command and be fine)
# it also plots the profiles
# requires cmssw-cc7 singularity environment to run, and cmsenv from a release, then can be used from anywhere
# should call with python3

import os, re, array, math
import argparse
## safe batch mode
import sys

sys.path.append(os.getcwd() + "/lumiStuff/")
from runPerEra import lumiForEra_UL_new as lfe

#pileupFile = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt"
pileupFile = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest.txt"
eras = list(lfe.keys())

eras.extend(["preVFP", "postVFP"])

for era in eras:
    cmd = f"pileupCalc.py -i lumiStuff/NanoAOD_json_Run2016{era}_GoldenJsonFilter.txt --inputLumiJSON {pileupFile} --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 pileupStuff/pileupProfileData_2016Legacy_Run{era}_04June2021.root"
    print(cmd)
