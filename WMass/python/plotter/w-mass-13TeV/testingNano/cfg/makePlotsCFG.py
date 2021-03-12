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

#line = "muon_pt_eta_pdf: Muon_pt[goodMuonsCharge][0]\:Muon_eta[goodMuonsCharge][0]: 48,-2.4,2.4,29,26,55; XTitle='muon #eta', YTitle='muon p_{T} (GeV)', AddWeight=\"LHEPdfWeight[X]\""

baseHistName = "muon_eta_pt_charge"

commonPart = "Muon_charge[goodMuonsCharge][0]\:Muon_pt[goodMuonsCharge][0]\:Muon_eta[goodMuonsCharge][0]: 48,-2.4,2.4,29,26,55,2,-2,2; XTitle='Muon #eta', YTitle='Muon p_{T} (GeV)', ZTitle='Muon charge'"

#nominal
print("%s: %s \n" % (baseHistName,commonPart))

# pdf
for i in range(1,103):
    print("%s_pdf%d: %s, AddWeight='LHEPdfWeight[%d]', ProcessRegexp='W.*|Z.*'\n" % (baseHistName,i,commonPart,i))
