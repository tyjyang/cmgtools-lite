#!/usr/bin/env python
from math import *
import re
import os, os.path
from array import array
import logging
from collections import OrderedDict
import numbaDefines

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gInterpreter.ProcessLine(".O3")
ROOT.EnableImplicitMT()

if "/functions_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/functions.cc")

rdf = ROOT.RDataFrame('Events', '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoAOD/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPreVFP/210319_193015/0000/NanoV8MCPreVFP_weightFix_313.root')
rdf.Define('genPrefsrLeps', 'Numba::prefsrLeptons(GenPart_status, GenPart_statusFlags, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt)')
rdf.Define('ptVgen', 'transversemomentum(GenPart_pt[genPrefsrLeps],GenPart_phi[genPrefsrLeps])')

rdf.Snapshot('Events', 'snapshot.root')