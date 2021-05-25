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

from copy import *
from cutsFile import *
from fakeRate import *
from mcCorrections import *

if "/functions_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/functions.cc")

if len(args) > 1:
    name = args[1]
else:
    name = 'test'
maxFiles = 500

withISR = True

# ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kViridis)
ROOT.TColor.InvertPalette()

def makePlot(name):
    print('making plot for ' + name)

    chain = ROOT.TChain('Events')

    if name == 'test':
        chain.Add('/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoAOD/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPreVFP/210319_193015/0000/NanoV8MCPreVFP_weightFix_313.root')

    else:
        path = paths[name]
        
        nFilesAdded = 0
        for f in os.listdir(path):
            if nFilesAdded > maxFiles:
                break
            if '.root' in f:
                chain.Add(path+f)
                nFilesAdded += 1

    if withISR:
        name = name + "_withISR"

    rdf = ROOT.RDataFrame(chain)
    # rdf = ROOT.RDataFrame('Events', '/afs/cern.ch/work/m/mseidel/WMass/MC/CMSSW_10_6_19_patch2/src/Configuration/WMassNanoGen/NanoGen.root')
    rdf = rdf.Define('genPrefsrLeps', 'Numba::prefsrLeptons(GenPart_status, GenPart_statusFlags, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt)')
    rdf = rdf.Define('ptVgen', 'transversemomentum(GenPart_pt[genPrefsrLeps],GenPart_phi[genPrefsrLeps])')
    # FIXME: handle pair production
    rdf = rdf.Define('ewSel', 'Numba::ewPhotonKinematicsSel(GenPart_status, GenPart_statusFlags, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, %s)' % str(withISR).lower())
    rdf = rdf.Define('sij', 'log10(invMLepPhotons(GenPart_pt[ewSel==1],GenPart_eta[ewSel==1],GenPart_phi[ewSel==1],0.105658,GenPart_pt[ewSel==3],GenPart_eta[ewSel==3],GenPart_phi[ewSel==3]))')
    rdf = rdf.Define('sik', 'log10(invMLepPhotons(GenPart_pt[ewSel==1],GenPart_eta[ewSel==1],GenPart_phi[ewSel==1],0.105658,GenPart_pt[ewSel==2],GenPart_eta[ewSel==2],GenPart_phi[ewSel==2]))')
    rdf = rdf.Define('sjk', 'log10(invMLepPhotons(GenPart_pt[ewSel==2],GenPart_eta[ewSel==2],GenPart_phi[ewSel==2],0.105658,GenPart_pt[ewSel==3],GenPart_eta[ewSel==3],GenPart_phi[ewSel==3]))')
    # rdf = rdf.Define('sijk', '(invMLepPhotons(GenPart_pt[ewSel>=2],GenPart_eta[ewSel>=2],GenPart_phi[ewSel>=2],0.105658,GenPart_pt[ewSel==3],GenPart_eta[ewSel==3],GenPart_phi[ewSel==3]))')

    # rdf.Snapshot('Events', 'snapshot_%s.root' % name, ['ptVgen', 'ewSel', 'nGenPart', 'GenPart_pdgId', 'GenPart_status', 'GenPart_statusFlags', 'GenPart_genPartIdxMother', 'GenPart_pt', 'GenPart_eta', 'GenPart_phi', 'sij', 'sik', 'sjk'])

    hist = rdf.Histo2D((name, ';log_{10} s(#mu- #sum #gamma);log_{10} s(#mu+ #sum #gamma)', 100, -1.5, 3.5, 100, -1.5, 3.5), 'sij', 'sjk')

    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.12)
    c.SetLeftMargin(0.15)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    c.SetLogz()
    c.cd()

    hist.Draw('colz')

    c.Print('plots/ewPhotonKinematics/ewPK_%s.pdf' % name)
    c.Print('plots/ewPhotonKinematics/ewPK_%s.png' % name)
    c.Print('plots/ewPhotonKinematics/ewPK_%s.root' % name)


paths = {
    'minnlo': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoAOD/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPreVFP/210319_193015/0000/',
    'horace-photos': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photos-isr-pythia/',
    'horace-pythia': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-pythia-isr-pythia/',
    'horace-exp': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-exp-fsr-off-isr-off/',
    'horace-exp-old': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia/',
    'horace-alpha-old': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-alpha-old-fsr-off-isr-pythia/',
    'horace-photoslow': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photoslow-isr-pythia/',
}

if name == 'all':
    for pathname in paths:
        makePlot(pathname)
else:
    makePlot(name)