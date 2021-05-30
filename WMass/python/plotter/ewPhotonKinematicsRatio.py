#!/usr/bin/env python
from math import *
import re
import os, os.path
from array import array
import logging
from collections import OrderedDict
import numbaDefines
from math import sqrt

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

withISR = True

ROOT.gStyle.SetOptStat(0)

def makeRatio(name1, name2):
    print('making plot for %s/%s' % (name1, name2))

    hist = hists[name1].Clone()
    hist.Divide(hists[name2])
    # hist.Smooth(1)
    
    for i in range(hist.GetNbinsX()+2):
        for j in range(hist.GetNbinsY()+2):
            val = hist.GetBinContent(i, j)
            if val == 0.:
                hist.SetBinContent(i, j, 1.)
                continue
            errWeight = 2.
            relErr = min(1., errWeight * hist.GetBinError(i, j) / val)
            ipolVal = (1.-relErr) * val + relErr
            hist.SetBinContent(i, j, ipolVal)

    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.12)
    c.SetLeftMargin(0.15)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    c.cd()

    ROOT.gStyle.SetPalette(ROOT.kRedBlue)
    ROOT.gStyle.SetNumberContours(100)
    ROOT.TColor.InvertPalette()

    hist.GetZaxis().SetRangeUser(0,3)
    hist.SetContour(100)
    hist.Draw('colz')

    c.Print('plots/ewPhotonKinematics/ewPKRatio_%s_%s.pdf'  % (name1, name2))
    c.Print('plots/ewPhotonKinematics/ewPKRatio_%s_%s.png'  % (name1, name2))
    c.Print('plots/ewPhotonKinematics/ewPKRatio_%s_%s.root' % (name1, name2))

def getHist(name):
    tfile = ROOT.TFile('plots/ewPhotonKinematics/ewPK_%s.root' % name)
    canvas = tfile.Get('c')
    hist = canvas.GetPrimitive(name)

    # merge no/low-emission bins
    maxBinLE = hist.GetXaxis().FindFixBin(0.)
    
    # integralLE = hist.Integral(11, maxBinLE, 11, maxBinLE)
    integralLE = 0.
    integralLEerr = 0.
    edgeRadius = 9
    for i in range(11, maxBinLE):
        for j in range(11, maxBinLE):
            dx = max(0., i-maxBinLE+edgeRadius)
            dy = max(0., j-maxBinLE+edgeRadius)
            if sqrt(dx*dx + dy*dy) < edgeRadius-0.5:
                integralLE += hist.GetBinContent(i, j)
                integralLEerr += hist.GetBinError(i, j)
    for i in range(11, maxBinLE):
        for j in range(11, maxBinLE):
            dx = max(0., i-maxBinLE+edgeRadius)
            dy = max(0., j-maxBinLE+edgeRadius)
            if sqrt(dx*dx + dy*dy) < edgeRadius-0.5:
                hist.SetBinContent(i, j, integralLE)
                hist.SetBinError(i, j, integralLEerr)

    hist.Scale(1./hist.Integral())
    return hist

paths = {
    'minnlo': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoAOD/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPreVFP/210319_193015/0000/',
    'horace-photos': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photos-isr-pythia/',
    'horace-pythia': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-pythia-isr-pythia/',
    'horace-exp': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-exp-fsr-off-isr-off/',
    'horace-exp-old': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia/',
    'horace-alpha-old': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-alpha-old-fsr-off-isr-pythia/',
    'horace-photoslow': '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photoslow-isr-pythia/',
}

name1 = args[1]
name2 = args[2]

hists = {}
# for name in paths:
for name in [name1, name2]:
    hists[name] = getHist(name + '_withISR')

makeRatio(name1, name2)