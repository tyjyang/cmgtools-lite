#!/usr/bin/env python

import os
import ROOT

from CMGTools.WMass.plotter.fakeRate import compileMacro

#if "/RooDoubleCBFast_cc.so" not in ROOT.gSystem.GetLibraries():
#   compileMacro("src/CMGTools/WMass/python/plotter/doubleCBfunctions/src/RooDoubleCBFast.cc")

ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/plotter/testDoubleCB/my_double_CB.cc+" % os.environ['CMSSW_BASE'])

