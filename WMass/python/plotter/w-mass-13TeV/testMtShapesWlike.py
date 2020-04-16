#!/bin/env python

# example
# python w-mass-13TeV/testMtShapesWlike.py -o plots/Wlike/TREE_4_WLIKE_MU/test_QCDbks_ss_onlyPromptMC/noDxy/compareShapes/compareMtMzInRegions_rebinMtMz_2_2/ --rebinMtMz 2 2 --palette 57

import ROOT, os, datetime, re, operator, math
from array import array

ROOT.gROOT.SetBatch(True)

import utilities
utilities = utilities.util()

import sys
sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

inputfolder = "plots/Wlike/TREE_4_WLIKE_MU/test_QCDbks_ss_onlyPromptMC/noDxy/"

files_ss = {"ss_passIso" : "MtvsMz_sameSignRegion_AllEvents_passIso/plus/",
            "ss_failIso" : "MtvsMz_sameSignRegion_AllEvents_failIso/plus/",
        }

files = files_ss.copy()
#files.update(files_os)

fname = "plots_zmm.root"
hname = "mt_wlike__zmass"


if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option('-o','--outdir',     dest='outdir',     default='',   type='string', help='outdput directory to save plots')
    parser.add_option(     '--palette'  , dest='palette',      default=55, type=int, help='Set palette: use a negative number to select a built-in one, otherwise the default is 55 (kRainbow)')
    parser.add_option(     '--rebinMtMz'  , dest='rebinMtMz',      default=(0,0), nargs=2, type=int, help='Rebinnign factor for mt-mz distribution. Default is none, equivalent to 1,1')
    (options, args) = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    outdir = options.outdir
    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    rebinMt = options.rebinMtMz[0]
    rebinMz = options.rebinMtMz[1]

    hists2 = {}
    hists1mz = {}
    hists1mt = {}

    adjustSettings_CMS_lumi()
    canvas2D = ROOT.TCanvas("canvas2D","",900,900)

    for f in files:
        fileFullName = inputfolder + files[f] + fname
        tf = ROOT.TFile.Open(fileFullName)
        if not tf:
            print "Error when opening file %s" % fileFullName
            quit()
        hname_full = hname
        hdata = tf.Get(hname_full + "_data")
        hbkg = tf.Get(hname_full + "_background")
        if not hdata:
            print "Error fetching hdata %s from file %s" % (hname_full,fileFullName)
            quit()
        hists2[f] = ROOT.TH2D("MtMz_{f}".format(f=f),"",
                              hdata.GetNbinsX(),hdata.GetXaxis().GetBinLowEdge(1),hdata.GetXaxis().GetBinLowEdge(1+hdata.GetNbinsX()),
                              hdata.GetNbinsY(),hdata.GetYaxis().GetBinLowEdge(1),hdata.GetYaxis().GetBinLowEdge(1+hdata.GetNbinsY()))
        hists2[f].SetDirectory(0)
        hists2[f].Add(hdata,hbkg,1.0,-1.0) 
        hists2[f].SetTitle(f)
        hists2[f].RebinX(rebinMt)
        hists2[f].RebinY(rebinMz)
        hists1mz[f] = hists2[f].ProjectionX("Mz_{f}".format(f=f),0,-1,"e")
        hists1mz[f].SetDirectory(0)
        hists1mt[f]  = hists2[f].ProjectionY("Mt_{f}".format(f=f),0,-1,"e")
        hists1mt[f].SetDirectory(0)
        tf.Close()
        drawCorrelationPlot(hists2[f],"dimuon invariant mass [GeV]","W-like transverse mass [GeV]","Events",
                            hists2[f].GetName(),"ForceTitle",outdir,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D)


    canvas1D = ROOT.TCanvas("canvas1D","",900,900)    

    h1sort_mt_ss = []
    h1sort_mz_ss = []
    legEntries_ss = []
    for f in sorted(files_ss.keys()):
        h1sort_mz_ss.append(hists1mz[f].Clone("clone_"+hists1mz[f].GetName()))
        h1sort_mt_ss.append(hists1mt[f].Clone("clone_"+hists1mt[f].GetName()))
        legEntries_ss.append(f)

    drawNTH1(h1sort_mt_ss,legEntries_ss,"W-like transverse mass [GeV]","Events","compareFakesVsMt_ss",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_mz_ss,legEntries_ss,"dimuon invariant mass [GeV]","Events","compareFakesVsMz_ss",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

    ## normalized plots

    for i,h in enumerate(h1sort_mz_ss):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_mt_ss):
        h.Scale(1./h.Integral())

    drawNTH1(h1sort_mt_ss,legEntries_ss,"W-like transverse mass [GeV]","arbitrary units","compareFakesVsPt_ss_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_mz_ss,legEntries_ss,"dimuon invariant mass [GeV]","arbitrary units::0.0,0.22","compareFakesVsEta_ss_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)
