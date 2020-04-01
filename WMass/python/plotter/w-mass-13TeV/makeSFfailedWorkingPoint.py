#!/bin/env python

# example
# python w-mass-13TeV/makeSFfailedWorkingPoint.py ../postprocessing/data/leptonSF/new2016_madeSummer2018/smoothEfficiency_muons_recoToSel_finerETA.root -o plots/scaleFactors_Final/muon/recoToSelection_pt_25_55_eta0p1_finalJuly2019/SF_failWorkingPoint/ --palette 57 -n scaleFactorSmooth_failIDandIso -t "fail ID+isolation"

import ROOT, os, datetime, re, operator, math
from array import array

ROOT.gROOT.SetBatch(True)

import utilities
utilities = utilities.util()

import sys
sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

def makeOneMinusEffHistogram(hnew,hold):
    for ix in range(1,1+hold.GetNbinsX()):
        for iy in range(1,1+hold.GetNbinsY()):
            hnew.SetBinContent(ix,iy,1-hold.GetBinContent(ix,iy))
            hnew.SetBinError(ix,iy,1-hold.GetBinError(ix,iy))

if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)
    from optparse import OptionParser
    parser = OptionParser(usage='%prog th2.root [options] ')
    parser.add_option('-o','--outdir',     dest='outdir',     default='',   type='string', help='outdput directory to save plots')
    parser.add_option(     '--palette'  , dest='palette',      default=55, type=int, help='Set palette: use a negative number to select a built-in one, otherwise the default is 55 (kRainbow)')
    parser.add_option('-n','--outname',     dest='outname',     default='scaleFactorFail',   type='string', help='name of output plot')
    parser.add_option('-t','--title',     dest='title',     default='fail working point',   type='string', help='title of histogram (written in canvas)')
    parser.add_option('-r','--zrange',  dest='zrange',  default=(0, 0),type="float", nargs=2, help="Min and max for the z axis of the plot")
    (options, args) = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    if len(args) < 1:
        parser.print_usage()
        quit()

    outdir = options.outdir
    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    tf = ROOT.TFile.Open(args[0])            
    effData = tf.Get("hdataSmoothCheck")
    effMC = tf.Get("hmcSmoothCheck")
    effData.SetDirectory(0)
    effMC.SetDirectory(0)
    tf.Close()

    effDataFail = effData.Clone("effDataFail")
    effMCFail   = effMC.Clone("effMCFail")
    makeOneMinusEffHistogram(effDataFail,effData)
    makeOneMinusEffHistogram(effMCFail,effMC)
    sfFail = effDataFail.Clone("sfFail")
    sfFail.Divide(effMCFail)
    zaxistitle = "scale factor"
    zmin = options.zrange[0]
    zmax = options.zrange[1]
    if zmin < zmax:
        zaxistitle = zaxistitle + "::" + str(zmin) + "," + str(zmax)
    sfFail.GetZaxis().SetTitle(zaxistitle)
    print "Min SF = %.3f" % sfFail.GetBinContent(sfFail.GetMinimumBin())
    print "Max SF = %.3f" % sfFail.GetBinContent(sfFail.GetMaximumBin())
    
    if options.title:
        sfFail.SetTitle(options.title)

    adjustSettings_CMS_lumi()
    canvas2D = ROOT.TCanvas("canvas2D","",900,900)
    drawCorrelationPlot(sfFail,sfFail.GetXaxis().GetTitle(),sfFail.GetYaxis().GetTitle(),sfFail.GetZaxis().GetTitle(),
                     options.outname,"ForceTitle",outdir,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D)

