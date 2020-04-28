#!/bin/env python

import ROOT, os, sys, re, array, math, copy
import time

from makeRoccorClosureWithStatUnc import setStyleTH2

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gROOT.SetBatch(True)

indir = "/afs/cern.ch/user/m/mciprian/www/wmass/13TeV/Wlike/test_zFitter/"
dirnomi = "newRoccor_fullSel_chargePlus_allEvents_trigMatchPlus_signEta_testFit_2CB_altPtBins/"
charge = "plus"
fname = "plot_dm_diff.root"
hname = "plot_dm_diff"
hvar = {"noSFandPU"    : "no scale factors and PU weight",
        "4xMassBins"   : "4x mass granularity",
        "upPtCut60"    : "p_{T} < 60 GeV"
}

if __name__ == "__main__":
            
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    (options, args) = parser.parse_args()

    #ROOT.TH1.StatOverflows(True) # to use under/overflow bins when computing rms of histogram
    ROOT.gStyle.SetNumberContours(51) # default is 20 
    ROOT.gStyle.SetPaintTextFormat(".3f")
    #ROOT.gStyle.SetOptStat(0) # uncomment to disable stat box, but in some cases it is needed
    ROOT.TColor.CreateGradientColorTable(3,
                                         array ("d", [0.00, 0.50, 1.00]),
                                         array ("d", [0.00, 1.00, 1.00]),
                                         array ("d", [0.34, 1.00, 0.65]),
                                         array ("d", [0.82, 1.00, 0.00]),
                                         255,  0.95)


    inputDirNomi = indir + dirnomi
    outdir = inputDirNomi + "comparisons/"
    createPlotDirAndCopyPhp(outdir)
    fnomi = inputDirNomi + fname
    fn = ROOT.TFile.Open(fnomi)
    hnomi = fn.Get(hname)
    hnomi.SetDirectory(0)
    fn.Close()

    adjustSettings_CMS_lumi()

    canvas = ROOT.TCanvas("canvas","",1800,1200)
    canvas.SetRightMargin(0.2)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    
    canvas1D = ROOT.TCanvas("canvas1D","",1000,1000)
    canvas.SetBottomMargin(0.14)

    h1D = ROOT.TH1D("h1D","",40,-5,5)
    histvars = {}
    for k in hvar.keys():
        canvas.cd()
        fvar = indir + dirnomi.rstrip('/') + "_" + k + "/" + fname
        fn = ROOT.TFile.Open(fvar)
        if not fn:
            print "Warning: could not open file %s" % fvar
            quit()
        histvars[k] = fn.Get(hname)
        if not histvars[k]:
            print "Warning: could not get histogram %s from file %s" % (hname,fvar)
            quit()
        histvars[k].SetDirectory(0)
        histvars[k].Add(hnomi,-1.0)

        # draw histograms
        histvars[k].SetStats(0)
        setStyleTH2(histvars[k],charge)
        ROOT.gStyle.SetPaintTextFormat(".3f")
        histvars[k].GetZaxis().SetRangeUser(-0.05,0.05)
        histvars[k].GetZaxis().SetTitle("difference (GeV): nomi - {v}".format(v=k))
        histvars[k].Draw("COLZ TEXTE")
        for ext in ['png', 'pdf']:
            canvas.SaveAs('{od}/{cn}_compare2D_{k}.{ext}'.format(od=outdir, cn=hname, k=k, ext=ext))

        h1D.Reset()
        for ix in range(1,1+histvars[k].GetNbinsX()):
            for iy in range(1,1+histvars[k].GetNbinsY()):
                h1D.Fill(histvars[k].GetBinContent(ix,iy)/histvars[k].GetBinError(ix,iy))

        h1D.SetTitle(hvar[k])
        canvas1D.SetTitle(hvar[k])
        text = "'nominal' vs '{x}'".format(x=hvar[k])
        drawTH1(h1D,"pull of closure","Entries",outdir,
                prefix="{cn}_compare1D_{k}".format(cn=hname, k=k),
                passCanvas=canvas1D,
                moreTextLatex=(text+"::0.15,0.95,0.08,0.045"))
