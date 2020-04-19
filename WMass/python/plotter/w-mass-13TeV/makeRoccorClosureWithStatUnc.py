#!/bin/env python

import ROOT, os, sys, re, array, math, copy
import time

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gROOT.SetBatch(True)

indir = "/afs/cern.ch/user/m/mciprian/www/wmass/13TeV/Wlike/test_zFitter/"
fplus = "newRoccor_fullSel_chargePlus_allEvents_trigMatchPlus_signEta_testFit_2CB_altPtBins/"
fminus = "newRoccor_fullSel_chargeMinus_allEvents_trigMatchMinus_signEta_testFit_2CB_altPtBins/"
fname = "plot_dm_diff.root"
hname = "plot_dm_diff"

def setStyleTH2(h2,ch):

    if ch == "plus": 
        charge = "positive"
    elif ch == "minus": 
        charge = "negative"

    h2.SetTitle("")
    h2.SetMarkerSize(1.2)
    # X
    h2.GetXaxis().SetTitle("%s muon #eta" % charge)
    h2.GetXaxis().SetTitleSize(0.05)
    h2.GetXaxis().SetLabelSize(0.04)
    h2.GetXaxis().SetTitleOffset(0.95)
    # Y
    h2.GetYaxis().SetTitle("%s muon p_{T} (GeV)" % charge)
    h2.GetYaxis().SetTitleSize(0.05)
    h2.GetYaxis().SetLabelSize(0.04)
    h2.GetYaxis().SetTitleOffset(0.95)
    # Z
    h2.GetZaxis().SetTitleSize(0.05)
    h2.GetZaxis().SetLabelSize(0.04)
    h2.GetZaxis().SetTitleOffset(1.2) # 1.4  

if __name__ == "__main__":
            
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    (options, args) = parser.parse_args()

    ROOT.TH1.StatOverflows(True) # to use under/overflow bins when computing rms of histogram
    ROOT.gStyle.SetNumberContours(51) # default is 20 
    ROOT.gStyle.SetPaintTextFormat(".3f")
    ROOT.gStyle.SetOptStat(0)
    ROOT.TColor.CreateGradientColorTable(3,
                                         array ("d", [0.00, 0.50, 1.00]),
                                         array ("d", [0.00, 1.00, 1.00]),
                                         array ("d", [0.34, 1.00, 0.65]),
                                         array ("d", [0.82, 1.00, 0.00]),
                                         255,  0.95)



    hrms = ROOT.TH1D("htmp","",200,-0.1,0.1)

    for charge in ["plus", "minus"]:

        inputDirNomi = indir + (fplus if charge == "plus" else fminus)
        fnomi = inputDirNomi + fname
        fn = ROOT.TFile.Open(fnomi)
        hnomi = fn.Get(hname)
        hnomi.SetDirectory(0)
        fn.Close()
        hTotStat = copy.deepcopy(hnomi.Clone("hTotStat_{c}".format(c=charge)))
        hTotStat.Reset()
        hFitUnc = copy.deepcopy(hnomi.Clone("hFitUnc_{c}".format(c=charge))) # to draw fit error only on TH2
        hFitUnc.Reset()
        hstat = []
        for i in range(100):
            fstat = indir + "newRoccor_charge{c}_stat{d}/".format(c="Plus" if charge == "plus" else "Minus",d=i) + fname
            fn = ROOT.TFile.Open(fstat)
            if not fn:
                print "Warning: could not open file %s" % fstat
                quit()
            hstat.append(fn.Get(hname))
            hstat[i].SetDirectory(0)
            hstat[i].Add(hnomi,-1.0)
            fn.Close()
        for ix in range(1,1+hnomi.GetNbinsX()):
            for iy in range(1,1+hnomi.GetNbinsY()):
                hrms.Reset()
                for i in range(100):
                    hrms.Fill(hstat[i].GetBinContent(ix,iy))
                hTotStat.SetBinContent(ix,iy,hnomi.GetBinContent(ix,iy))
                totErr = hrms.GetStdDev() * hrms.GetStdDev() + hnomi.GetBinError(ix,iy) * hnomi.GetBinError(ix,iy)
                hTotStat.SetBinError(ix,iy,ROOT.TMath.Sqrt(totErr))
                #hFitUnc.SetBinContent(ix,iy,hnomi.GetBinError(ix,iy))
                hFitUnc.SetBinContent(ix,iy,hrms.GetStdDev())


        # draw histograms
        #adjustSettings_CMS_lumi()
        canvas = ROOT.TCanvas("canvas","",1800,1200)
        canvas.SetRightMargin(0.2)
        canvas.SetTickx(1)
        canvas.SetTicky(1)
        canvas.cd()
        setStyleTH2(hTotStat,charge)
        hTotStat.SetBarOffset(0.2)
        hTotStat.Draw("COLZ TEXTE")
        hFitUnc.SetMarkerColor(ROOT.kRed+2)
        hFitUnc.SetBarOffset(-0.2)
        hFitUnc.Draw("TEXT SAME")
        for ext in ['png', 'pdf']:
            canvas.SaveAs('{od}/{cn}_withStatUnc.{ext}'.format(od=inputDirNomi, cn=hname, ext=ext))
