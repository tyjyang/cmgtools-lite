#!/bin/env python

import ROOT, os, sys, re, array, math, copy
import time

## safe batch mode 
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

indir = "/afs/cern.ch/user/m/mciprian/www/wmass/13TeV/Wlike/test_zFitter/"
#fplus = "newRoccor_fullSel_chargePlus_allEvents_trigMatchPlus_signEta_testFit_2CB_altPtBins/"
#fminus = "newRoccor_fullSel_chargeMinus_allEvents_trigMatchMinus_signEta_testFit_2CB_altPtBins/"
#fplus = "newRoccor_fullSel_chargePlus_allEvents_trigMatchPlus_signEta_testFit_2CB_1ptBin23to60_0p4eta/"
#fminus = "newRoccor_fullSel_chargeMinus_allEvents_trigMatchMinus_signEta_testFit_2CB_1ptBin23to60_0p4eta/"
fdir = "newRoccorValidation_nominal_templateFitScaleNoSmear_pt14_eta0p4"
#fname = "plot_dm_diff.root"
#hname = "plot_dm_diff"
fname = "plot_dm.root"
hname = "plot_dm"

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
    parser.add_option("-c", "--charge", dest="charge",   type="string", default="plus,minus", help="Charges to run closure on");
    (options, args) = parser.parse_args()

    ROOT.TH1.StatOverflows(True) # to use under/overflow bins when computing rms of histogram
    ROOT.gStyle.SetNumberContours(51) # default is 20 
    ROOT.gStyle.SetPaintTextFormat(".5f")
    ROOT.gStyle.SetOptStat(0)
    ROOT.TColor.CreateGradientColorTable(3,
                                         array ("d", [0.00, 0.50, 1.00]),
                                         array ("d", [0.00, 1.00, 1.00]),
                                         array ("d", [0.34, 1.00, 0.65]),
                                         array ("d", [0.82, 1.00, 0.00]),
                                         255,  0.95)



    hrms = ROOT.TH1D("htmp","",200,-0.1,0.1)

    charges = [x for x in options.charge.split(',')]
    for charge in charges:

        inputDirNomi = indir + fdir + "/" + charge + "/"
        fnomi = inputDirNomi + fname
        fn = ROOT.TFile.Open(fnomi)
        hnomi = fn.Get(hname)
        hnomi.SetDirectory(0)
        fn.Close()
        hTotUnc = copy.deepcopy(hnomi.Clone("hTotUnc_{c}".format(c=charge)))
        hTotUnc.Reset()
        hStatUncOnly = copy.deepcopy(hnomi.Clone("hStatUncOnly_{c}".format(c=charge))) # to draw fit error only on TH2
        hStatUncOnly.Reset()
        hstat = []

        nbinsClosure = hnomi.GetNbinsX()*hnomi.GetNbinsY()
        gdiff = ROOT.TGraphErrors(nbinsClosure)
        gdiff.SetName("graph_{c}".format(c=charge))
        gdiffStatUncOnly = ROOT.TGraphErrors(hnomi.GetNbinsX()*hnomi.GetNbinsY())
        gdiffStatUncOnly.SetName("graph_statOnly_{c}".format(c=charge))

        for i in range(100):
            fstat = inputDirNomi + "RoccorStatVar/charge{c}_stat{d}/{ch}/".format(c="Plus" if charge == "plus" else "Minus",d=i,ch=charge) + fname
            fn = ROOT.TFile.Open(fstat)
            if not fn:
                print "Warning: could not open file %s" % fstat
                quit()
            hstat.append(fn.Get(hname))
            hstat[i].SetDirectory(0)
            hstat[i].Add(hnomi,-1.0)
            fn.Close()
        ig = 0
        for iy in range(1,1+hnomi.GetNbinsY()):
            for ix in range(1,1+hnomi.GetNbinsX()):
                hrms.Reset()
                for i in range(100):
                    hrms.Fill(hstat[i].GetBinContent(ix,iy))
                hTotUnc.SetBinContent(ix,iy,hnomi.GetBinContent(ix,iy))
                totErr = hrms.GetStdDev() * hrms.GetStdDev() + hnomi.GetBinError(ix,iy) * hnomi.GetBinError(ix,iy)
                hTotUnc.SetBinError(ix,iy,ROOT.TMath.Sqrt(totErr))
                #hFitUnc.SetBinContent(ix,iy,hnomi.GetBinError(ix,iy))
                hStatUncOnly.SetBinContent(ix,iy,hrms.GetStdDev())
                
                # prepare graph
                gdiff.SetPoint(ig,      hTotUnc.GetBinContent(ix,iy), ig+0.5)
                gdiff.SetPointError(ig, hTotUnc.GetBinError(ix,iy)  , 0.0)
                gdiff.SetLineWidth(4)
                gdiff.SetMarkerStyle(ROOT.kFullCircle)
                gdiff.SetLineColor(ROOT.kAzure+7)
                gdiff.SetMarkerColor(ROOT.kAzure+7)
                # partial uncertainty
                gdiffStatUncOnly.SetPoint(ig,      hTotUnc.GetBinContent(ix,iy),        ig+0.5)
                gdiffStatUncOnly.SetPointError(ig, hStatUncOnly.GetBinContent(ix,iy)  , 0.0)
                gdiffStatUncOnly.SetLineWidth(4)
                gdiffStatUncOnly.SetMarkerStyle(ROOT.kFullCircle)
                gdiffStatUncOnly.SetLineColor(ROOT.kPink-6)
                gdiffStatUncOnly.SetMarkerColor(ROOT.kPink-6)
                gdiffStatUncOnly.SetMarkerSize(1.5)
                ig += 1


        # draw histograms
        #adjustSettings_CMS_lumi()
        canvas = ROOT.TCanvas("canvas","",2200,1200)
        canvas.SetRightMargin(0.2)
        canvas.SetTickx(1)
        canvas.SetTicky(1)
        canvas.SetTopMargin(0.05)
        canvas.cd()
        setStyleTH2(hTotUnc,charge)
        hTotUnc.SetBarOffset(0.2)
        hTotUnc.Draw("COLZ0 TEXT45E")
        hStatUncOnly.SetMarkerColor(ROOT.kRed+2)
        hStatUncOnly.SetBarOffset(-0.2)
        hStatUncOnly.Draw("TEXT45 SAME")
        for ext in ['png', 'pdf']:
            canvas.SaveAs('{od}/{cn}_withStatUnc.{ext}'.format(od=inputDirNomi, cn=hname, ext=ext))


        # draw histograms

        # dummy histogram to draw frame of graph summary
        xmax = max(gdiff.GetX()[i] + 1.3*gdiff.GetErrorX(i) for i in xrange(nbinsClosure))
        xmin = min(gdiff.GetX()[i] - 1.3*gdiff.GetErrorX(i) for i in xrange(nbinsClosure))
        dx = 0.1*(xmax-xmin)
        frame = ROOT.TH2D("frame","", 100, xmin-dx, xmax-dx, nbinsClosure, 0., nbinsClosure)
        #frame.GetXaxis().SetTitle("#Deltam - #Deltam_{MC}  (GeV)");
        frame.GetXaxis().SetTitle("scale");
        frame.GetXaxis().SetNdivisions(505)
        ig = 1
        for iy in range(1,1+hTotUnc.GetNbinsY()):
            for ix in range(1,1+hTotUnc.GetNbinsX()):
                # Y axis label
                label = "%g < p_{T} < %g   %.1f < #eta < %.1f" % (hTotUnc.GetYaxis().GetBinLowEdge(iy),
                                                                  hTotUnc.GetYaxis().GetBinLowEdge(iy+1),
                                                                  hTotUnc.GetXaxis().GetBinLowEdge(ix),
                                                                  hTotUnc.GetXaxis().GetBinLowEdge(ix+1))
                frame.GetYaxis().SetBinLabel(ig, label)            
                ig += 1

        c2 = ROOT.TCanvas("c2","",1000,1200)
        c2.SetTickx(1)
        c2.SetTicky(1)
        c2.SetGrid()
        c2.SetRightMargin(0.04)
        c2.SetLeftMargin(0.3)
        c2.SetTopMargin(0.05)
        c2.cd()
        frame.Draw()
        gdiff.Draw("PZ SAME")
        gdiffStatUncOnly.Draw("PZ SAME")
        frame.GetXaxis().SetRangeUser(-0.001,0.001)
        gdiff.GetXaxis().SetRangeUser(-0.001,0.001)
        lat = ROOT.TLatex()
        lat.SetNDC()
        lat.SetTextFont(42)
        lat.SetTextSize(0.03)
        lat.DrawLatex(0.3, 0.96, '#bf{CMS} #it{Preliminary}')
        lat.DrawLatex(0.75, 0.96, '36.3 fb^{-1} (13 TeV)')
        leg = ROOT.TLegend(0.01,0.01,0.55,0.06)
        leg.SetBorderSize(0)
        leg.SetNColumns(2)
        leg.AddEntry(gdiff,"fit #oplus stat. unc.","PL")
        leg.AddEntry(gdiffStatUncOnly,"stat. unc.","PL")
        leg.Draw("same")
        c2.RedrawAxis("sameaxis")
        for ext in ['png', 'pdf']:
            c2.SaveAs('{od}/{cn}_withStatUnc_summary.{ext}'.format(od=inputDirNomi, cn=hname, ext=ext))
