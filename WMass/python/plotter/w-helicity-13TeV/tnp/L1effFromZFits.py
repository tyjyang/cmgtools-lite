#!/usr/bin/env python
import sys, os
import numpy as np
from math import *
from array import array
from CMGTools.WMass.plotter.mcPlots import doSpam,SAFE_COLOR_LIST

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

SAFE_STYLE_LIST = [ROOT.kOpenCircle,ROOT.kOpenSquare,ROOT.kOpenTriangleUp,ROOT.kOpenDiamond,ROOT.kOpenTriangleDown]

def doTinyCmsPrelimStandalone(textLeft="_default_",textRight="_default_",hasExpo=False,textSize=0.033,lumi=None, xoffs=0, options=None):
    if textLeft  == "_default_": textLeft  = "#bf{CMS} #it{Preliminary}"
    if textRight == "_default_": textRight = "%(lumi) (13 TeV)"
    if lumi      == None       : lumi      = 35.9
    if   lumi > 3.54e+1: lumitext = "%.0f fb^{-1}" % lumi
    elif lumi > 3.54e+0: lumitext = "%.1f fb^{-1}" % lumi
    elif lumi > 3.54e-1: lumitext = "%.2f fb^{-1}" % lumi
    elif lumi > 3.54e-2: lumitext = "%.0f pb^{-1}" % (lumi*1000)
    elif lumi > 3.54e-3: lumitext = "%.1f pb^{-1}" % (lumi*1000)
    else               : lumitext = "%.2f pb^{-1}" % (lumi*1000)
    lumitext = "%.1f fb^{-1}" % lumi
    textLeft = textLeft.replace("%(lumi)",lumitext)
    textRight = textRight.replace("%(lumi)",lumitext)
    if textLeft not in ['', None]:
        doSpam(textLeft, (.28 if hasExpo else .17)+xoffs, .955, .60+xoffs, .995, align=12, textSize=textSize)
    if textRight not in ['', None]:
        doSpam(textRight,.68+xoffs, .955, .99+xoffs, .995, align=12, textSize=textSize)

if __name__ == '__main__':

    if len(args) != 3:
        print "You must provide two tables: num and denom"
        exit(0)

    ROOT.gROOT.ProcessLine(".x ~/cpp/tdrstyle.cc")
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    outFileName = "plots/l1EG_eff.root"
    effFile = ROOT.TFile.Open(outFileName,'recreate')

    tfbxm1 = ROOT.TFile.Open(args[1])
    hm1 = tfbxm1.Get('bxm1_f_signal')
    hm1_eta = hm1.ProjectionX()
    tfbx0 = ROOT.TFile.Open(args[2])
    h0 = tfbx0.Get('bx0_f_signal')
    h0_eta = h0.ProjectionX()
    h0.Add(hm1)
    h0_eta.Add(hm1_eta)

    pEff = ROOT.TEfficiency(hm1,h0)
    pEff.SetStatisticOption(ROOT.TEfficiency.kFCP);
    pEff_eta = ROOT.TEfficiency(hm1_eta,h0_eta)
    pEff_eta.SetStatisticOption(ROOT.TEfficiency.kFCP);

    c1 = ROOT.TCanvas("c1","",600,600)
    c1.SetLogy()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(1);
    ROOT.gStyle.SetPaintTextFormat(".3f")
    pEff.Draw("COLZ0 TEXTE")
    doTinyCmsPrelimStandalone(lumi=35.9)
    c1.SaveAs("l1eff.pdf"); c1.SaveAs("l1eff.png")

    pEff_eta.SetTitle(";#eta;BX-1 probability")
    pEff_eta.Draw("AP")    
    ROOT.gPad.Update()
    graph = pEff_eta.GetPaintedGraph()
    graph.SetMinimum(0.001)
    graph.SetMaximum(1.5)
    func_plus  = ROOT.TF1("func_plus"  , "expo", 1.5, 3.0);
    func_minus = ROOT.TF1("func_minus" , "expo",-3.0,-1.5);
    graph.Fit('func_plus','R','sames')
    graph.Fit('func_minus','RS+','sames')

    x=ROOT.Double(0); ymin = ROOT.Double(0);  ymax = ROOT.Double(0);
    graph.GetPoint(0,x,ymin)
    graph.GetPoint(graph.GetN()-1,x,ymax)

    ROOT.gPad.Update(); 
    doTinyCmsPrelimStandalone(lumi=35.9)
    c1.SaveAs("plots/l1eff_eta.pdf"); c1.SaveAs("plots/l1eff_eta.png")

    effFile.cd()

    eff_th2 = ROOT.TH2D("l1EG_eff","",120,-3.0,3.0,1,30,100)
    for ieta in xrange(120):
        hi = eff_th2.GetXaxis().GetBinUpEdge(ieta+1)
        center = eff_th2.GetXaxis().GetBinCenter(ieta+1)
        if hi < -1.5:
            prefireProb = min(ymin,func_minus(center))
        elif hi < 1.5:
            prefireProb = 0
        else:
            prefireProb = min(ymax,func_plus(center))
        eff_th2.SetBinContent(ieta+1, 1, 1.0-prefireProb)
    pEff.Write()
    eff_th2.Write()
    effFile.Close()

    print "Smoothed pre-fire SF saved in ",outFileName
