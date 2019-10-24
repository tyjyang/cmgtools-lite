import os, sys, math
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

f_etaEE0p2 = "../postprocessing/data/leptonSF/new2016_madeSummer2018/systEff_trgel_pt30to55_etaEB0p1_etaEE0p2.root"
f_etaEE0p1 = "../postprocessing/data/leptonSF/new2016_madeSummer2018/systEff_trgel_pt30to45_etaEE0p1.root"
outdir = "plots/scaleFactors_Final/effSyst_fromRooFitResult_onlyStatUnc_el_merged/"
#f_merged   = "../postprocessing/data/leptonSF/new2016_madeSummer2018/systEff_trgel.root"
f_merged   = outdir + "systEff_trgel.root"

createPlotDirAndCopyPhp(outdir)

tf1 = ROOT.TFile.Open(f_etaEE0p2)
if not tf1 or not tf1.IsOpen():
    print "Error opening file %s" % f_etaEE0p2
    quit()
tf2 = ROOT.TFile.Open(f_etaEE0p1)
if not tf2 or not tf2.IsOpen():
    print "Error opening file %s" % f_etaEE0p1
    quit()

tfout = ROOT.TFile.Open(f_merged,"recreate")
if not tfout or not tfout.IsOpen():
    print "Error opening file %s" % f_merged
    quit()

h1 = {}
h2 = {}
hmerge = {}

for name in ["p0", "p1", "p2"]:
    h1[name] = tf1.Get(name)
    h1[name].SetDirectory(0)
    h2[name] = tf2.Get(name)
    h2[name].SetDirectory(0)
    npt = h1[name].GetNbinsY()
    ptl = h1[name].GetYaxis().GetBinLowEdge(1)
    pth = h1[name].GetYaxis().GetBinLowEdge(npt+1)
    hmerge[name] = ROOT.TH2F("hmerge_%s" % name,"Nuisance for Erf %s" % name, 50, -2.5, 2.5, npt, ptl, pth)
    hmerge[name].Reset("ICESM")
    for ix in range(1,hmerge[name].GetNbinsX()+1):
        for iy in range(1,hmerge[name].GetNbinsY()+1):
            if abs(hmerge[name].GetXaxis().GetBinCenter(ix)) > 1.7:
                content = h2[name].GetBinContent(ix,iy)
            else:
                content = h1[name].GetBinContent(ix,iy)
            hmerge[name].SetBinContent(ix, iy, content)
    drawCorrelationPlot(hmerge[name],"electron #eta","electron p_{T} [GeV]::30,55",
                        "variation / nominal" + ("::-0.004,0.001" if name == "p0" else "::-0.01,0.01"),
                        hmerge[name].GetName(),"ForceTitle",outdir,0,0,False,False,False,
                        1,palette=55)
    hmerge[name].Write(name)


tfout.Close()
print "Created file %s" % f_merged
tf2.Close()
tf1.Close()

