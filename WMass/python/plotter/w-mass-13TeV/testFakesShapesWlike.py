#!/bin/env python

# example
# python w-mass-13TeV/testFakesShapesWlike.py -o plots/Wlike/TREE_4_WLIKE_MU/test_QCDbks_ss_onlyPromptMC/noDxy/compareShapes/compareEtaPtInRegions_rebinEtaPt_2_2/ --palette 57 --rebinEtaPt 2 2

import ROOT, os, datetime, re, operator, math
from array import array

ROOT.gROOT.SetBatch(True)

import utilities
utilities = utilities.util()

import sys
sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

inputfolder = "plots/Wlike/TREE_4_WLIKE_MU/test_QCDbks_ss_onlyPromptMC/noDxy/"
# key folder
#files_ss = {"ss_passIso_passMtMz" : "test_sameSignRegion_passIsoMtMz/",
#            "ss_failIso_failMtMz" : "test_sameSignRegion_failIso_failMtorMz/",
#            "ss_passIso_failMtMz" : "test_sameSignRegion_passIso_failMtorMz/",
#            "ss_failIso_passMtMz" : "test_sameSignRegion_failIso_passMtMz/",
#        }
files_ss = {"ss_passIso_passMtMz" : "testMt_sameSignRegion_AllEvents_passMtMz_passIso/plus/",
            "ss_failIso_failMtMz" : "testMt_sameSignRegion_AllEvents_failMtMz_failIso/plus/",
            "ss_passIso_failMtMz" : "testMt_sameSignRegion_AllEvents_failMtMz_passIso/plus/",
            "ss_failIso_passMtMz" : "testMt_sameSignRegion_AllEvents_passMtMz_failIso/plus/",
        }

files_os = {"os_failIso_failMtMz" : "testMt_oppositeSignRegion_oddEvents_failMtMz_failIso/plus/",
            "os_failIso_passMtMz" : "testMt_oppositeSignRegion_oddEvents_passMtMz_failIso/plus/",
            "os_passIso_failMtMz" : "testMt_oppositeSignRegion_oddEvents_failMtMz_passIso/plus/",
            #"os_passIso_passMtMz" : "test_oppositeSignRegion_passIso_passMtMz/" subtraction won't work, data<MC
        }

# need python >= 3.5
#files = {**files_ss, **files_os}
files = files_ss.copy()
files.update(files_os)

fname = "plots_zmm.root"
hname = "pt__eta"


if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option('-o','--outdir',     dest='outdir',     default='',   type='string', help='outdput directory to save plots')
    parser.add_option(     '--palette'  , dest='palette',      default=55, type=int, help='Set palette: use a negative number to select a built-in one, otherwise the default is 55 (kRainbow)')
    parser.add_option(     '--rebinEtaPt'  , dest='rebinEtaPt',      default=(0,0), nargs=2, type=int, help='Rebinnign factor for eta-pt distribution. Default is none, equivalent to 1,1')
    (options, args) = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    outdir = options.outdir
    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    rebinEta = options.rebinEtaPt[0]
    rebinPt = options.rebinEtaPt[1]

    hists2 = {}
    hists1eta = {}
    hists1pt = {}

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
        if not hbkg:
            print "Error fetching hbkg %s from file %s" % (hname_full,fileFullName)
            quit()
        hists2[f] = ROOT.TH2D("etaPt_{f}".format(f=f),"",
                              hdata.GetNbinsX(),hdata.GetXaxis().GetBinLowEdge(1),hdata.GetXaxis().GetBinLowEdge(1+hdata.GetNbinsX()),
                              hdata.GetNbinsY(),hdata.GetYaxis().GetBinLowEdge(1),hdata.GetYaxis().GetBinLowEdge(1+hdata.GetNbinsY()))
        hists2[f].SetDirectory(0)
        hists2[f].Add(hdata,hbkg,1.0,-1.0) 
        hists2[f].SetTitle(f)
        hists2[f].RebinX(rebinEta)
        hists2[f].RebinY(rebinPt)
        hists1eta[f] = hists2[f].ProjectionX("eta_{f}".format(f=f),0,-1,"e")
        hists1eta[f].SetDirectory(0)
        hists1pt[f]  = hists2[f].ProjectionY("pt_{f}".format(f=f),0,-1,"e")
        hists1pt[f].SetDirectory(0)
        tf.Close()
        drawCorrelationPlot(hists2[f],"muon #eta","muon p_{T} [GeV]","Events",
                            hists2[f].GetName(),"ForceTitle",outdir,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D)


    canvas1D = ROOT.TCanvas("canvas1D","",900,900)    
        
    h1sort_pt_ss = []
    h1sort_eta_ss = []
    legEntries_ss = []
    for f in sorted(files_ss.keys()):
        h1sort_eta_ss.append(hists1eta[f].Clone("clone_"+hists1eta[f].GetName()))
        h1sort_pt_ss.append(hists1pt[f].Clone("clone_"+hists1pt[f].GetName()))
        legEntries_ss.append(f)

    h1sort_pt_os = []
    h1sort_eta_os = []
    legEntries_os = []
    for f in sorted(files_os.keys()):
        h1sort_eta_os.append(hists1eta[f].Clone("clone_"+hists1eta[f].GetName()))
        h1sort_pt_os.append(hists1pt[f].Clone("clone_"+hists1pt[f].GetName()))
        legEntries_os.append(f)

    h1sort_pt_failIso = []
    h1sort_eta_failIso = []
    legEntries_failIso = []
    for f in sorted(files.keys()):
        if "_failIso" not in f: continue
        h1sort_eta_failIso.append(hists1eta[f].Clone("clone_failIso_"+hists1eta[f].GetName()))
        h1sort_pt_failIso.append(hists1pt[f].Clone("clone_failIso_"+hists1pt[f].GetName()))
        legEntries_failIso.append(f)

    h1sort_pt_failMtorMz = []
    h1sort_eta_failMtorMz = []
    legEntries_failMtorMz = []
    for f in sorted(files.keys()):
        if "_failMtMz" not in f: continue
        h1sort_eta_failMtorMz.append(hists1eta[f].Clone("clone_failMtorMz_"+hists1eta[f].GetName()))
        h1sort_pt_failMtorMz.append(hists1pt[f].Clone("clone_failMtorMz_"+hists1pt[f].GetName()))
        legEntries_failMtorMz.append(f)

    h1sort_pt_passIso = []
    h1sort_eta_passIso = []
    legEntries_passIso = []
    for f in sorted(files.keys()):
        if "_passIso" not in f: continue
        h1sort_eta_passIso.append(hists1eta[f].Clone("clone_passIso_"+hists1eta[f].GetName()))
        h1sort_pt_passIso.append(hists1pt[f].Clone("clone_passIso_"+hists1pt[f].GetName()))
        legEntries_passIso.append(f)

    h1sort_pt_passMtorMz = []
    h1sort_eta_passMtorMz = []
    legEntries_passMtorMz = []
    for f in sorted(files.keys()):
        if "_passMtMz" not in f: continue
        h1sort_eta_passMtorMz.append(hists1eta[f].Clone("clone_passMtorMz_"+hists1eta[f].GetName()))
        h1sort_pt_passMtorMz.append(hists1pt[f].Clone("clone_passMtorMz_"+hists1pt[f].GetName()))
        legEntries_passMtorMz.append(f)

    drawNTH1(h1sort_pt_ss,legEntries_ss,"muon p_{T} [GeV]","Events","compareFakesVsPt_ss",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_ss,legEntries_ss,"muon #eta","Events","compareFakesVsEta_ss",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_pt_os,legEntries_os,"muon p_{T} [GeV]","Events","compareFakesVsPt_os",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_os,legEntries_os,"muon #eta","Events","compareFakesVsEta_os",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_pt_failIso,legEntries_failIso,"muon p_{T} [GeV]","Events","compareFakesVsPt_failIso",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_failIso,legEntries_failIso,"muon #eta","Events","compareFakesVsEta_failIso",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_pt_failMtorMz,legEntries_failMtorMz,"muon p_{T} [GeV]","Events",
             "compareFakesVsPt_failMtorMz",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_failMtorMz,legEntries_failMtorMz,"muon #eta","Events",
             "compareFakesVsEta_failMtorMz",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_pt_passIso,legEntries_passIso,"muon p_{T} [GeV]","Events","compareFakesVsPt_passIso",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_passIso,legEntries_passIso,"muon #eta","Events","compareFakesVsEta_passIso",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_pt_passMtorMz,legEntries_passMtorMz,"muon p_{T} [GeV]","Events",
             "compareFakesVsPt_passMtorMz",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_passMtorMz,legEntries_passMtorMz,"muon #eta","Events",
             "compareFakesVsEta_passMtorMz",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

    ## normalized plots

    for i,h in enumerate(h1sort_eta_ss):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_pt_ss):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_eta_os):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_pt_os):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_eta_failIso):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_pt_failIso):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_eta_failMtorMz):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_pt_failMtorMz):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_eta_passIso):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_pt_passIso):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_eta_passMtorMz):
        h.Scale(1./h.Integral())
    for i,h in enumerate(h1sort_pt_passMtorMz):
        h.Scale(1./h.Integral())

    drawNTH1(h1sort_pt_ss,legEntries_ss,"muon p_{T} [GeV]","arbitrary units","compareFakesVsPt_ss_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_ss,legEntries_ss,"muon #eta","arbitrary units::0.0,0.16","compareFakesVsEta_ss_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_pt_os,legEntries_os,"muon p_{T} [GeV]","arbitrary units","compareFakesVsPt_os_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_os,legEntries_os,"muon #eta","arbitrary units::0.0,0.16","compareFakesVsEta_os_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_pt_failIso,legEntries_failIso,"muon p_{T} [GeV]","arbitrary units",
             "compareFakesVsPt_failIso_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_failIso,legEntries_failIso,"muon #eta","arbitrary units::0.0,0.16",
             "compareFakesVsEta_failIso_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_pt_failMtorMz,legEntries_failMtorMz,"muon p_{T} [GeV]","arbitrary units",
             "compareFakesVsPt_failMtorMz_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_failMtorMz,legEntries_failMtorMz,"muon #eta","arbitrary units::0.0,0.16",
             "compareFakesVsEta_failMtorMz_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)
      
    drawNTH1(h1sort_pt_passIso,legEntries_passIso,"muon p_{T} [GeV]","arbitrary units",
             "compareFakesVsPt_passIso_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_passIso,legEntries_passIso,"muon #eta","arbitrary units::0.0,0.16",
             "compareFakesVsEta_passIso_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_pt_passMtorMz,legEntries_passMtorMz,"muon p_{T} [GeV]","arbitrary units",
             "compareFakesVsPt_passMtorMz_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.35,0.85,0.7,0.9",passCanvas=canvas1D)

    drawNTH1(h1sort_eta_passMtorMz,legEntries_passMtorMz,"muon #eta","arbitrary units::0.0,0.16",
             "compareFakesVsEta_passMtorMz_norm",outdir,
             draw_both0_noLog1_onlyLog2=1,labelRatioTmp="X/first::0.1,2.0",
             legendCoords="0.2,0.85,0.7,0.9",passCanvas=canvas1D)

########
########
