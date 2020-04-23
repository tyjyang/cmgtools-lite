#!/bin/env python

# make 2D and 1D plots for fakes from same-sign lepton region.
# needs as input a root file created with utilityMacros/src/loopNtuples_fakes.C


import ROOT, os, datetime, re, operator, math
from array import array

ROOT.gROOT.SetBatch(True)

import utilities
utilities = utilities.util()

import sys
sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)
    from optparse import OptionParser
    parser = OptionParser(usage='%prog th3.root [options] ')
    parser.add_option('-o','--outdir',     dest='outdir',     default='',   type='string', help='outdput directory to save plots')
    parser.add_option(     '--palette'  , dest='palette',      default=55, type=int, help='Set palette: use a negative number to select a built-in one, otherwise the default is 55 (kRainbow)')
    parser.add_option(     '--rebinEtaPt'  , dest='rebinEtaPt',      default=(0,0), nargs=2, type=int, help='Rebinnign factor for eta-pt distribution. Default is none, equivalent to 1,1')
    (options, args) = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    if len(args) < 1:
        parser.print_usage()
        quit()

    outdir = options.outdir
    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    # if options.outdir:
    #     ROOT.gROOT.SetBatch()
    #     if not os.path.isdir(options.outdir):
    #         os.system('mkdir -p {od}'.format(od=options.outdir))
    #     os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/m/mciprian/public/index.php',od=options.outdir))

    h2_parity_Mt_Mz_passIso = {}
    h2_parity_eta_pt_passIso = {}
    h2_parity_Mt_Mz_failIso = {}
    h2_parity_eta_pt_failIso = {}

    tf = ROOT.TFile.Open(args[0])            

    for parity in ["odd", "even"]:

        h3_MtMzIso = tf.Get("h3_%s_Mt_Mz_iso" % parity)
        nx    = h3_MtMzIso.GetNbinsX()
        xlow  = h3_MtMzIso.GetXaxis().GetBinLowEdge(1)
        xhigh = h3_MtMzIso.GetXaxis().GetBinLowEdge(1+nx)
        ny    = h3_MtMzIso.GetNbinsY()
        ylow  = h3_MtMzIso.GetYaxis().GetBinLowEdge(1)
        yhigh = h3_MtMzIso.GetYaxis().GetBinLowEdge(1+ny)
        nz    = h3_MtMzIso.GetNbinsZ() 
        h2_parity_Mt_Mz_passIso[parity] = ROOT.TH2F("h2_{p}_Mt_Mz_passIso".format(p=parity),"",
                                                    nx,xlow,xhigh,ny,ylow,yhigh)
        h2_parity_Mt_Mz_passIso[parity].SetDirectory(0)
        h2_parity_Mt_Mz_failIso[parity] = ROOT.TH2F("h2_{p}_Mt_Mz_failIso".format(p=parity),"",
                                                    nx,xlow,xhigh,ny,ylow,yhigh)
        h2_parity_Mt_Mz_failIso[parity].SetDirectory(0)
        binIso0p15 = h3_MtMzIso.GetZaxis().FindFixBin(0.1499) # slightly lower value to get bin from 0.14 to 0.15
        fillTH2fromTH3zrange(h2_parity_Mt_Mz_passIso[parity],h3_MtMzIso,0,binIso0p15)
        fillTH2fromTH3zrange(h2_parity_Mt_Mz_failIso[parity],h3_MtMzIso,1+binIso0p15,1+nz) # include overflow bin

        for whatmt in ["pass", "fail"]:
            for whatmz in ["pass", "fail"]:
 
                h3_etaPtIso = tf.Get("h3_{p}_eta_pt_iso_{mt}Mt{mz}Mz".format(p=parity,mt=whatmt,mz=whatmz))
                nx    = h3_etaPtIso.GetNbinsX()
                xlow  = h3_etaPtIso.GetXaxis().GetBinLowEdge(1)
                xhigh = h3_etaPtIso.GetXaxis().GetBinLowEdge(1+nx)
                ny    = h3_etaPtIso.GetNbinsY()
                ylow  = h3_etaPtIso.GetYaxis().GetBinLowEdge(1)
                yhigh = h3_etaPtIso.GetYaxis().GetBinLowEdge(1+ny)
                nz    = h3_etaPtIso.GetNbinsZ() 
                h2_parity_eta_pt_passIso[(parity,whatmt,whatmz)] = ROOT.TH2F("h2_{p}_eta_pt_passIso{mt}Mt{mz}Mz".format(p=parity,mt=whatmt,mz=whatmz),"",nx,xlow,xhigh,ny,ylow,yhigh)
                h2_parity_eta_pt_passIso[(parity,whatmt,whatmz)].SetDirectory(0)
                h2_parity_eta_pt_failIso[(parity,whatmt,whatmz)] = ROOT.TH2F("h2_{p}_eta_pt_failIso{mt}Mt{mz}Mz".format(p=parity,mt=whatmt,mz=whatmz),"",nx,xlow,xhigh,ny,ylow,yhigh)
                h2_parity_eta_pt_failIso[(parity,whatmt,whatmz)].SetDirectory(0)
                binIso0p15 = h3_etaPtIso.GetZaxis().FindFixBin(0.1499) # slightly lower value to get bin from 0.14 to 0.15
                fillTH2fromTH3zrange(h2_parity_eta_pt_passIso[(parity,whatmt,whatmz)],h3_etaPtIso,0,binIso0p15)
                fillTH2fromTH3zrange(h2_parity_eta_pt_failIso[(parity,whatmt,whatmz)],h3_etaPtIso,1+binIso0p15,1+nz) # include overflow bin
    tf.Close()    

    adjustSettings_CMS_lumi()
    canvas2D = ROOT.TCanvas("canvas2D","",900,900)


    for k in h2_parity_Mt_Mz_passIso.keys():
        print k
        title = "{p} pass_isolation".format(p=k)
        h2_parity_Mt_Mz_passIso[k].SetTitle(title)
        drawCorrelationPlot(h2_parity_Mt_Mz_passIso[k],"W-like transverse mass [GeV]","m_{ll} [GeV]","Events",
                            h2_parity_Mt_Mz_passIso[k].GetName(),"ForceTitle",outdir,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D)

    for k in h2_parity_Mt_Mz_failIso.keys():
        title = "{p} fail_isolation".format(p=k)
        h2_parity_Mt_Mz_failIso[k].SetTitle(title)
        drawCorrelationPlot(h2_parity_Mt_Mz_failIso[k],"W-like transverse mass [GeV]","m_{ll} [GeV]","Events",
                            h2_parity_Mt_Mz_failIso[k].GetName(),"ForceTitle",outdir,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D)


    rebinEta = options.rebinEtaPt[0]
    rebinPt = options.rebinEtaPt[1]

    for k in h2_parity_eta_pt_passIso.keys():
        title = "{p} {mt}_m_{{T}} {mz}_m_{{Z}} pass_iso".format(p=k[0],mt=k[1],mz=k[2])
        h2_parity_eta_pt_passIso[k].SetTitle(title)
        drawCorrelationPlot(h2_parity_eta_pt_passIso[k],"muon #eta","muon p_{T} [GeV]","Events",
                            h2_parity_eta_pt_passIso[k].GetName(),"ForceTitle",outdir,rebinEta,rebinPt,False,False,False,1,palette=options.palette,passCanvas=canvas2D)

    for k in h2_parity_eta_pt_failIso.keys():
        title = "{p} {mt}_m_{{T}} {mz}_m_{{Z}} fail_iso".format(p=k[0],mt=k[1],mz=k[2])
        h2_parity_eta_pt_failIso[k].SetTitle(title)
        drawCorrelationPlot(h2_parity_eta_pt_failIso[k],"muon #eta","muon p_{T} [GeV]","Events",
                            h2_parity_eta_pt_failIso[k].GetName(),"ForceTitle",outdir,rebinEta,rebinPt,False,False,False,1,palette=options.palette,passCanvas=canvas2D)

    for parity in ["even", "odd"]:
        h2_parity_eta_pt_passIso_failMtOrMz = h2_parity_eta_pt_passIso[(parity,"fail","fail")].Clone("h2_{p}_eta_pt_passIso_failMtOrMz".format(p=parity))
        h2_parity_eta_pt_passIso_failMtOrMz.Add(h2_parity_eta_pt_passIso[(parity,"fail","pass")])
        h2_parity_eta_pt_passIso_failMtOrMz.Add(h2_parity_eta_pt_passIso[(parity,"pass","fail")])
        title = "{p} fail_m_{{T}}_or_m_{{Z}} pass_iso".format(p=parity)
        h2_parity_eta_pt_passIso_failMtOrMz.SetTitle(title)
        drawCorrelationPlot(h2_parity_eta_pt_passIso_failMtOrMz,"muon #eta","muon p_{T} [GeV]","Events",
                            h2_parity_eta_pt_passIso_failMtOrMz.GetName(),"ForceTitle",outdir,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D)

        h2_parity_eta_pt_failIso_failMtOrMz = h2_parity_eta_pt_failIso[(parity,"fail","fail")].Clone("h2_{p}_eta_pt_failIso_failMtOrMz".format(p=parity))
        h2_parity_eta_pt_failIso_failMtOrMz.Add(h2_parity_eta_pt_failIso[(parity,"fail","pass")])
        h2_parity_eta_pt_failIso_failMtOrMz.Add(h2_parity_eta_pt_failIso[(parity,"pass","fail")])
        title = "{p} fail_m_{{T}}_or_m_{{Z}} fail_iso".format(p=parity)
        h2_parity_eta_pt_failIso_failMtOrMz.SetTitle(title)
        drawCorrelationPlot(h2_parity_eta_pt_failIso_failMtOrMz,"muon #eta","muon p_{T} [GeV]","Events",
                            h2_parity_eta_pt_failIso_failMtOrMz.GetName(),"ForceTitle",outdir,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D)
