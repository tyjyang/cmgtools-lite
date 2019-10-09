#!/bin/env python

import ROOT, os, sys, re, array, math

from make_diff_xsec_cards import getDiffXsecBinning
from make_diff_xsec_cards import templateBinning

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gROOT.SetBatch(True)

import utilities
utilities = utilities.util()

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('--input-muon', dest='inputMuon', default='', type='string', help='input root files with all histograms for muons')
    parser.add_option('--input-electron', dest='inputElectron', default='', type='string', help='input root files with all histograms for electrons')
    parser.add_option('--input-combination', dest='inputCombination', default='', type='string', help='input root files with all histograms for combination')
    parser.add_option('-o','--outdir', dest='outdir', default='', type='string', help='output directory to save things')
    parser.add_option('-b','--binning', dest='binning', default='', type='string', help='File with binning')
    (options, args) = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH1.StatOverflows(True)

    outname = ""
    if options.outdir:
        outname = options.outdir
        addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
        createPlotDirAndCopyPhp(outname)
    else:
        print "Error: you should specify an output folder using option -o <name>. Exit"
        quit()

    if not options.binning:
        print "Error: you should specify the binning file using option -b <file>.Exit"
        quit()
    etaPtBinningVec = getDiffXsecBinning(options.binning, "gen")
    genBins  = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])
    netabins = genBins.Neta
    nptbins  = genBins.Npt

    names = ["hChAsymm1Dpt_LEPTON",
             "hDiffXsec_1Dpt_LEPTON_plus",
             "hDiffXsec_1Dpt_LEPTON_minus",
             "hChAsymm1Deta_LEPTON",
             "hDiffXsec_1Deta_LEPTON_plus",
             "hDiffXsec_1Deta_LEPTON_minus",
             "hDiffXsec_LEPTON_plus_unrollTo1D_pt",
             "hDiffXsec_LEPTON_minus_unrollTo1D_pt",
            "hDiffXsec_LEPTON_plus_unrollTo1D_eta",
             "hDiffXsec_LEPTON_minus_unrollTo1D_eta"] 

    fitTypes = ["muon", "electron", "lepton"]

    files = {"muon"     : options.inputMuon,
             "electron" : options.inputElectron,
             "lepton"   : options.inputCombination}
    
    hMuon = {}
    hElectron = {}
    hLepton = {}

    allHist = {"muon"     : hMuon, 
               "electron" : hElectron,
               "lepton"   : hLepton}

    adjustSettings_CMS_lumi()

    for lep in fitTypes:
        if files[lep] == "": continue
        print "="*30
        print "Getting histograms for %s" % lep
        print "-"*30
        f = ROOT.TFile(files[lep], 'read')
        if not f:
            print "Error: could not read file %s. Abort" % files[lep]
            quit()
        for name in names:
            newname = name.replace("LEPTON",lep)
            print newname
            allHist[lep][name] = f.Get(newname)
            if allHist[lep][name]: 
                allHist[lep][name].SetDirectory(0)
            else:
                print "Error in getting histogram %s. Abort" % newname
                quit()
        f.Close()
        print ""


    canvas = 0
    canvUnroll = ROOT.TCanvas("canvUnroll","",4500,1500) # 3000,1500
    canvas1D = ROOT.TCanvas("canvas1D","",800,800)

    leftMargin = 0.15 
    rightMargin = 0.04
    additionalText = ""
    additionalTextLatex = ""
    ptRangeText = "p_{{T}}^{{{l}}} #in [{ptmin:3g}, {ptmax:3g}] GeV".format(l="l", ptmin=genBins.ptBins[0], ptmax=genBins.ptBins[-1])
    etaRangeText = "|#eta|^{{{l}}} #in [{etamin:3g}, {etamax:3g}]".format(l=" l", etamin=genBins.etaBins[0], etamax=genBins.etaBins[-1])


    for name in names:
        chargeSign = "+" if "plus" in name else "-" if "minus" in name else ""
        xaxisTitle = 'gen lepton p_{T} [GeV]' if "1Dpt" in name else 'gen lepton |#eta|'
        yaxisTitle = ""
        if "ChAsym" in name: yaxisTitle = "Asymmetry::0,0.25"
        elif "DiffXsec" in name:
            yaxisTitle = "d#sigma / d%s " % "dp_{T}" if "1Dpt" in name else "|#eta|"
            if "Norm" in name: yaxisTitle += " / #sigma_{tot}"
            else             : yaxisTitle += " [pb/GeV]::0,340" if "1Dpt" in name else " [pb]"
        cname = name
        canvas = canvas1D
        legendCoords = ""
        texCoord = ""
        if "1Dpt" in name:
            legendCoords = "0.62,0.92,0.7,0.9"
            texCoord = "0.2,0.5" if "plus" in name else "0.2,0.5" if "minus" in name else "0.2,0.5"
        else:
            texCoord = "0.2,0.8" if "plus" in name else "0.2,0.5" if "minus" in name else "0.2,0.8"
            legendCoords = "0.62,0.92,0.4,0.6" if "ChAsym" in name else "0.62,0.92,0.7,0.9"

        additionalText = ""
        additionalTextLatex = "W^{{{chs}}} #rightarrow l#nu;{vartext}::{txc},0.08,0.04".format(chs=" "+chargeSign,
                                                                                               vartext=etaRangeText if "1Dpt" in name else ptRangeText,
                                                                                               txc=texCoord)

        varBinRanges = []
        vertLinesArg = ""

        if "unrollTo1D" in name:
            leftMargin =  0.06
            rightMargin = 0.02
            legendCoords = "0.2,0.5,0.84,0.90;3"
            if "ChAsym" in name: 
                yaxisTitle = "Asymmetry::0,0.25"
            elif "DiffXsec" in name:
                yaxisTitle = "d^{2}#sigma / d|#eta|dp_{T}"
                if "Norm" in name: 
                    yaxisTitle += " / #sigma_{tot} [1/GeV]"
                else: 
                    yaxisTitle += " [pb/GeV]::0,250"

            additionalText = "W^{{{chs}}} #rightarrow l#nu::0.8,0.84,0.9,0.9".format(chs=" "+chargeSign)
            additionalTextLatex = ""

            canvas = canvUnroll
            if "_eta" in name:
                #xaxisTitle = xaxisTitle + " = 1 + ipt + ieta * %d; ipt in [%d,%d], ieta in [%d,%d]" % (nptbins-1,0,nptbins-1,0,netabins-1)
                xaxisTitle = "cross section unrolled along |#eta|: |#eta| #in [%.1f, %.1f]" % (genBins.etaBins[0], genBins.etaBins[-1])
                vertLinesArg = "{a},{b}".format(a=genBins.Npt,b=genBins.Neta)
                for ipt in range(0,genBins.Npt):
                    varBinRanges.append("p_{{T}} #in [{ptmin:3g}, {ptmax:.3g}]".format(ptmin=genBins.ptBins[ipt], ptmax=genBins.ptBins[ipt+1]))
            else:
                #xaxisTitle = xaxisTitle + " = 1 + ieta + ipt * %d; ipt in [%d,%d], ieta in [%d,%d]" % (netabins-1,0,nptbins-1,0,netabins-1)
                xaxisTitle = "cross section unrolled along p_{T}: p_{T} #in [%.3g, %.3g] GeV" % (genBins.ptBins[0], genBins.ptBins[-1])
                vertLinesArg = "{a},{b}".format(a=genBins.Neta,b=genBins.Npt)
                for ieta in range(0,genBins.Neta):
                    varBinRanges.append("|#eta| #in [{etamin:.1f}, {etamax:.1f}]".format(etamin=genBins.etaBins[ieta], etamax=genBins.etaBins[ieta+1]))


        drawMuElComparison(hLepton[name], hMuon[name], hElectron[name],xaxisTitle,yaxisTitle,cname,
                           outname, labelRatioTmp="lep/comb.::0.8,1.2", legendCoords=legendCoords,
                           draw_both0_noLog1_onlyLog2=1, leftMargin=leftMargin, rightMargin=rightMargin, passCanvas=canvas, lumi=35.9,
                           drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.35, 
                           moreTextLatex=additionalTextLatex, moreText=additionalText)
