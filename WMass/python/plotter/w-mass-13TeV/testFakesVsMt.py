#!/usr/bin/env python3

# example for 2 mT bins with border at 40 GeV, and validation in signal region (needs existing plots)
# python w-mass-13TeV/testFakesVsMt.py plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion//fakeRateRegion_postVFP_plus_lowIso//plots_fakerate.root plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion//fakeRateRegion_postVFP_plus_highIso//plots_fakerate.root plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion/testFakesVsMt_0to40to60/ --palette 87 --rebin-x 3 --rebin-y 5 --mt-bin-edges 0,40,60 --test-file plots/testNanoAOD/WmassPlots_jetEta2p4_fixMassWeight_splitW/fakeRateRegion_postVFP_plus_systTH3/postprocessing/distributions_signalRegion/plots.root

import os, re, array, math
import time
import argparse

## safe batch mode                                 
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

sys.path.append(os.getcwd())
from cropNegativeTemplateBins import cropNegativeContent


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfileLowIso", type=str, nargs=1)
    parser.add_argument("inputfileHighIso", type=str, nargs=1)
    parser.add_argument("outputfolder",   type=str, nargs=1)
    #parser.add_argument("--hname", default="mt_pt_eta", help="Root of histogram name inside root file")
    parser.add_argument("-x", "--x-axis-name", dest="xAxisName", default="Muon #eta", help="x axis name")
    parser.add_argument("-y", "--y-axis-name", dest="yAxisName", default="Muon p_{T} (GeV)", help="y axis name")
    parser.add_argument("-z", "--z-axis-name", dest="zAxisName", default="m_{T} (GeV)", help="z axis name")
    parser.add_argument("--mt-bin-edges", dest="mtEdges", default="0,10,20,30,40,50,60", type=str, help="Comma-separated list of bin edges for mT")
    parser.add_argument("--mt-nominal-range", dest="mtNominalRange", default="0,40", type=str, help="Comma-separated list of 2 bin edges for mT, representing the nominal range, used to derive the correction using also option --mt-value-correction")
    parser.add_argument("--mt-value-correction", dest="mtValueCorrection", default=55.0, type=float, help="Value at high mT where to evaluate correction with respect to nominal range passed with --mt-nominal-range")
    parser.add_argument(     "--rebin-x", dest="rebinEta", default=1, type=int, help="To rebin x axis (eta)")
    parser.add_argument(     "--rebin-y", dest="rebinPt", default=1, type=int, help="To rebin y axis (pt)")
    parser.add_argument(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (let it be odd, so the central strip is white if not using --abs-value and the range is symmetric)')
    parser.add_argument(     '--palette'  , dest='palette',      default=55, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette' , default=False , action='store_true',   help='Inverte color ordering in palette')
    parser.add_argument("--test-file", dest="testFile", type=str, default=None, help="Optional file to test the correction: it has the prefit shapes in the signal region, as obtained from plotFakesTemplate.py")
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    #outfolder = args.outputfolder[0]
    outfolder = args.outputfolder[0]
    createPlotDirAndCopyPhp(outfolder)

    xAxisName = args.xAxisName
    yAxisName = args.yAxisName
    zAxisName = args.zAxisName
    
    histoLowIso = None
    fileLowIso = safeOpenFile(args.inputfileLowIso[0])
    datakey = [key.GetName() for key in fileLowIso.GetListOfKeys() if "data" in key.GetName()][0]
    histoLowIso = safeGetObject(fileLowIso, datakey)
    #dataLowIso1D = histoLowIso.ProjectionZ("data_mt_lowIso", 0, -1, 0, -1, "e")
    #mcLowIso1D = []
    for key in fileLowIso.GetListOfKeys():
        name = key.GetName()
        if "data" in name: continue
        obj = key.ReadObj()
        obj.SetDirectory(0)
        histoLowIso.Add(obj, -1.0)
        #mcLowIso1D.append(obj.ProjectionZ(f"{key}_mt_lowIso", 0, -1, 0, -1, "e"))
    fileLowIso.Close()
    
    fileHighIso = safeOpenFile(args.inputfileHighIso[0])
    datakey = [key.GetName() for key in fileHighIso.GetListOfKeys() if "data" in key.GetName()][0]
    histoHighIso = safeGetObject(fileHighIso, datakey)
    #dataHighIso1D = histoHighIso.ProjectionZ("data_mt_HighIso", 0, -1, 0, -1, "e")
    #mcHighIso1D = []
    for key in fileHighIso.GetListOfKeys():
        name = key.GetName()
        if "data" in name: continue
        obj = key.ReadObj()
        obj.SetDirectory(0)
        histoHighIso.Add(obj, -1.0)
        #mcHighIso1D.append(obj.ProjectionZ(f"{key}_mt_highIso", 0, -1, 0, -1, "e"))
    fileHighIso.Close()

    histoLowIso.Rebin3D(args.rebinEta, args.rebinPt)
    histoHighIso.Rebin3D(args.rebinEta, args.rebinPt)

    cropNegativeContent(histoLowIso)
    cropNegativeContent(histoHighIso)
    
    adjustSettings_CMS_lumi()    
    canvas = ROOT.TCanvas("canvas","",800,800)
    canvas1D = ROOT.TCanvas("canvas1D","",800,700)
    
    mtEdges = [round(int(x),1) for x in args.mtEdges.split(',')] 
    nMtBins = len(mtEdges) -1
    ratio = []
    for imt in range(nMtBins):
        lowEdge = mtEdges[imt]
        highEdge = mtEdges[imt+1]
        binStart = histoLowIso.GetZaxis().FindFixBin(lowEdge)
        binEnd = histoLowIso.GetZaxis().FindFixBin(highEdge+0.001) - 1 # bin up edges belong to "next" bin
        h2LowIso = getTH2fromTH3(histoLowIso, f"pt_eta_mt{lowEdge}to{highEdge}_lowIso", binStart, binEnd)
        h2HighIso = getTH2fromTH3(histoHighIso, f"pt_eta_mt{lowEdge}to{highEdge}_highIso", binStart, binEnd)
        cropNegativeContent(h2LowIso)
        cropNegativeContent(h2HighIso)

        h2LowIso.SetTitle("Low Iso: m_{T} #in [%d, %d]" % (lowEdge, highEdge))
        drawCorrelationPlot(h2LowIso,
                            xAxisName,
                            yAxisName,
                            "Events (data - MC)",
                            h2LowIso.GetName(), plotLabel="ForceTitle", outdir=outfolder,
                            draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                            invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)

        h2HighIso.SetTitle("High Iso: m_{T} #in [%d, %d]" % (lowEdge, highEdge))
        drawCorrelationPlot(h2HighIso,
                            xAxisName,
                            yAxisName,
                            "Events (data - MC)",
                            h2HighIso.GetName(), plotLabel="ForceTitle", outdir=outfolder,
                            draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                            invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)

        ratio.append(h2LowIso.Clone(f"fakerateFactor_mt{lowEdge}to{highEdge}"))
        ratio[imt].SetTitle("m_{T} #in [%d, %d]" % (lowEdge, highEdge))
        ratio[imt].Divide(h2HighIso)
        drawCorrelationPlot(ratio[imt],
                            xAxisName,
                            yAxisName,
                            "fakerate factor: N(iso) / N(not-iso)::0,3",
                            ratio[imt].GetName(), plotLabel="ForceTitle", outdir=outfolder,
                            draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                            invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)


    if args.mtNominalRange:
        lowEdge, highEdge = map(int, args.mtNominalRange.split(','))
        binStart = histoLowIso.GetZaxis().FindFixBin(lowEdge)
        binEnd = histoLowIso.GetZaxis().FindFixBin(highEdge+0.001) - 1 # bin up edges belong to "next" bin
        h2LowIso = getTH2fromTH3(histoLowIso, f"pt_eta_nominalmt{lowEdge}to{highEdge}_lowIso", binStart, binEnd)
        h2HighIso = getTH2fromTH3(histoHighIso, f"pt_eta_nominalmt{lowEdge}to{highEdge}_highIso", binStart, binEnd)
        cropNegativeContent(h2LowIso)
        cropNegativeContent(h2HighIso)
        nominalFakerateFactor = h2LowIso.Clone(f"nominalFakerateFactor_mt{lowEdge}to{highEdge}")
        nominalFakerateFactor.SetTitle("m_{T} #in [%d, %d]" % (lowEdge, highEdge))        
        nominalFakerateFactor.Divide(h2HighIso)

    # for ieta in range(1, 2+h2LowIso.GetNbinsX()):
    #     print(round(h2LowIso.GetXaxis().GetBinLowEdge(ieta),1))

    etaLow = round(0.01 + h2LowIso.GetXaxis().GetBinLowEdge(1), 1)
    etaHigh = round(0.01 + h2LowIso.GetXaxis().GetBinLowEdge(1+h2LowIso.GetNbinsX()), 1)
    ptLow = round(0.01 + h2LowIso.GetYaxis().GetBinLowEdge(1), 1)
    ptHigh = round(0.01 + h2LowIso.GetYaxis().GetBinLowEdge(1+h2LowIso.GetNbinsY()), 1)

    hFakerateFactorCorrection = ROOT.TH2D("fakerateFactorCorrection", "m_{T} > %d GeV" % int(args.mtNominalRange.split(',')[1]),
                                          h2LowIso.GetNbinsX(), round(etaLow,1), round(etaHigh,1),
                                          h2LowIso.GetNbinsY(), round(ptLow,1), round(ptHigh,1))
    
    # now preparing a summary for each pt bin
    ptCentralBin = h2LowIso.GetYaxis().FindFixBin(39.5)
    for ipt in range(1, 1+h2LowIso.GetNbinsY()):
        ptBinLow = int(h2LowIso.GetYaxis().GetBinLowEdge(ipt))
        ptBinHigh = int(h2LowIso.GetYaxis().GetBinLowEdge(ipt+1))
        fakerateFactor_vs_etaMt = ROOT.TH2D("fakerateFactor_vs_etaMt_pt%dto%d" % (ptBinLow, ptBinHigh),
                                            "Muon p_{T} #in [%d, %d] GeV" % (ptBinLow, ptBinHigh),
                                            h2LowIso.GetNbinsX(), round(etaLow,1), round(etaHigh,1),
                                            nMtBins, array("d", mtEdges)
                                           )

        outfolder1D = outfolder + "fakerateFactor_fits_pt%dto%d/" % (ptBinLow, ptBinHigh)
        createPlotDirAndCopyPhp(outfolder1D)

        
        for ieta in range(1, 1+fakerateFactor_vs_etaMt.GetNbinsX()):

            etaBinLowNoRound = fakerateFactor_vs_etaMt.GetXaxis().GetBinLowEdge(ieta)
            etaBinHighNoRound = fakerateFactor_vs_etaMt.GetXaxis().GetBinLowEdge(ieta+1)
            etaBinLow =  round(0.01 + etaBinLowNoRound, 1)
            etaBinHigh = round(0.01 + etaBinHighNoRound, 1)
            # print(f"ieta = {ieta}    {ptBinLow} < pt < {ptBinHigh}     {etaBinLow} < eta < {etaBinHigh}    {etaBinLow} < etaNoRound < {etaBinHigh}")
            hFRfactorVsMt = ROOT.TH1D(f"hFRfactorVsMt_ieta{ieta}_pt{ptBinLow}to{ptBinHigh}",
                                      "%.1f < #eta < %.1f, p_{T} #in [%d, %d] GeV" % (etaBinLow, etaBinHigh, ptBinLow, ptBinHigh),
                                      nMtBins, array("d", mtEdges))

            # to make easier computation of correction factor below
            hTmp = []
            
            for imt in range(1, 1+fakerateFactor_vs_etaMt.GetNbinsY()):
                binContent = ratio[imt-1].GetBinContent(ieta, ipt)
                binError = ratio[imt-1].GetBinError(ieta, ipt)
                fakerateFactor_vs_etaMt.SetBinContent(ieta, imt, binContent)
                fakerateFactor_vs_etaMt.SetBinError(ieta, imt, binError)
                hFRfactorVsMt.SetBinContent(imt, binContent)
                hFRfactorVsMt.SetBinError(  imt, binError)
                if nMtBins == 2:
                    hTmp.append(ROOT.TH1D(f"hTmp{imt}","",1,0,1))
                    hTmp[imt-1].SetBinContent(1, max(0.0, binContent))
                    hTmp[imt-1].SetBinError(  1, binError)
                
            textLatex = "%.1f < #eta < %.1f;p_{T} #in [%d, %d] GeV::0.2,0.3,0.1,0.045" % (etaBinLow, etaBinHigh, ptBinLow, ptBinHigh)
            if nMtBins > 2:
                fitFunc = drawSingleTH1withFit(hFRfactorVsMt, zAxisName, "Fakerate factor: N(iso) / N(not-iso)",
                                               hFRfactorVsMt.GetName(),
                                               outfolder1D, lowerPanelHeight=0.0, passCanvas=canvas1D, moreTextLatex=textLatex,
                                               legendCoords="0.64,0.96,0.69,0.93", fitRange="0,40", fitOptions="WLMFS+")
                valHighMt = fitFunc.Eval(args.mtValueCorrection)
                if valHighMt < 0:
                    printLine(marker=" ")
                    printLine()
                    print(f"Warning: ieta = {ieta},   ipt = {ipt},   FRF(mt={args.mtValueCorrection}) = {valHighMt}")
                    print("Setting FRF to 0.05!")
                    printLine()
                    printLine(marker=" ")
                    valHighMt = 0.05
                hFakerateFactorCorrection.SetBinContent(ieta, ipt, valHighMt / nominalFakerateFactor.GetBinContent(ieta, ipt))
                
            elif nMtBins == 2:
                hTmp[1].Divide(hTmp[0])
                hFakerateFactorCorrection.SetBinContent(ieta, ipt, hTmp[1].GetBinContent(1))
                hFakerateFactorCorrection.SetBinError(  ieta, ipt, hTmp[1].GetBinError(1))
                drawSingleTH1(hFRfactorVsMt, zAxisName, "Fakerate factor: N(iso) / N(not-iso)", hFRfactorVsMt.GetName(),
                              outfolder1D, lowerPanelHeight=0.0, passCanvas=canvas1D, moreTextLatex=textLatex,
                              legendCoords="0.64,0.96,0.77,0.93")
                
        drawCorrelationPlot(fakerateFactor_vs_etaMt,
                            yAxisName,
                            zAxisName,
                            "fakerate factor: N(iso) / N(not-iso)",
                            fakerateFactor_vs_etaMt.GetName(), plotLabel="ForceTitle", outdir=outfolder,
                            draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                            invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)

        # drawCorrelationPlot(fakerateFactor_vs_etaMt,
        #                     yAxisName,
        #                     zAxisName,
        #                     "abs. uncertainty on fakerate factor",
        #                     fakerateFactor_vs_etaMt.GetName()+"_absUnc", plotLabel="ForceTitle", outdir=outfolder,
        #                     draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
        #                     invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True, plotError=True)

        drawCorrelationPlot(fakerateFactor_vs_etaMt,
                            yAxisName,
                            zAxisName,
                            "rel. uncertainty on fakerate factor",
                            "relUnc_"+fakerateFactor_vs_etaMt.GetName(), plotLabel="ForceTitle", outdir=outfolder,
                            draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                            invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True, plotRelativeError=True)


        drawCorrelationPlot(hFakerateFactorCorrection,
                            xAxisName,
                            yAxisName,
                            "QCD template correction",
                            hFakerateFactorCorrection.GetName(), plotLabel="ForceTitle", outdir=outfolder,
                            draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                            invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)
        drawCorrelationPlot(hFakerateFactorCorrection,
                            xAxisName,
                            yAxisName,
                            "rel. unc. on QCD template correction",
                            "relUnc_"+hFakerateFactorCorrection.GetName(), plotLabel="ForceTitle", outdir=outfolder,
                            draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                            invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True, plotRelativeError=True)



        if args.testFile:

            outfolderCheck = outfolder + "checkCorrection_signalRegion/"
            createPlotDirAndCopyPhp(outfolderCheck)

            fshape = safeOpenFile(args.testFile)
            hdata = safeGetObject(fshape, "data2D")
            hmc   = safeGetObject(fshape, "sigAndBkgNoFakes2D")
            hqcd  = safeGetObject(fshape, "data_fakes")
            fshape.Close()

            # apply correction to qcd
            for ix in range(1, 1+hqcd.GetNbinsX()):
                for iy in range(1, 1+hqcd.GetNbinsY()):
                    xCen = hqcd.GetXaxis().GetBinCenter(ix)
                    yCen = hqcd.GetYaxis().GetBinCenter(iy)
                    xBinCorr = max(1, min(hFakerateFactorCorrection.GetNbinsX(), hFakerateFactorCorrection.GetXaxis().FindFixBin(xCen)))
                    yBinCorr = max(1, min(hFakerateFactorCorrection.GetNbinsY(), hFakerateFactorCorrection.GetYaxis().FindFixBin(yCen)))
                    corr = hFakerateFactorCorrection.GetBinContent(xBinCorr, yBinCorr)
                    hqcd.SetBinContent(ix, iy, hqcd.GetBinContent(ix, iy) * corr)
                    hqcd.SetBinError(ix, iy, hqcd.GetBinError(ix, iy) * corr)

            hqcd.SetMarkerSize(0)
            hmc.SetMarkerSize(0)
            hdata_eta = hdata.ProjectionX("data_eta", 0, -1, "e")
            hmc_eta   = hmc.ProjectionX("mc_eta", 0, -1, "e")
            hqcd_eta   = hqcd.ProjectionX("qcd_eta", 0, -1, "e")
            hdata_pt  = hdata.ProjectionY("data_pt", 0, -1, "e")
            hmc_pt    = hmc.ProjectionY("mc_pt", 0, -1, "e")
            hqcd_pt    = hqcd.ProjectionY("qcd_pt", 0, -1, "e")

            stack_eta = ROOT.THStack("stack_eta", "signal and backgrounds")
            hmc_eta.SetFillColor(ROOT.kRed+2)
            hmc_eta.SetLineColor(ROOT.kBlack)
            hmc_eta.SetMarkerSize(0)
            stack_eta.Add(hmc_eta)
            hqcd_eta.SetFillColor(ROOT.kGray)
            hqcd_eta.SetLineColor(ROOT.kBlack)
            hqcd_eta.SetMarkerSize(0)
            stack_eta.Add(hqcd_eta)

            stack_pt = ROOT.THStack("stack_pt", "signal and backgrounds")
            hmc_pt.SetFillColor(ROOT.kRed+2)
            hmc_pt.SetLineColor(ROOT.kBlack)
            hmc_pt.SetMarkerSize(0)
            stack_pt.Add(hmc_pt)
            hqcd_pt.SetFillColor(ROOT.kGray)
            hqcd_pt.SetLineColor(ROOT.kBlack)
            hqcd_pt.SetMarkerSize(0)
            stack_pt.Add(hqcd_pt)

            legend = ROOT.TLegend(0.2,0.72,0.95,0.92)
            legend.SetFillColor(0)
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            legend.SetNColumns(2)
            legend.AddEntry(hdata_eta, "Data", "EP")
            legend.AddEntry(hmc_eta, "MC prompt", "F")
            legend.AddEntry(hqcd_eta, "QCD multijet", "F")

            canvas1Dproj = ROOT.TCanvas("canvas1Dproj", "", 800, 900)

            drawTH1dataMCstack(hdata_eta, stack_eta, xAxisName, "Events", "muon_eta_withCorr",
                               outfolderCheck, legend, ratioPadYaxisNameTmp="Data/MC::0.92,1.08", passCanvas=canvas1Dproj,
                               drawLumiLatex=True, xcmsText=0.3, noLegendRatio=True
            )
            drawTH1dataMCstack(hdata_pt, stack_pt, yAxisName, "Events", "muon_pt_withCorr",
                               outfolderCheck, legend, ratioPadYaxisNameTmp="Data/MC::0.92,1.08", passCanvas=canvas1Dproj,
                               drawLumiLatex=True, xcmsText=0.3, noLegendRatio=True
            )


