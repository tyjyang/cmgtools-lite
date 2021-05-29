#!/usr/bin/env python3

# manipulate plots for MC truth efficiency (needs to make the ratio of numerator and denominator to derive the efficiencies)
# e.g. using folders in http://mciprian.web.cern.ch/mciprian/WMassAnalysis/testNanoAOD/MCtruthEfficiency/testW_perEra/plus/
# and also making ratios of different eras
# histogram is named <prefix>__<workingPoint>_<process>_<eraVFP> for the numerator
#                and <prefix>_<process>_<eraVFP> for the denominator        

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

def getMinMaxForSameWorkingPoint(histDict, args, excludeMax=None):
    ret = {}
    for (era,wp) in histDict.keys():            
        if args.effRange[0] < args.effRange[1]:
            ret[wp] = (args.effRange[0], args.effRange[1])
        else:
            hmin,hmax = getMinMaxHisto(histDict[(era,wp)], sumError=False, excludeMax=excludeMax)
            if wp in ret.keys():
                ret[wp] = (min(hmin,ret[wp][0]), max(hmax,ret[wp][1]))
            else:
                ret[wp] = (hmin, hmax)
    return ret


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfolder", type=str, nargs=1)
    parser.add_argument("outputfolder",   type=str, nargs=1)
    parser.add_argument("--hname", default="bareMuon_pt_eta", help="Root of histogram name inside root file")
    parser.add_argument("-w", "--working-points", dest="workingPoints", type=str, default="accept,idip,trig,trigNoBit,iso,idipANDtrig,idipANDisonotrig,idipANDtrigANDiso,idipANDtrigNoBit", help="Comma separated list of working points to fetch input histograms")
    parser.add_argument("-e", "--era",    type=str, default="B,C,D,E,F,BToF,G,H,GToH", help="Comma separated list of eras, which identify the input subfolder")
    parser.add_argument("-p", "--process", default="Wmunu_plus", type=str, help="Process to pick histogram")
    parser.add_argument("-r", "--eff-range", dest="effRange", default=(0,-1), type=float, nargs=2, help="Z axis range for efficiency plot")
    parser.add_argument("--pt-range", dest="ptRange", default=(0,-1), type=float, nargs=2, help="Pt axis range for efficiency plot")
    parser.add_argument(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (let it be odd, so the central strip is white if not using --abs-value and the range is symmetric)')
    parser.add_argument(     '--palette'  , dest='palette',      default=58, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette' , default=False , action='store_true',   help='Inverte color ordering in palette')
    parser.add_argument(     '--plot-ratio-era', dest='plotRatioToEra', type=str, default=None,   help='Plot ratios with respect to this era')
    parser.add_argument(     '--tnp-file', dest='tnpFile', type=str, default=None,   help='Root file with MC efficiencies from tag-and-probe')
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    inputfolder = args.inputfolder[0] 
    outfolder = args.outputfolder[0]
    createPlotDirAndCopyPhp(outfolder)

    workingPoints = [str(x) for x in args.workingPoints.split(',')]
    eras = [str(x) for x in args.era.split(',')]

    adjustSettings_CMS_lumi()    
    canvas = ROOT.TCanvas("canvas","",800,800)

    mcEff = {}
    for era in eras:
        
        eraVFP = "postVFP" if era in ["G", "H", "GToH"] else "preVFP"
        rfname = f"{inputfolder}/{era}/plots_test.root"
        rf = safeOpenFile(rfname)
        hnomi = safeGetObject(rf, f"{args.hname}_{args.process}_{eraVFP}")
        hwp = {}
        for wp in workingPoints:
            hwp[wp] = safeGetObject(rf, f"{args.hname}__{wp}_{args.process}_{eraVFP}")
            mcEff[(era,wp)] = copy.deepcopy(hwp[wp].Clone(f"mcTruthEff_{wp}_{era}"))
            mcEff[(era,wp)].Divide(hnomi)
            mcEff[(era,wp)].SetTitle(f"{args.process}: {era}")
            #print(f"{mcEff.GetName()} {mcEff.Integral()}")
        if rf.IsOpen():    
            rf.Close()
            
    xAxisName = "bare muon #eta"
    yAxisName = "bare muon p_{T} (GeV)"
    if args.ptRange[0] < args.ptRange[1]: 
        ymin = args.ptRange[0]
        ymax = args.ptRange[1]
        yAxisName += f"::{ymin},{ymax}"
        
    minmax = getMinMaxForSameWorkingPoint(mcEff, args)       

    rfoutname = f"{outfolder}/mcTruthEff.root"
    rfout = safeOpenFile(rfoutname, mode="RECREATE")
    
    for era in eras:
        for wp in workingPoints:
            zAxisName = f"{wp} MC efficiency::{minmax[wp][0]},{minmax[wp][1]}"    
            drawCorrelationPlot(mcEff[(era,wp)],
                                xAxisName,
                                yAxisName,
                                zAxisName,
                                mcEff[(era,wp)].GetName(), plotLabel="ForceTitle", outdir=outfolder,
                                draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)
            mcEff[(era,wp)].Write()

        
    if args.plotRatioToEra:

        refEra = args.plotRatioToEra
        mcEffRatio = {}
        for (era,wp) in mcEff.keys():            
            if era == refEra:
                continue
            ratioName = f"mcTruthEffRatio_{wp}_{era}over{refEra}"
            mcEffRatio[(era,wp)] = copy.deepcopy(mcEff[(era,wp)].Clone(ratioName))
            mcEffRatio[(era,wp)].Divide(mcEff[(refEra,wp)])
            
        minmax = getMinMaxForSameWorkingPoint(mcEffRatio, args, excludeMax=1.1)       
        for era in eras:
            if era == refEra:
                continue
            for wp in workingPoints:
                mcEffRatio[(era,wp)].SetTitle(f"{args.process}: {era}/{refEra}")
                zAxisName = f"{wp} MC eff({era})/eff({refEra})::{minmax[wp][0]},{minmax[wp][1]}"    
                drawCorrelationPlot(mcEffRatio[(era,wp)],
                                    xAxisName,
                                    yAxisName,
                                    zAxisName,
                                    mcEffRatio[(era,wp)].GetName(), plotLabel="ForceTitle", outdir=outfolder,
                                    draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                    invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)
                mcEffRatio[(era,wp)].Write()


    tnpEras = ["B", "C", "D", "E", "F_preVFP", "BtoF"] # to fetch histograms in TnP file
    if args.tnpFile:

        wpToTNP = {"idip" : "idip",
                   "idipANDtrig" : "trigger",
                   "idipANDtrigANDiso"  : "iso",
                   "idipANDisonotrig"  : "isonotrig"
        }
        
        rf = safeOpenFile(args.tnpFile)
        tnpEff = {}
        for tnpera in tnpEras:
            for wp in wpToTNP.keys():
                charge = "both" if wpToTNP[wp] != "trigger" else "plus" if "/plus/" in inputfolder else "minus"
                eraName = "F" if tnpera == "F_preVFP" else "BToF" if tnpera == "BtoF" else tnpera
                tnpEff[(eraName,wp)] = safeGetObject(rf, f"effMC_{wpToTNP[wp]}_{tnpera}_{charge}")
                #tnpEff[(eraName,wp)] = ROOT.TH2D(f"effMC_{wpToTNP[wp]}_{tnpera}_{charge}_rebin", "",
                #                                 48, -2.4, 2.4, 40, 25, 65)
                #mcTruth = ROOT.TH2D(f"mcTruth_{wpToTNP[wp]}_{tnpera}_{charge}_rebin", "",
                #                    48, -2.4, 2.4, 40, 25, 65)
                mcTruth = copy.deepcopy(tnpEff[(eraName,wp)].Clone(f"mcTruth_{wpToTNP[wp]}_{tnpera}_{charge}_rebin"))
                for ix in range(1, 1 + mcTruth.GetNbinsX()):
                    for iy in range(1, 1 +  mcTruth.GetNbinsY()):
                        ptval = tnpEff[(eraName,wp)].GetYaxis().GetBinCenter(iy)
                        ybin = mcEff[(eraName,wp)].GetYaxis().FindFixBin(ptval+0.001)
                        mcTruth.SetBinContent(ix, iy, mcEff[(eraName,wp)].GetBinContent(ix, ybin))
                        # need to transform TnP histograms into new ones with same binning as the MC truth ones
                ratioTnpOverMCtruth = copy.deepcopy(tnpEff[(eraName,wp)].Clone(f"ratioTnpOverMCtruth_{wp}_{eraName}"))
                ratioTnpOverMCtruth.Divide(mcTruth)
                ratioTnpOverMCtruth.SetTitle(f"TnP / MC truth: {eraName}")
                zAxisName = f"{wp} eff(TnP) / eff(MCtruth)"    
                drawCorrelationPlot(ratioTnpOverMCtruth,
                                    xAxisName,
                                    yAxisName,
                                    zAxisName,
                                    ratioTnpOverMCtruth.GetName(), plotLabel="ForceTitle", outdir=outfolder,
                                    draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                    invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)
                rfout.cd()
                ratioTnpOverMCtruth.Write()
        
    if rfout.IsOpen():
        rfout.Close()
        print()
        print(f"All histograms saved in {rfout.GetName()}")
        print()
