#!/usr/bin/env python3

# manipulate plots for MC truth efficiency (needs to make the ratio of numerator and denominator to derive the efficiencies)
# e.g. using folders in http://mciprian.web.cern.ch/mciprian/WMassAnalysis/testNanoAOD/MCtruthEfficiency/testW_perEra/plus/
# and also making ratios of different eras
# histogram is named <prefix>__<workingPoint>_<process>_<eraVFP> for the numerator
#                and <prefix>_<process>_<eraVFP> for the denominator        

# versus pt-eta (with PU weights and all eras)
# python w-mass-13TeV/makeMCtruthEffStudy.py plots/testNanoAOD/MCtruthEfficiency/W_perEra_noRecoAccept_finePt/plus/ --palette 87 --plot-ratio-era H --tnp-file testMuonSF/2021-05-31_allSFs_nodz_dxybs.root --rebin-pt 4 -w 'trackerOrGlobal,tracker,standalone,global,accept,idip,trig,trigNoBit,iso,idipANDtrig,idipANDisonotrig,idipANDtrigANDiso,idipANDtrigNoBit'

# versus pt-eta (no PU weights here, comparing only pre and postVFP)
# python w-mass-13TeV/makeMCtruthEffStudy.py plots/testNanoAOD/MCtruthEfficiency/W_perEra_noRecoAccept_finePt_noPUweights/plus/ --palette 87 --plot-ratio-era GToH --tnp-file testMuonSF/2021-05-31_allSFs_nodz_dxybs.root --rebin-pt 4 -w 'trackerOrGlobalAndStandalone,trackerOrGlobal,tracker,standalone,global,accept,idip,trig,trigNoBit,iso,idipANDtrig,idipANDisonotrig,idipANDtrigANDiso,idipANDtrigNoBit' -e "BToF,GToH"

# versus phi-eta (no PU weights here, comparing only pre and postVFP)
# python w-mass-13TeV/makeMCtruthEffStudy.py plots/testNanoAOD/MCtruthEfficiency/W_perEra_noRecoAccept_EtaPhi_noPUweights/plus/ --palette 87 --plot-ratio-era GToH -w 'trackerOrGlobalAndStandalone,trackerOrGlobal,tracker,standalone,global,accept,idip,trig,trigNoBit,iso,idipANDtrig,idipANDisonotrig,idipANDtrigANDiso,idipANDtrigNoBit' -e "BToF,GToH" -y "bare muon #phi" --hname "bareMuon_phi_eta"

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

def getMinMaxForSameWorkingPoint(histDict, args, excludeMax=None, excludeMin=None):
    ret = {}
    for (era,wp) in histDict.keys():            
        if args.effRange[0] < args.effRange[1]:
            ret[wp] = (args.effRange[0], args.effRange[1])
        else:
            hmin,hmax = getMinMaxHisto(histDict[(era,wp)], sumError=False, excludeMax=excludeMax, excludeMin=excludeMin)
            if wp in ret.keys():
                ret[wp] = (min(hmin,ret[wp][0]), max(hmax,ret[wp][1]))
            else:
                ret[wp] = (hmin, hmax)
    return ret


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfolder", type=str, nargs=1)
    # output folder is better chosen as input plus a subfolder
    #parser.add_argument("outputfolder",   type=str, nargs=1)
    parser.add_argument("--hname", default="bareMuon_pt_eta", help="Root of histogram name inside root file")
    parser.add_argument("-y", "--y-axis-name", dest="yAxisName", default="bare muon p_{T} (GeV)", help="y axis name")
    parser.add_argument("--postfix", default="", help="Postfix for output folder")
    parser.add_argument("-w", "--working-points", dest="workingPoints", type=str, default="trackerOrGlobalAndStandalone,trackerOrGlobal,tracker,standalone,global,accept,idip,trig,trigNoBit,iso,idipANDtrig,idipANDisonotrig,idipANDtrigANDiso,idipANDtrigNoBit", help="Comma separated list of working points to fetch input histograms")
    parser.add_argument("-e", "--era",    type=str, default="B,C,D,E,F,BToF,G,H,GToH", help="Comma separated list of eras, which identify the input subfolder")
    parser.add_argument("-p", "--process", default="Wmunu_plus", type=str, help="Process to pick histogram")
    parser.add_argument("-r", "--eff-range", dest="effRange", default=(0,-1), type=float, nargs=2, help="Z axis range for efficiency plot")
    parser.add_argument("--pt-range", dest="ptRange", default=(0,-1), type=float, nargs=2, help="Pt axis range for efficiency plot")
    parser.add_argument(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (let it be odd, so the central strip is white if not using --abs-value and the range is symmetric)')
    parser.add_argument(     '--palette'  , dest='palette',      default=58, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette' , default=False , action='store_true',   help='Inverte color ordering in palette')
    parser.add_argument(     '--plot-ratio-era', dest='plotRatioToEra', type=str, default=None,   help='Plot ratios with respect to this era')
    parser.add_argument(     '--tnp-file', dest='tnpFile', type=str, default=None,   help='Root file with MC efficiencies from tag-and-probe')
    #parser.add_argument(     '--tnp-ratio-eff-range', dest="tnpRatioEffRange", default=(0,-1), type=float, nargs=2, help="Z axis range for efficiency plot")
    #parser.add_argument(     "--rebin-pt", dest="rebinPt", default="10,12,14,16,18,20,22,24,25,27.5,30,32,34,36,38,40,42,44,47,50,55,65", type=str, help="To rebin pt axis")
    parser.add_argument(     "--rebin-pt", dest="rebinPt", default=1, type=int, help="To rebin pt axis")
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    #newBinEdges = [float(x) for x in args.rebinPt.split(",")] if len(args.rebinPt) else None
        
    inputfolder = args.inputfolder[0] 
    #outfolder = args.outputfolder[0]
    outfolder = args.inputfolder[0] + (f"/results_{args.postfix}/" if args.postfix else "/results/")
    createPlotDirAndCopyPhp(outfolder)

    workingPoints = [str(x) for x in args.workingPoints.split(',')]
    eras = [str(x) for x in args.era.split(',')]

    adjustSettings_CMS_lumi()    
    canvas = ROOT.TCanvas("canvas","",800,800)

    # to compare to TNP I need to reproduce P(A|B), while for MC truth I had P(A & B|bareMuon) or P(A|bareMuon)
    # so e.g. for idipANDtrigANDiso I can normalize to idipANDtrig to emulate P(iso|trigger & idip)
    # right hand side is the tnp keyword, left hand side is the way it can be obtained from available MC truth histograms
    wpToTNP = {"idip/trackerOrGlobal" : "idip",
               "idipANDtrig/idip" : "trigger",
               "idipANDtrigANDiso/idipANDtrig"  : "iso",
               "idipANDisonotrig/idip"  : "isonotrig"
    }
    
    mcEff = {}
    for era in eras:
        
        eraVFP = "postVFP" if era in ["F_postVFP", "G", "H", "GToH"] else "preVFP"
        rfname = f"{inputfolder}/{era}/plots_test.root"
        rf = safeOpenFile(rfname)
        hnomi = safeGetObject(rf, f"{args.hname}_{args.process}_{eraVFP}")
        #if newBinEdges != None:
        #    hnomi.RebinY(len(newBinEdges)-1,"",array('d',newBinEdges))
        hnomi.RebinY(args.rebinPt)
        hwp = {}
        for wp in workingPoints:
            hwp[wp] = safeGetObject(rf, f"{args.hname}__{wp}_{args.process}_{eraVFP}")
            mcEff[(era,wp)] = copy.deepcopy(hwp[wp].Clone(f"mcTruthEff_{wp}_{era}"))
            print(f"{args.hname}__{wp}_{args.process}_{eraVFP}")
            #if newBinEdges != None:
            #    mcEff[(era,wp)].RebinY(len(newBinEdges)-1,"",array('d',newBinEdges))
            mcEff[(era,wp)].RebinY(args.rebinPt)
            mcEff[(era,wp)].Divide(hnomi)
            mcEff[(era,wp)].SetTitle(f"{args.process}: {era}")
            #print(f"{mcEff.GetName()} {mcEff.Integral()}")    
        #
        # add other working points according to TnP steps, so e.g.
        # P(iso) means P(iso|trigger&idip)
        # P(idip) means P(idip|global OR tracker)
        print()
        if all(x in workingPoints for x in ["trackerOrGlobalAndStandalone", "standalone"]):
            print("Adding new working point (isTracker OR isGlobal | isStandalone)")
            mcEff[(era,"trackOrGlobGivenStandAlone")] = copy.deepcopy(mcEff[(era,"trackerOrGlobalAndStandalone")].Clone(f"mcTruthEff_trackOrGlobGivenStandAlone_{era}"))
            mcEff[(era,"trackOrGlobGivenStandAlone")].Divide(mcEff[(era,"standalone")])
            print()
        for wp in wpToTNP.keys():
            if wpToTNP[wp] not in workingPoints:
                continue
            num = str(wp.split('/')[0])
            den = str(wp.split('/')[1])
            tnpStep = wpToTNP[wp]
            hnum = safeGetObject(rf, f"{args.hname}__{num}_{args.process}_{eraVFP}")
            hden = safeGetObject(rf, f"{args.hname}__{den}_{args.process}_{eraVFP}")
            hnum.RebinY(args.rebinPt)
            hden.RebinY(args.rebinPt)
            mcEff[(era,f"{tnpStep}StepOfTnP")] = copy.deepcopy(hnum.Clone(f"mcTruthEff_{tnpStep}StepOfTnP_{era}"))
            mcEff[(era,f"{tnpStep}StepOfTnP")].Divide(hden)
            mcEff[(era,f"{tnpStep}StepOfTnP")].SetTitle(f"{args.process}: {era}")
        if rf.IsOpen():    
            rf.Close()

    workingPoints = [str(wp) for era,wp in mcEff.keys() if era == eras[0]] # add the new ones for the TnP steps
    
    xAxisName = "bare muon #eta"
    yAxisName = args.yAxisName
    if args.ptRange[0] < args.ptRange[1]: 
        ymin = args.ptRange[0]
        ymax = args.ptRange[1]
        yAxisName += f"::{ymin},{ymax}"
        
    minmax = getMinMaxForSameWorkingPoint(mcEff, args, excludeMax=1.0)       
    if "isoStepOfTnP" in minmax:
        minmax["isoStepOfTnP"] = (0.85, minmax["isoStepOfTnP"][1]) # customize some working points

    rfoutname = f"{outfolder}/mcTruthEff.root"
    rfout = safeOpenFile(rfoutname, mode="RECREATE")

    outfolder_eff = f"{outfolder}/efficiency/" 
    createPlotDirAndCopyPhp(outfolder_eff)

    mcTrigBitOverNoBit = {}
    for era in eras:
        for wp in workingPoints: 
            if "StepOfTnP" in wp:
                zAxisName = f"TnP step: {wp} MC efficiency::{minmax[wp][0]},{minmax[wp][1]}"    
            else:
                zAxisName = f"{wp} MC efficiency::{minmax[wp][0]},{minmax[wp][1]}"    
            drawCorrelationPlot(mcEff[(era,wp)],
                                xAxisName,
                                yAxisName,
                                zAxisName,
                                mcEff[(era,wp)].GetName(), plotLabel="ForceTitle", outdir=outfolder_eff,
                                draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)
            mcEff[(era,wp)].Write()
        if all((era, x) in mcEff.keys() for x in ["trig", "trigNoBit"]):
            mcTrigBitOverNoBit[(era,"trigBitOverNoBit")] = copy.deepcopy(mcEff[(era,"trig")].Clone(f"mcTruthEffRatio_trigBitOverNoBit_{era}"))
            mcTrigBitOverNoBit[(era,"trigBitOverNoBit")].Divide(mcEff[(era,"trigNoBit")])

        
    if all((era, x) in mcEff.keys() for x in ["trig", "trigNoBit"]):
        outfolder_effCheckTrigBit = f"{outfolder}/efficiency_checkTriggerBit/" 
        createPlotDirAndCopyPhp(outfolder_effCheckTrigBit)
        minmax = getMinMaxForSameWorkingPoint(mcTrigBitOverNoBit, args, excludeMin=0.99, excludeMax=1.01)    
        thisminz = minmax["trigBitOverNoBit"][0]
        thismaxz = minmax["trigBitOverNoBit"][1]
        for era in eras:
            zAxisName = f"MC eff(trig)/eff(trigNoBit)::{thisminz},{thismaxz}"    
            mcTrigBitOverNoBit[(era,"trigBitOverNoBit")].SetTitle(f"{args.process}: {era}")
            drawCorrelationPlot(mcTrigBitOverNoBit[(era,"trigBitOverNoBit")],
                                xAxisName,
                                yAxisName,
                                zAxisName,
                                mcTrigBitOverNoBit[(era,"trigBitOverNoBit")].GetName(), plotLabel="ForceTitle",
                                outdir=outfolder_effCheckTrigBit,
                                draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)


    if args.plotRatioToEra:

        refEra = args.plotRatioToEra

        mcEffRatio = {}
        for (era,wp) in mcEff.keys():            
            if era == refEra:
                continue
            ratioName = f"mcTruthEffRatio_{wp}_{era}over{refEra}"
            mcEffRatio[(era,wp)] = copy.deepcopy(mcEff[(era,wp)].Clone(ratioName))
            mcEffRatio[(era,wp)].Divide(mcEff[(refEra,wp)])

        outfolder_effRatioToEra = f"{outfolder}/efficiency_MCtruthRatioTo{refEra}/" 
        createPlotDirAndCopyPhp(outfolder_effRatioToEra)
        
        minmax = getMinMaxForSameWorkingPoint(mcEffRatio, args, excludeMax=1.1)       
        # small hack because numerator histograms for trigger efficiency are almost empty for pt < 24
        if "idipANDtrigANDiso" in minmax: minmax["idipANDtrigANDiso"] = (0.85, minmax["idipANDtrigANDiso"][1]) 
        if "idipANDtrig" in minmax:       minmax["idipANDtrig"]       = (0.85, minmax["idipANDtrig"][1]) 
        if "idipANDtrigNoBit" in minmax:  minmax["idipANDtrigNoBit"]  = (0.85, minmax["idipANDtrigNoBit"][1])
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
                                    mcEffRatio[(era,wp)].GetName(), plotLabel="ForceTitle", outdir=outfolder_effRatioToEra,
                                    draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                    invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)
                mcEffRatio[(era,wp)].Write()


    tnpErasTmp = ["B", "C", "D", "E", "F_preVFP", "G", "H", "BtoF", "GtoH"] # to fetch histograms in TnP file
    tnpEras = [x for x in tnpErasTmp if x.replace("_preVFP","F").replace("to","To") in eras]
    if args.tnpFile:
        
        rf = safeOpenFile(args.tnpFile)
        tnpEff = {}
        ratioTnpOverMCtruth = {}
        for tnpera in tnpEras:
            for wptmp in wpToTNP.keys():
                wp = f"{wpToTNP[wptmp]}StepOfTnP"
                charge = "both" if wpToTNP[wptmp] != "trigger" else "plus" if "/plus/" in inputfolder else "minus"
                eraName = "F" if tnpera == "F_preVFP" else "BToF" if tnpera == "BtoF" else "GToH" if tnpera == "GtoH" else tnpera
                tnpEff[(eraName,wp)] = safeGetObject(rf, f"effMC_{wpToTNP[wptmp]}_{tnpera}_{charge}")
                # clone from tnp histogram to preserve binning but then refill later
                mcTruth = copy.deepcopy(tnpEff[(eraName,wp)].Clone(f"mcTruth_{wpToTNP[wptmp]}StepOfTnP_{tnpera}_{charge}_rebin"))
                for ix in range(1, 1 + mcTruth.GetNbinsX()):
                    for iy in range(1, 1 +  mcTruth.GetNbinsY()):
                        ptval = tnpEff[(eraName,wp)].GetYaxis().GetBinCenter(iy)
                        ybin = mcEff[(eraName,wp)].GetYaxis().FindFixBin(ptval+0.001)
                        mcTruth.SetBinContent(ix, iy, mcEff[(eraName,wp)].GetBinContent(ix, ybin))
                        # need to transform TnP histograms into new ones with same binning as the MC truth ones
                ratioTnpOverMCtruth[(eraName,wp)] = copy.deepcopy(tnpEff[(eraName,wp)].Clone(f"ratioTnpOverMCtruth_{wp}_{eraName}"))
                ratioTnpOverMCtruth[(eraName,wp)].Divide(mcTruth)
                ratioTnpOverMCtruth[(eraName,wp)].SetTitle(f"TnP / MC truth: {eraName}")

                # if wpnorm != None:
                #     mcTruthNorm = copy.deepcopy(tnpEff[(eraName,wp)].Clone(f"mcTruth_{wpToTNP[wptmp]}_{tnpera}_{charge}_rebin_norm"))
                #     for ix in range(1, 1 + mcTruthNorm.GetNbinsX()):
                #         for iy in range(1, 1 +  mcTruthNorm.GetNbinsY()):
                #             ptval = tnpEff[(eraName,wp)].GetYaxis().GetBinCenter(iy)
                #             ybin = mcEff[(eraName,wpnorm)].GetYaxis().FindFixBin(ptval+0.001)
                #             mcTruthNorm.SetBinContent(ix, iy, mcEff[(eraName,wpnorm)].GetBinContent(ix, ybin))
                #             # need to transform TnP histograms into new ones with same binning as the MC truth ones
                #     ratioTnpOverMCtruth[(eraName,wp)].Multiply(mcTruthNorm)
                
        minmax = getMinMaxForSameWorkingPoint(ratioTnpOverMCtruth, args, excludeMax=2)       
        for wptmp in wpToTNP.keys():

            wp = f"{wpToTNP[wptmp]}StepOfTnP"
            outfolder_effRatioTnpOverMC = f"{outfolder}/efficiency_ratioTNPoverMCtruth/stepTnP_{wpToTNP[wptmp]}" 
            createPlotDirAndCopyPhp(outfolder_effRatioTnpOverMC)
            
            for tnpera in tnpEras:    
                eraName = "F" if tnpera == "F_preVFP" else "BToF" if tnpera == "BtoF" else "GToH" if tnpera == "GtoH" else tnpera
                zAxisName = f"TnP step = {wpToTNP[wptmp]}: eff(TnP) / eff(MCtruth)::{minmax[wp][0]},{minmax[wp][1]}"    
                drawCorrelationPlot(ratioTnpOverMCtruth[(eraName,wp)],
                                    xAxisName,
                                    yAxisName,
                                    zAxisName,
                                    ratioTnpOverMCtruth[(eraName,wp)].GetName(), plotLabel="ForceTitle",
                                    outdir=outfolder_effRatioTnpOverMC,
                                    draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                    invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)
                rfout.cd()
                ratioTnpOverMCtruth[(eraName,wp)].Write()

                if args.plotRatioToEra:    
                    refEra = args.plotRatioToEra
                    if eraName != refEra:
                        ratioTnpOverMCtruth_ratioOverRefEra = copy.deepcopy(ratioTnpOverMCtruth[(eraName,wp)])
                        ratioTnpOverMCtruth_ratioOverRefEra.SetTitle(f"TnP / MC truth: {eraName}/{refEra}")
                        ratioTnpOverMCtruth_ratioOverRefEra.Divide(ratioTnpOverMCtruth[(refEra,wp)])
                        zAxisName = f"TnP step {wpToTNP[wptmp]}: ratio TnP/MC {eraName}/{refEra}::0.985,1.015"    

                        outfolder_effRatioTnpOverMC_doubleRatio = f"{outfolder}/efficiency_ratioTNPoverMCtruth_eraOver{refEra}/stepTnP_{wpToTNP[wptmp]}" 
                        createPlotDirAndCopyPhp(outfolder_effRatioTnpOverMC_doubleRatio)

                        drawCorrelationPlot(ratioTnpOverMCtruth_ratioOverRefEra,
                                            xAxisName,
                                            yAxisName,
                                            zAxisName,
                                            ratioTnpOverMCtruth_ratioOverRefEra.GetName(), plotLabel="ForceTitle",
                                            outdir=outfolder_effRatioTnpOverMC_doubleRatio,
                                            draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                            invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)
                        
                
    if rfout.IsOpen():
        rfout.Close()
        print()
        print(f"All histograms saved in {rfout.GetName()}")
        print()
