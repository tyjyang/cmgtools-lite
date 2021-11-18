#!/bin/env python

# script to make PU weights from existing PU profiles in 2016 UL for Wmass analysis
# this is not fully equivalent to the weights obtained in postprocessing
# because of the slightly different treatment of large weights and profile renormalization
# the profiles being used are currently hardcoded inside this script
# for MC we use the sum of Pileup_nTrueInt_Wmunu_preVFP and Pileup_nTrueInt_Wmunu_postVFP
# because the PU profile is the same, so we get more statistics, especially for high PU bins

# it was verified that the PU weights modify the MC normalization only by less than a few 1e-7 (where "a few" depends on how high weights are cropped)




#  >>>>>>  NOTE <<<<<<<<<

# THIS SCRIPT IS OBSOLETE, I RECOMMEND USING pileupStuff/makePUweightPerEra.py

# >>>>>>>>>>>>>>>>>>>>>>>

import os, re, array, math
import time
import argparse
import copy

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

# from helicity analysis, used for comparisons
oldWeights = [0.3505407355600995, 0.8996968628890968, 1.100322319466069, 0.9562526765089195, 1.0366251229154624, 1.0713954619016586, 0.7593488199769544, 0.47490309461978414, 0.7059895997695581, 0.8447022252423783, 0.9169159386164522, 1.0248924033173097, 1.0848877947714115, 1.1350984224561655, 1.1589888429954602, 1.169048420382294, 1.1650383018054549, 1.1507200023444994, 1.1152571438041776, 1.0739529436969637, 1.0458014000030829, 1.032500407707141, 1.0391236062781293, 1.041283620738903, 1.0412963370894526, 1.0558823002770783, 1.073481674823461, 1.0887053272606795, 1.1041701696801014, 1.123218903738397, 1.1157169321377927, 1.1052520327174429, 1.0697489590429388, 1.0144652740600584, 0.9402657069968621, 0.857142825520793, 0.7527112615290031, 0.6420618248685722, 0.5324755829715156, 0.4306470627563325, 0.33289171600176093, 0.24686361729094983,0.17781595237914027, 0.12404411884835284, 0.08487088505600057, 0.056447805688061216, 0.03540829360547507, 0.022412461576677457, 0.013970541270658443, 0.008587896629717911, 0.004986410514292661, 0.00305102303701641, 0.001832072556146534, 0.0011570757619737708, 0.0008992999249003301, 0.0008241241729452477, 0.0008825716073180279, 0.001187003960081393, 0.0016454104270429153, 0.0022514113879764414, 0.003683196037880878, 0.005456695951503178, 0.006165248770884191, 0.007552675218762607, 0.008525338219226993, 0.008654690499815343, 0.006289068906974821, 0.00652551838513972, 0.005139581024893171, 0.005115751962934923, 0.004182527768384693, 0.004317593022028565, 0.0035749335962533355, 0.003773660372937113, 0.002618732319396435, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--doDataVSdata", action="store_true", help="Compare PU profiles in data (pre and post VFP (overrides isPreVFP)")
    parser.add_argument("--doInclusiveData", action="store_true", help="Use full 2016 dataset for comparison (overrides isPreVFP, but overridden by doDataVSdata)")
    parser.add_argument("--isPreVFP", action="store_true", help="Make PU profiles using preVFP MC")
    parser.add_argument("--normalizeIntegralUpToBin", type=int, default=100, help="Use integral up to this bin to normalize PU profiles")
    parser.add_argument("--maxWeightWarning", type=float, default=5.0, help="At the end issue a warning if weight > this value")
    parser.add_argument("--cropHighWeight", type=float, default=5.0, help="If larger than 0, crop values larger than maxWeightWarning, setting them to this value (i.e. use 1 for no reweighting, or same as maxWeightWarning)")
    args = parser.parse_args()

    newWeights = []

    normalizeIntegralUpToBin = args.normalizeIntegralUpToBin # this would usually correspond to the max bin where there is enough stat in data and MC, above which weird weights can be observed. Those large weights are set to cropHighWeight, and we don't care about renormalizing the distributions again to correct the weights
    # It might make sense to run the script a first time with the full integral to get a clue about the ideal value, and then run again with customized range 
    doDataVSdata = args.doDataVSdata # if True it overrides isPreVFP below
    doInclusiveData = args.doInclusiveData # overrides isPreVFP, but overridden by doDataVSdata
    isPreVFP = args.isPreVFP
    maxWeightWarning = args.maxWeightWarning  # at the end issue a warning if weight > this value
    cropHighWeight = args.cropHighWeight      # crop values larger than maxWeightWarning, setting them to this value (i.e. use 1 for no reweighting, or same as maxWeightWarning)
    #mcfile   = "mcPU_2016_W94XnoCuts_onfirst10extParts_nTrueInt.root" # no genWeight
    #mcname   = "h_nTrueInt"
    #
    if doDataVSdata:
        mcfile   = "pileupStuff/MyDataPileupHistogram_2016Legacy_FpostHIPMandGH.root"
        mcname   = "pileup"
        mcLabel  = "data FtoH"
        #
        datafile = "pileupStuff/MyDataPileupHistogram_2016Legacy_upTo2016FwithHIPM.root"
        dataname = "pileup"
        dataLabel = "data BtoF"
        #
        outdir   = "plots/testNanoAOD/PU_weights/2016_dataVSdata/"
        lumi     = "19.3" 
        ratioLabel = "pre/post VFP::0.0,1.5"
    elif doInclusiveData:
        #mcfile   = "pileupStuff/MyMCPileupHistogram_2016Legacy_Wplus_preVFP.root"  # made with genWeights
        mcfile   = "pileupStuff/MyMCPileupHistogram_2016Legacy_noGenWeights_preAndPostVFP.root"    
        mcname   = "Pileup_nTrueInt_Wmunu_preVFP,Pileup_nTrueInt_Wmunu_postVFP" # or Pileup_nTrueInt_Wmunu_preVFP or Pileup_nTrueInt_Wmunu_postVFP, the profile is the same
        mcLabel  = "W MC"    
        #
        datafile = "pileupStuff/MyDataPileupHistogram_2016Legacy_all2016.root"
        dataname = "pileup"
        dataLabel = "data (all)"
        #
        outdir   = "plots/testNanoAOD/PU_weights/2016_allData/"
        lumi     = "35.9" 
        ratioLabel = "data/MC::0.0,1.5"
    elif isPreVFP:
        #mcfile   = "pileupStuff/MyMCPileupHistogram_2016Legacy_Wplus_preVFP.root"  # with genWeight
        mcfile   = "pileupStuff/MyMCPileupHistogram_2016Legacy_noGenWeights_preAndPostVFP.root"    
        mcname   = "Pileup_nTrueInt_Wmunu_preVFP,Pileup_nTrueInt_Wmunu_postVFP"
        #mcfile   = "/afs/cern.ch/user/m/mciprian/www/wmass/13TeV/testNanoAOD/PU_weights/PUprofilesMC_allEventsWmunu/plots_test.root"    
        #mcname   = "Pileup_nTrueInt_Wmunu_preVFP"
        mcLabel  = "W MC"    
        #
        datafile = "pileupStuff/MyDataPileupHistogram_2016Legacy_upTo2016FwithHIPM.root"
        dataname = "pileup"
        dataLabel = "data"
        #
        outdir   = "plots/testNanoAOD/PU_weights/2016_preVFP/"
        lumi     = "19.3" 
        ratioLabel = "data/MC::0.0,1.5"
    else:
        #mcfile   = "pileupStuff/MyMCPileupHistogram_2016Legacy_Z_postVFP.root" # with genWeights
        mcfile   = "pileupStuff/MyMCPileupHistogram_2016Legacy_noGenWeights_preAndPostVFP.root"    
        mcname   = "Pileup_nTrueInt_Wmunu_preVFP,Pileup_nTrueInt_Wmunu_postVFP"
        mcLabel  = "W MC"
        #
        datafile = "pileupStuff/MyDataPileupHistogram_2016Legacy_FpostHIPMandGH.root"
        dataname = "pileup"
        dataLabel = "data"
        #
        outdir   = "plots/testNanoAOD/PU_weights/2016_postVFP/"
        lumi     = "16.6" 
        ratioLabel = "data/MC::0.0,1.5"

    tf = ROOT.TFile.Open(datafile)
    datahist = tf.Get(dataname)
    datahist.SetDirectory(0)
    tf.Close()

    mchist = None
    tf = ROOT.TFile.Open(mcfile)
    for i,mcn in enumerate(mcname.split(',')):
        tmp = tf.Get(mcn)
        if i == 0:
            mchist = copy.deepcopy(tmp.Clone("mchist"))
            mchist.SetDirectory(0)
        else:
            mchist.Add(tmp)
    tf.Close()

    if not datahist or not mchist:
        print("Error getting histograms")
        quit()

    if datahist.GetNbinsX() != mchist.GetNbinsX():
        print("Warning: histograms have different number of bins")
        quit()
    nbins = datahist.GetNbinsX()
    if datahist.GetXaxis().GetBinLowEdge(1) != mchist.GetXaxis().GetBinLowEdge(1):
        print("Warning: histograms have different low edge")
        quit()
    if datahist.GetXaxis().GetBinLowEdge(1+nbins) != mchist.GetXaxis().GetBinLowEdge(1+nbins):
        print("Warning: histograms have different up edge")
        quit()

    #if not ROOT.TH1.CheckConsistency(datahist,mchist):
    #    print("Warning: histograms not compatible"
    #    quit()

    integralDataFull = datahist.Integral(0,1+datahist.GetNbinsX()) 
    integralMCFull = mchist.Integral(0,1+mchist.GetNbinsX()) 
    integralDataRange = datahist.Integral(0,normalizeIntegralUpToBin) 
    integralMCRange = mchist.Integral(0,normalizeIntegralUpToBin) 
    
    mchist.Scale(datahist.Integral(0,normalizeIntegralUpToBin)/mchist.Integral(0,normalizeIntegralUpToBin))
    integralMCafterNormData = mchist.Integral(0,1+mchist.GetNbinsX())
    ratio = datahist.Clone("puWeight_2016_UL")
    oldratio = datahist.Clone("puWeight_2016_ReReco")
    ratio.Reset("ICESM")
    oldratio.Reset("ICESM")

    for i in range(1,1+nbins):
        puratio = 1.0
        if mchist.GetBinContent(i) == 0:
            ratio.SetBinContent(i,puratio)
        else:
            puratio = datahist.GetBinContent(i)/mchist.GetBinContent(i)
            ratio.SetBinContent(i,puratio)
        newWeights.append(puratio)
        oldratio.SetBinContent(i,oldWeights[i-1])
    #ratio.Divide(datahist,mchist)

    #drawTH1(ratio,"number of true interactions","PU weight ()",outdir,ratio.GetName(),"","900,800",None,"")

    hists = [datahist,mchist]
    legEntries = [dataLabel,mcLabel]
    drawNTH1(hists,legEntries,"number of true interactions","Events","comparePU_data_MC",outdir,
             labelRatioTmp=ratioLabel,legendCoords="0.6,0.9,0.76,0.9",canvasSize="900,900",
             lumi=lumi)

    hists = [ratio,oldratio]
    legEntries = ["new","SMP-18-012"]
    drawNTH1(hists,legEntries,"number of true interactions","PU weight","comparePU_UL_ReReco",outdir,
             labelRatioTmp="new/old::0.5,1.5",legendCoords="0.5,0.95,0.76,0.9",canvasSize="900,900",
             lumi=lumi)

    mchistWeight = copy.deepcopy(mchist.Clone("mchistWeight"))
    mchistWeightCrop = copy.deepcopy(mchist.Clone("mchistWeightCrop"))
    for i in range(1,1+nbins):
        mchistWeight.SetBinContent(i,newWeights[i-1]*mchistWeight.GetBinContent(i))
    integralMCafterWeight = mchistWeight.Integral(0,1+mchistWeight.GetNbinsX())
    hists = [datahist,mchistWeight]
    legEntries = [dataLabel,mcLabel+" (weight)"]
    drawNTH1(hists,legEntries,"number of true interactions","Events","comparePU_data_MCweight",outdir,
             labelRatioTmp=ratioLabel,legendCoords="0.5,0.9,0.76,0.9",canvasSize="900,900",
             lumi=lumi)

    print("")
    print("-"*30)
    for i in range(0,nbins):
        if newWeights[i] > maxWeightWarning:
            print("nPU = %d -> weight = %.3f %s" % (i+1,newWeights[i], f"cropping to {cropHighWeight}" if cropHighWeight > 0 else ""))
            if cropHighWeight > 0:
                newWeights[i] = cropHighWeight
        mchistWeightCrop.SetBinContent(i+1,newWeights[i]*mchistWeightCrop.GetBinContent(i+1))
    integralMCafterWeightCrop = mchistWeightCrop.Integral(0,1+mchistWeightCrop.GetNbinsX())

    print("-"*30)
    print("")

    print("")
    print("Printing new PU weights %s" % f"(cropped to {cropHighWeight} if > {maxWeightWarning})" if cropHighWeight > 0 else "")
    print("-"*30)
    print(newWeights)
    print("-"*30)
    print("Summary of integrals in full range")
    print(f"Data initial              : {integralDataFull}")
    print(f"MC initial                : {integralMCFull}")
    print(f"MC after norm   to data   : {integralMCafterNormData}")
    print(f"MC after PU weights       : {integralMCafterWeight}")
    print(f"MC after PU weights (crop): {integralMCafterWeightCrop}")
    
