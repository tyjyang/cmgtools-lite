#!/usr/bin/env python3

from shutil import copyfile
import re, sys, os, os.path, subprocess, json, ROOT
import numpy as np
import argparse

# processes to group
def getProcessGroup(proc, era):
    # check whether Wmunu is split by charge or includes both samples, same for Wtaunu
    if any(x in proc for x in ["Wmunu", "Wtaunu"]) and all(c not in proc for c in ["plus", "minus"]):
        if era == "all":
            ret = f"--pg '{proc} := {proc}_plus_preVFP,{proc}_plus_postVFP,{proc}_minus_preVFP,{proc}_minus_postVFP'"
        else:
            ret = f"--pg '{proc} := {proc}_plus_{era},{proc}_minus_{era}'"
    else:
        if era == "all":
            ret = f"--pg '{proc} := {proc}_preVFP,{proc}_postVFP'"
        else:
            ret = f"--pg '{proc} := {proc}_{era}'"
    return ret


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dry-run', dest="dryRun", action='store_true', help='Print, but do not execute')
    parser.add_argument('-e', '--era', type=str, default="all", choices=["all","preVFP","postVFP"], help='Era')
    args = parser.parse_args()

    era = args.era

    ######################################################################
    ## Configurations below
    ######################################################################
    lumi = 16.8 if era == "postVFP" else 19.5 if era == "preVFP" else 36.3

    # cfg files
    cfgFolder = "w-mass-13TeV/testingNano/cfg/"
    mca    = cfgFolder + "mca-wmass.txt"
    cut    = cfgFolder + "test/cuts_fakerate.txt"
    plot   = cfgFolder + "plots_fakerate.txt"
    define = cfgFolder + "test/rdfDefine_fakerate.txt"

    # additional cuts
    addcut = ""
    #addcut = "-A awayjet pfRelIso04 'Muon_pfRelIso04_all[goodMuons][0] < 0.15' "

    # additional defines and aliases
    otherDefines = " --rdf-alias 'goodMuonsCharge: goodMuons:.*'"

    # input samples
    samples = "/data/shared/originalNANO/"

    # output
    plotdir = f"plots/testNanoAOD/WmassPlots/fakeRateRegion_TEST_{era}/"

    # histograms to make (use .* to activate all those in plot file)
    hists = "muon_pt,muon_eta_fine,mt_MET,MET_pt,deltaPhi_MuMet,nJetClean,leadJetClean_pt,leadJetClean_eta,deta_jetmu,dR_jetmu"

    # processes and related options (also for ratio plots when customizing names)
    processes = "data,Zmumu,Ztautau,Wmunu,Wtaunu,Top,Diboson"
    procGroups = " ".join([getProcessGroup(p,era) for p in processes.split(',')])
    procOptions = f"-p '{processes}' {procGroups}"

    # ratio (settings, while customization of processes should go in procOptions)
    ratio = "--showRatio --maxRatioRange 0.5 1.5 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel 'Data/MC'"

    # event weight (global one for MC)
    weight = f"puw_2016UL_era(Pileup_nTrueInt,eraVFP)*_get_fullMuonSF(Muon_pt[goodMuonsCharge][0],Muon_eta[goodMuonsCharge][0],Muon_charge[goodMuonsCharge][0],-1,-1,eraVFP,Muon_pfRelIso04_all[goodMuons][0]<0.15)*_get_MuonPrefiringSF(Muon_eta,Muon_pt,Muon_looseId,eraVFP)"

    # gen weight customization
    genweight = "--max-genWeight-procs 'W|Z' '50118.72' --clip-genWeight-toMax"

    # whatever with legend
    legOptions = "--legendFontSize 0.042 --allProcInLegend --n-column-legend 2 --setLegendCoordinates 0.2,0.76,0.9,0.92"

    # general options not in a specific group
    general = "-f  -j 8  --nanoaod-tree -v 3 --skipPlot"

    ######################################################################
    ## Finally the command
    ######################################################################
    
    command = f"python mcPlots.py -l {lumi} {mca} {cut} {plot} --noCms -P {samples} --sP '{hists}'   -W '{weight}' --pdir {plotdir} {procOptions} {legOptions}  {genweight}  --rdf-define-file {define} {otherDefines} {addcut} {ratio} {general}"

    if args.dryRun:
        print(command)
    else:
        os.system(command)
