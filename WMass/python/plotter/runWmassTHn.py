#!/usr/bin/env python3

# charge plus
# python runFakeRate.py -o plots/testNanoAOD/WmassPlots/Wnnpdf30_noPDFonZ/ -e postVFP --variables ".*" --plot-file "plots_fakerate_systTH3_nnpdf30.txt" --options " --skipPlot " -c plus -s
# charge minus
# python runFakeRate.py -o plots/testNanoAOD/WmassPlots/Wnnpdf30_noPDFonZ/ -e postVFP --variables ".*" --plot-file "plots_fakerate_systTH3_nnpdf30.txt" --options " --skipPlot " -c minus -s

from shutil import copyfile
import re, sys, os, os.path, subprocess, json, ROOT
import numpy as np
import argparse

from lumiStuff.runPerEra import lumiForEra_UL_new as lpe
from wmassUtilities.theory import pdfMap

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import createPlotDirAndCopyPhp

# processes to group
def getProcessGroup(proc, era, subEra=None):
    # check whether Wmunu is split by charge or includes both samples, same for Wtaunu
    if "data" in proc and subEra != None:
        era = subEra
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


def main(args):

    era = args.era
    subEra = args.subEra
    folderEra = era if subEra == None else subEra
    postfix = args.postfix
    if len(args.postfix) and not postfix.startswith('_'):
        postfix = "_" + postfix
    outdir = args.outdir
        
    ######################################################################
    ## Configurations below
    ######################################################################
    lumi = 16.8 if era == "postVFP" else 19.5 if era == "preVFP" else 36.3
    if subEra != None:
        lumi = lpe[subEra]
        #postfix += f"_2016{subEra}"
        
    # cfg files
    cfgFolder = args.cfgFolder
    mca    = cfgFolder + "mca-wmass.txt"
    cut    = cfgFolder + "test/cuts_fakerate.txt"
    plot   = cfgFolder + args.plotFile
    define = cfgFolder + "test/rdfDefine_fakerate.txt"

    # input samples
    samples = "/scratch/shared/originalNANO/"

    # output
    plotdir = f"{outdir}/fakeRateRegion{postfix}/{folderEra}/allTHn/"

    # histograms to make (use .* to activate all those in plot file)
    hists = args.variables
    #hists = "muon_pt_eta_isoMtRegions"
    
    # processes and related options (also for ratio plots when customizing names)
    processes = "data,Zmumu,Ztautau,Top,Diboson,Wmunu_plus,Wtaunu_plus,Wmunu_minus,Wtaunu_minus"
    if args.useVptWeight:
        originalprocs = [str(x) for x in processes.split(',')]
        tmpprocs = []
        for x in originalprocs:
            if any(x.startswith(k) for k in ["W", "Z"]):
                tmpprocs.append(x + "_vpt")
            else:
                tmpprocs.append(x)
        processes = ",".join(tmpprocs)

    print('-'*30)
    processesList = processes.split(',')
    print(f"Using these processes ({len(processesList)} in total)")
    for i,p in enumerate(processesList):
        print(f"{i}) {p}")
    print('-'*30)
        
    procGroups = " ".join([getProcessGroup(p, era, subEra) for p in processes.split(',')])
    procOptions = f"-p '{processes}' {procGroups}"

    # ratio (settings, while customization of processes should go in procOptions)
    ratio = "--showRatio --maxRatioRange 0.5 1.5 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel 'Data/MC'"

    # event weight (global one for MC)
    weight = ""
    if subEra != None:
        #
        # no reco*tracking here for now!
        weight = f"puw_2016UL_era(Pileup_nTrueInt,{subEra})*_get_fullMuonSF_perDataEra(Muon_pt[goodMuons][0],Muon_eta[goodMuons][0],Muon_charge[goodMuons][0],-1,-1,{subEra},Muon_pfRelIso04_all[goodMuons][0]<0.15)*_get_newMuonPrefiringSF(Muon_eta,Muon_pt,Muon_phi,Muon_looseId,{subEra})"
        ####
    else:
        weight = f"puw_2016UL_era(Pileup_nTrueInt,eraVFP)*_get_fullMuonSF(Muon_pt[goodMuons][0],Muon_eta[goodMuons][0],Muon_charge[goodMuons][0],-1,-1,eraVFP,Muon_pfRelIso04_all[goodMuons][0]<0.15)*_get_newMuonPrefiringSF(Muon_eta,Muon_pt,Muon_phi,Muon_looseId,eraVFP)*_get_tnpRecoSF(Muon_pt[goodMuons][0],Muon_eta[goodMuons][0],Muon_charge[goodMuons][0],-1,-1,eraVFP,0,reco)*_get_tnpTrackingSF(Muon_pt[goodMuons][0],Muon_eta[goodMuons][0],Muon_charge[goodMuons][0],-1,-1,eraVFP)"


    pdfInfo = pdfMap[args.pdf]
    
    pdfMod = ""    
    processesWithPDFs = "W.*" if pdfInfo["onlyW"] else "W.*|Z.*"
    if args.pdf != "nnpdf31":
        pdfMod = f"--rdf-define 'nominalWeight : {pdfInfo['weight']}[0] : {processesWithPDFs} : 1.0'"
        weight = f"{weight}*nominalWeight"

    # SF file
    sfOpt = f" --scale-factor-file '{args.scaleFactorFile}'"
    # the following is needed for SF made before June 2021
    if args.oldSFname:
        sfOpt += " --old-sf-name"
    
    # whatever with legend
    legOptions = "--legendFontSize 0.042 --allProcInLegend --n-column-legend 2 --setLegendCoordinates 0.2,0.76,0.9,0.92"

    # general options not in a specific group
    general = "--nanoaod-tree -v 3"

    ######################################################################
    ## Finally the command
    ######################################################################
    
    command = f"python mcPlots.py -l {lumi} {mca} {cut} {plot} --noCms -P {samples} --sP '{hists}' -W '{weight}' --pdir {plotdir} {procOptions} {legOptions} --rdf-define-file {define} {ratio} {general} {sfOpt} {pdfMod} {args.options}"
    if args.subEra:
        command += f" --lumi-weight {lumi} "

    if args.checkTime:
        command = "/usr/bin/time -v " + command
        
    if args.dryRun:
        print(command)
    else:
        print("Going to run the following mcPlots.py command")
        print()
        print(command)
        print("-"*30)
        os.system(command)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dry-run', dest="dryRun", action='store_true', help='Print, but do not execute')
    parser.add_argument(      '--vpt-weight', dest="useVptWeight", action='store_true', help='Apply Wpt reweighting to central value using SCETlib correction')
    parser.add_argument('-e', '--era',     type=str, default="all", choices=["all","preVFP","postVFP"], help='Era')
    parser.add_argument('--sub-era', dest="subEra",  type=str, default=None, choices=["B","C","D","E","F","G","H"], help='Optional, sub era (needs dedicated SF and appropriate root file as input)')
    parser.add_argument('-p', '--postfix', type=str, default="", help='Postfix to output folder name')
    parser.add_argument('-o', '--outdir', required=True, type=str, default="", help='Output folder')
    parser.add_argument(      '--options', type=str, default="", help='Other options to pass to command, if not already present')
    parser.add_argument(      '--variables', type=str, default=".*", help='Histograms to make from the configuration txt file (.* for all)')
    parser.add_argument(      '--cfg-folder', dest="cfgFolder", type=str, default="w-mass-13TeV/testingNano/cfg/", help='Folder where cfg files are taken. Can leave this default')
    parser.add_argument(      '--plot-file', dest="plotFile", type=str, default="plots_fakerate.txt", help='File with histogram definition (inside cfgFolder)')
    parser.add_argument("--scale-factor-file", dest="scaleFactorFile", type=str, default="./testMuonSF/scaleFactorProduct_28Oct2021_nodz_dxybs_genMatchDR01.root", help="File to be used for scale factors")
    parser.add_argument("--old-sf-name", dest="oldSFname", action="store_true", help="To use old SF in file, whose names did not have nominal or dataAltSig (but note that eff.syst requires the new version)");
    parser.add_argument("-t", "--time", dest="checkTime", action="store_true", help="Check time of execution using /usr/bin/time -v");
    parser.add_argument('--pdf', default="nnpdf31", choices=pdfMap.keys(), help='PDF set to use')
    args = parser.parse_args()
    
    createPlotDirAndCopyPhp(args.outdir)

    print("="*40)
    print("="*40)
    print("="*40)
    print()
    #print("WARNING: using old PU weights (before June 2021, older PU json was used)")
    #print("Please update this script when you have the new ones")
    print()
    print("="*40)
    print("="*40)
    print("="*40)
        
    if args.useVptWeight:
        args.postfix += "_vptWeight"
        
    args.postfix += "_systTHn"
    # added back to cut file
    #regionCut = "transverseMass >= 40.0 || Sum(goodCleanJets)>=1"
    #args.options += f" -A trigMatch regionCut '{regionCut}' "
    print("")
    main(args)
    print('-'*30)
    print("")
