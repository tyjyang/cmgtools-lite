#!/usr/bin/env python3

import re
import os, os.path
import logging
import argparse
import shutil

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

sys.path.append(os.getcwd())
from lumiStuff.runPerEra import lumiForEra_UL_new as lpe

if __name__ == "__main__":

    #eras = ["B", "C", "D", "E", "F", "BToF", "G", "H", "GToH"] # would have been better to have GtoH and BtoF but fine 
    #eras = ["F", "G"]
    eras = ["BToF", "GToH"]

    WorZ = "Z" # W or Z
    
    for era in eras:

        if era == "BToF":
            lumi = 19.514703
            totalLumi = 19.514703
            process = "Wmunu_plus_preVFP" if WorZ == "W" else "Zmumu_preVFP"
        elif era == "GToH":
            lumi = 16.810813
            totalLumi = 16.810813
            process = "Wmunu_plus_postVFP" if WorZ == "W" else "Zmumu_postVFP"
        else:
            lumi = lpe[era]
            if era in ["G", "H"]:
                totalLumi = 16.810813
                process = "Wmunu_plus_postVFP" if WorZ == "W" else "Zmumu_postVFP"
            else:
                totalLumi = 19.514703
                process = "Wmunu_plus_preVFP" if WorZ == "W" else "Zmumu_preVFP"

        lumiOverTotalLumi = lumi/totalLumi
                
        ### nominal case
        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /scratch/shared/originalNANO/ --sP 'bareMuon_pt_eta,bareMuon_pt_eta__accept,bareMuon_pt_eta__trackerOrGlobal,bareMuon_pt_eta__trackerOrGlobalAndStandalone,bareMuon_pt_eta__tracker,bareMuon_pt_eta__global,bareMuon_pt_eta__standalone,bareMuon_pt_eta__idip,bareMuon_pt_eta__trig,bareMuon_pt_eta__trigNoBit,bareMuon_pt_eta__idipANDisonotrig,bareMuon_pt_eta__iso,bareMuon_pt_eta__idipANDtrig,bareMuon_pt_eta__idipANDtrigANDiso,bareMuon_pt_eta__idipANDtrigNoBit,bareMuon_pt_eta__veto' --lumi-weight {lumi}  -W 'puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/{WorZ}_customNano_newPU_ptReco15andSameCharge/plus/{era}/ -p '{process}'   --nanoaod-tree  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"

        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_Z.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /scratch/shared/originalNANO/ --sP 'bareMuon_pt_eta,bareMuon_pt_eta__accept,bareMuon_pt_eta__trackerOrGlobal,bareMuon_pt_eta__trackerOrGlobalAndStandalone,bareMuon_pt_eta__tracker,bareMuon_pt_eta__global,bareMuon_pt_eta__standalone,bareMuon_pt_eta__idip,bareMuon_pt_eta__trig,bareMuon_pt_eta__trigNoBit,bareMuon_pt_eta__idipANDisonotrig,bareMuon_pt_eta__iso,bareMuon_pt_eta__idipANDtrig,bareMuon_pt_eta__idipANDtrigANDiso,bareMuon_pt_eta__idipANDtrigNoBit,bareMuon_pt_eta__veto' --lumi-weight {lumi}  -W 'puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/{WorZ}_perEra_noRecoAccept_finePt_newPUweightsFrom04June_1bareMuon_FG/plus/{era}/ -p '{process}'   --nanoaod-tree  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"

        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /scratch/shared/customNanoForEfficiency_addTrackInfo/ --sP 'bareMuon_pt_eta,bareMuon_pt_eta__accept,bareMuon_pt_eta__trackerOrGlobal,bareMuon_pt_eta__trackerOrGlobalAndStandalone,bareMuon_pt_eta__tracker,bareMuon_pt_eta__global,bareMuon_pt_eta__standalone,bareMuon_pt_eta__idip,bareMuon_pt_eta__trig,bareMuon_pt_eta__trigNoBit,bareMuon_pt_eta__idipANDisonotrig,bareMuon_pt_eta__iso,bareMuon_pt_eta__idipANDtrig,bareMuon_pt_eta__idipANDtrigANDiso,bareMuon_pt_eta__idipANDtrigNoBit,bareMuon_pt_eta__veto,bareMuon_pt_eta__generalTrack,bareMuon_pt_eta__generalTrackAnyCharge,bareMuon_pt_eta__basicTrackMatchedToTrackerOrGlobal,bareMuon_pt_eta__basicTrackMatchedToTrackerOrGlobalAnyCharge,bareMuon_pt_eta__basicTrackMatchedToStandalone,bareMuon_pt_eta__standaloneMatchedToGlobal,bareMuon_pt_eta__standaloneAndGlobal,bareMuon_pt_eta__basicTrackAllMatchedToStandalone,bareMuon_pt_eta__generalTrackAll' --lumi-weight {lumi}  -W 'puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/{WorZ}_customNanoWithTrackInfo_newPU_ptReco15andSameCharge_globalMuonWhenRelevant_genTagSelection/plus/{era}/ -p '{process}'   --nanoaod-tree  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' --rdf-alias 'bareMuonsOther: bareMuonsMinus:.*'  -v 3 --scale-factor-file './testMuonSF/scaleFactorProduct_28Oct2021_nodz_dxybs_genMatchDR01.root' -A genaccept tagSelection 'Sum(genTagMuons2p4)>0' "  # sf file not needed, but as we were changing histogram names the code was crashing because it still tries to load them
        #
        # Z
        command = f"python mcPlots.py -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /scratch/shared/customNanoForEfficiency_addTrackInfo/ --sP 'bareMuon_upar_pt,bareMuon_upar_pt__trackerOrGlobal,bareMuon_upar_pt__trackerOrGlobalAndStandalone,bareMuon_upar_pt__tracker,bareMuon_upar_pt__global,bareMuon_upar_pt__standalone,bareMuon_upar_pt__idip,bareMuon_upar_pt__trig,bareMuon_upar_pt__idipANDisonotrig,bareMuon_upar_pt__iso,bareMuon_upar_pt__idipANDtrig,bareMuon_upar_pt__idipANDtrigANDiso,bareMuon_upar_pt__veto,bareMuon_upar_pt__generalTrack,bareMuon_upar_pt__basicTrackMatchedToTrackerOrGlobal,bareMuon_upar_pt__basicTrackMatchedToStandalone,bareMuon_upar_pt__standaloneMatchedToGlobal,bareMuon_upar_pt__standaloneAndGlobal,bareMuon_upar_pt__basicTrackAllMatchedToStandalone,bareMuon_upar_pt__generalTrackAll' --lumi-weight {lumi}  -W 'puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/{WorZ}_customNanoWithTrackInfo_testRecoil/genTagSelection/plus/{era}/ -p '{process}'   --nanoaod-tree  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' --rdf-alias 'bareMuonsOther: bareMuonsMinus:.*'  -v 3 --scale-factor-file './testMuonSF/scaleFactorProduct_28Oct2021_nodz_dxybs_genMatchDR01.root' -A genaccept tagSelection 'Sum(genTagMuons2p4)>0' "  # sf file not needed, but as we were changing histogram names the code was crashing because it still tries to load them
        #
        # W
        #command = f"python mcPlots.py -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /scratch/shared/originalNANO_newWithAltPDF/ --sP 'bareMuon_upar_pt,bareMuon_upar_pt__trackerOrGlobal,bareMuon_upar_pt__tracker,bareMuon_upar_pt__global,bareMuon_upar_pt__idip,bareMuon_upar_pt__trig,bareMuon_upar_pt__idipANDisonotrig,bareMuon_upar_pt__iso,bareMuon_upar_pt__idipANDtrig,bareMuon_upar_pt__idipANDtrigANDiso' --lumi-weight {lumi}  -W 'puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/{WorZ}_testRecoil/noTagSelection/plus/{era}/ -p '{process}'   --nanoaod-tree  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff_notCustomNano.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' --rdf-alias 'bareMuonsOther: bareMuonsMinus:.*'  -v 3 --scale-factor-file './testMuonSF/scaleFactorProduct_28Oct2021_nodz_dxybs_genMatchDR01.root' --n-thread 128 "  # sf file not needed, but as we were changing histogram names the code was crashing because it still tries to load them

        
        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /scratch/shared/customNanoForEfficiency_addTrackInfo/ --sP 'bareMuon_genZ_eta,bareMuon_genZ_eta__trackerOrGlobal,bareMuon_genZ_eta__generalTrack,bareMuon_genZ_eta__trackerOrGlobalAndStandalone,bareMuon_genZ_eta__tracker,bareMuon_genZ_eta__global,bareMuon_genZ_eta__standalone,bareMuon_genZ_eta__idip,bareMuon_genZ_eta__idipANDisonotrig,bareMuon_genZ_eta__idipANDtrig,bareMuon_genZ_eta__idipANDtrigANDiso,bareMuon_genZ_eta__generalTrack,bareMuon_genZ_eta__basicTrackMatchedToTrackerOrGlobal,bareMuon_genZ_eta__basicTrackMatchedToTrackerOrGlobalAnyCharge,bareMuon_genZ_eta__basicTrackMatchedToStandalone,bareMuon_genZ_eta__standaloneMatchedToGlobal,bareMuon_genZ_eta__standaloneAndGlobal,bareMuon_genZ_eta__basicTrackAllMatchedToStandalone,bareMuon_genZ_eta__generalTrackAll' --lumi-weight {lumi}  -W 'puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/{WorZ}_customNanoWithTrackInfo_newPU_ptReco15andSameCharge_globalMuonWhenRelevant_EtaVsGenZ_extendGenEta/plus/{era}/ -p '{process}'   --nanoaod-tree  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3 --scale-factor-file './testMuonSF/scaleFactorProduct_09Oct2021_nodz_dxybs_genMatchDR01.root' "
        
        print(f"### ERA {era}")
        print(command)
        print()
