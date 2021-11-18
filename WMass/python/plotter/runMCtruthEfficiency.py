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
                
        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP 'bareMuon_pt_eta,bareMuon_pt_eta__accept,bareMuon_pt_eta__trackerOrGlobal,bareMuon_pt_eta__tracker,bareMuon_pt_eta__global,bareMuon_pt_eta__standalone,bareMuon_pt_eta__idip,bareMuon_pt_eta__trig,bareMuon_pt_eta__trigNoBit,bareMuon_pt_eta__idipANDisonotrig,bareMuon_pt_eta__iso,bareMuon_pt_eta__idipANDtrig,bareMuon_pt_eta__idipANDtrigANDiso,bareMuon_pt_eta__idipANDtrigNoBit'   -W '({lumiOverTotalLumi})*puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/W_perEra_noRecoAccept_finePt/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree   --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"

        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP 'bareMuon_pt_eta,bareMuon_pt_eta__accept,bareMuon_pt_eta__trackerOrGlobalAndStandalone,bareMuon_pt_eta__trackerOrGlobal,bareMuon_pt_eta__tracker,bareMuon_pt_eta__global,bareMuon_pt_eta__standalone,bareMuon_pt_eta__idip,bareMuon_pt_eta__trig,bareMuon_pt_eta__trigNoBit,bareMuon_pt_eta__idipANDisonotrig,bareMuon_pt_eta__iso,bareMuon_pt_eta__idipANDtrig,bareMuon_pt_eta__idipANDtrigANDiso,bareMuon_pt_eta__idipANDtrigNoBit'   -W '({lumiOverTotalLumi})' --pdir plots/testNanoAOD/MCtruthEfficiency/W_perEra_noRecoAccept_finePt_noPUweights/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree   --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"

        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP 'bareMuon_phi_eta,bareMuon_phi_eta__accept,bareMuon_phi_eta__trackerOrGlobalAndStandalone,bareMuon_phi_eta__trackerOrGlobal,bareMuon_phi_eta__tracker,bareMuon_phi_eta__global,bareMuon_phi_eta__standalone,bareMuon_phi_eta__idip,bareMuon_phi_eta__trig,bareMuon_phi_eta__trigNoBit,bareMuon_phi_eta__idipANDisonotrig,bareMuon_phi_eta__iso,bareMuon_phi_eta__idipANDtrig,bareMuon_phi_eta__idipANDtrigANDiso,bareMuon_phi_eta__idipANDtrigNoBit'   -W '({lumiOverTotalLumi})' --pdir plots/testNanoAOD/MCtruthEfficiency/W_perEra_noRecoAccept_EtaPhi_noPUweights/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree   --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"

        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP 'bareMuon_pt_eta_fullPtRange,bareMuon_pt_eta_fullPtRange__veto'   -W '({lumiOverTotalLumi})*puw_2016UL_era_OLD(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/testLeptonVetoNoPtCut/W_perEra_noRecoAccept_finePt/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree   --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"
        # --filter-proc-files '.*' '.*_[1-5].root'


        ### nominal case
        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP 'bareMuon_pt_eta,bareMuon_pt_eta__accept,bareMuon_pt_eta__trackerOrGlobal,bareMuon_pt_eta__trackerOrGlobalAndStandalone,bareMuon_pt_eta__tracker,bareMuon_pt_eta__global,bareMuon_pt_eta__standalone,bareMuon_pt_eta__idip,bareMuon_pt_eta__trig,bareMuon_pt_eta__trigNoBit,bareMuon_pt_eta__idipANDisonotrig,bareMuon_pt_eta__iso,bareMuon_pt_eta__idipANDtrig,bareMuon_pt_eta__idipANDtrigANDiso,bareMuon_pt_eta__idipANDtrigNoBit,bareMuon_pt_eta__veto' --lumi-weight {lumi}  -W 'puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/{WorZ}_customNano_newPU_ptReco15andSameCharge/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"

        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_Z.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP 'bareMuon_pt_eta,bareMuon_pt_eta__accept,bareMuon_pt_eta__trackerOrGlobal,bareMuon_pt_eta__trackerOrGlobalAndStandalone,bareMuon_pt_eta__tracker,bareMuon_pt_eta__global,bareMuon_pt_eta__standalone,bareMuon_pt_eta__idip,bareMuon_pt_eta__trig,bareMuon_pt_eta__trigNoBit,bareMuon_pt_eta__idipANDisonotrig,bareMuon_pt_eta__iso,bareMuon_pt_eta__idipANDtrig,bareMuon_pt_eta__idipANDtrigANDiso,bareMuon_pt_eta__idipANDtrigNoBit,bareMuon_pt_eta__veto' --lumi-weight {lumi}  -W 'puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/{WorZ}_perEra_noRecoAccept_finePt_newPUweightsFrom04June_1bareMuon_FG/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"

        command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/customNanoForEfficiency_addTrackInfo/ --sP 'bareMuon_pt_eta,bareMuon_pt_eta__accept,bareMuon_pt_eta__trackerOrGlobal,bareMuon_pt_eta__trackerOrGlobalAndStandalone,bareMuon_pt_eta__tracker,bareMuon_pt_eta__global,bareMuon_pt_eta__standalone,bareMuon_pt_eta__idip,bareMuon_pt_eta__trig,bareMuon_pt_eta__trigNoBit,bareMuon_pt_eta__idipANDisonotrig,bareMuon_pt_eta__iso,bareMuon_pt_eta__idipANDtrig,bareMuon_pt_eta__idipANDtrigANDiso,bareMuon_pt_eta__idipANDtrigNoBit,bareMuon_pt_eta__veto,bareMuon_pt_eta__generalTrack,bareMuon_pt_eta__generalTrackAnyCharge,bareMuon_pt_eta__basicTrackMatchedToTrackerOrGlobal,bareMuon_pt_eta__basicTrackMatchedToTrackerOrGlobalAnyCharge,bareMuon_pt_eta__basicTrackMatchedToStandalone,bareMuon_pt_eta__standaloneMatchedToGlobal,bareMuon_pt_eta__standaloneAndGlobal,bareMuon_pt_eta__basicTrackAllMatchedToStandalone,bareMuon_pt_eta__generalTrackAll' --lumi-weight {lumi}  -W 'puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/{WorZ}_customNanoWithTrackInfo_newPU_ptReco15andSameCharge_globalMuonWhenRelevant_genTagSelection/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' --rdf-alias 'bareMuonsOther: bareMuonsMinus:.*'  -v 3 --scale-factor-file './testMuonSF/scaleFactorProduct_09Oct2021_nodz_dxybs_genMatchDR01.root' -A genaccept tagSelection 'Sum(genTagMuons)>0' "  # sf file not needed, but as we were changing histogram names the code was crashing because it still tries to load them
        #

        
        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/customNanoForEfficiency_addTrackInfo/ --sP 'bareMuon_genZ_eta,bareMuon_genZ_eta__trackerOrGlobal,bareMuon_genZ_eta__generalTrack,bareMuon_genZ_eta__trackerOrGlobalAndStandalone,bareMuon_genZ_eta__tracker,bareMuon_genZ_eta__global,bareMuon_genZ_eta__standalone,bareMuon_genZ_eta__idip,bareMuon_genZ_eta__idipANDisonotrig,bareMuon_genZ_eta__idipANDtrig,bareMuon_genZ_eta__idipANDtrigANDiso,bareMuon_genZ_eta__generalTrack,bareMuon_genZ_eta__basicTrackMatchedToTrackerOrGlobal,bareMuon_genZ_eta__basicTrackMatchedToTrackerOrGlobalAnyCharge,bareMuon_genZ_eta__basicTrackMatchedToStandalone,bareMuon_genZ_eta__standaloneMatchedToGlobal,bareMuon_genZ_eta__standaloneAndGlobal,bareMuon_genZ_eta__basicTrackAllMatchedToStandalone,bareMuon_genZ_eta__generalTrackAll' --lumi-weight {lumi}  -W 'puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/{WorZ}_customNanoWithTrackInfo_newPU_ptReco15andSameCharge_globalMuonWhenRelevant_EtaVsGenZ_extendGenEta/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3 --scale-factor-file './testMuonSF/scaleFactorProduct_09Oct2021_nodz_dxybs_genMatchDR01.root' "
        
        print(f"### ERA {era}")
        print(command)
        print()
