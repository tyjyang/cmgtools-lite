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

    eras = ["B", "C", "D", "E", "F", "BToF", "G", "H", "GToH"]
    #eras = ["BToF", "GToH"]
    
    for era in eras:

        if era == "BToF":
            lumi = 19.514703
            totalLumi = 19.514703
            process = "Wmunu_plus_preVFP"
        elif era == "GToH":
            lumi = 16.810813
            totalLumi = 16.810813
            process = "Wmunu_plus_postVFP"
        else:
            lumi = lpe[era]
            if era in ["G", "H"]:
                totalLumi = 16.810813
                process = "Wmunu_plus_postVFP"
            else:
                totalLumi = 19.514703
                process = "Wmunu_plus_preVFP"

        lumiOverTotalLumi = lumi/totalLumi
                
        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP 'bareMuon_pt_eta,bareMuon_pt_eta__accept,bareMuon_pt_eta__trackerOrGlobal,bareMuon_pt_eta__tracker,bareMuon_pt_eta__global,bareMuon_pt_eta__standalone,bareMuon_pt_eta__idip,bareMuon_pt_eta__trig,bareMuon_pt_eta__trigNoBit,bareMuon_pt_eta__idipANDisonotrig,bareMuon_pt_eta__iso,bareMuon_pt_eta__idipANDtrig,bareMuon_pt_eta__idipANDtrigANDiso,bareMuon_pt_eta__idipANDtrigNoBit'   -W '({lumiOverTotalLumi})*puw_2016UL_era(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/W_perEra_noRecoAccept_finePt/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree   --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"

        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP 'bareMuon_pt_eta,bareMuon_pt_eta__accept,bareMuon_pt_eta__trackerOrGlobalAndStandalone,bareMuon_pt_eta__trackerOrGlobal,bareMuon_pt_eta__tracker,bareMuon_pt_eta__global,bareMuon_pt_eta__standalone,bareMuon_pt_eta__idip,bareMuon_pt_eta__trig,bareMuon_pt_eta__trigNoBit,bareMuon_pt_eta__idipANDisonotrig,bareMuon_pt_eta__iso,bareMuon_pt_eta__idipANDtrig,bareMuon_pt_eta__idipANDtrigANDiso,bareMuon_pt_eta__idipANDtrigNoBit'   -W '({lumiOverTotalLumi})' --pdir plots/testNanoAOD/MCtruthEfficiency/W_perEra_noRecoAccept_finePt_noPUweights/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree   --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"

        #command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP 'bareMuon_phi_eta,bareMuon_phi_eta__accept,bareMuon_phi_eta__trackerOrGlobalAndStandalone,bareMuon_phi_eta__trackerOrGlobal,bareMuon_phi_eta__tracker,bareMuon_phi_eta__global,bareMuon_phi_eta__standalone,bareMuon_phi_eta__idip,bareMuon_phi_eta__trig,bareMuon_phi_eta__trigNoBit,bareMuon_phi_eta__idipANDisonotrig,bareMuon_phi_eta__iso,bareMuon_phi_eta__idipANDtrig,bareMuon_phi_eta__idipANDtrigANDiso,bareMuon_phi_eta__idipANDtrigNoBit'   -W '({lumiOverTotalLumi})' --pdir plots/testNanoAOD/MCtruthEfficiency/W_perEra_noRecoAccept_EtaPhi_noPUweights/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree   --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"

        command = f"python mcPlots.py -f -l {lumi} w-mass-13TeV/testingNano/cfg/test/mca_MCtruthEff.txt w-mass-13TeV/testingNano/cfg/test/cuts_MCtruthEff_W.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP 'bareMuon_pt_eta_fullPtRange,bareMuon_pt_eta_fullPtRange__veto'   -W '({lumiOverTotalLumi})*puw_2016UL_era_OLD(Pileup_nTrueInt,{era})' --pdir plots/testNanoAOD/MCtruthEfficiency/testLeptonVetoNoPtCut/W_perEra_noRecoAccept_finePt/plus/{era}/ -p '{process}'  -j 8  --nanoaod-tree   --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_MCtruthEff.txt --rdf-alias 'bareMuonsCharge: bareMuonsPlus:.*' -v 3"
        # --filter-proc-files '.*' '.*_[1-5].root'

        print(f"### ERA {era}")
        print(command)
        print()
