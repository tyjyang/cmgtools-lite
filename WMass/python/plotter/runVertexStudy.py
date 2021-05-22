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

if __name__ == "__main__":

    cuts = ["alwaystrue", "onemuon", "muonID", "trigAndMatch", "pfRelIso04"]
    weight = "1.0"
    #weight = "puw_2016UL_era(Pileup_nTrueInt,eraVFP)*_get_MuonPrefiringSF(Muon_eta,Muon_pt,Muon_looseId,eraVFP)"
    
    for cut in cuts:

        command = f"python mcPlots.py -f -l 1.0 w-mass-13TeV/testingNano/cfg/mca-vertexStudy.txt w-mass-13TeV/testingNano/cfg/test/cuts_vertexStudy.txt w-mass-13TeV/testingNano/cfg/plots_vertexStudy.txt --noCms -P /data/shared/originalNANO/ --sP 'dzGenRecoVtx_wpt,dzGenRecoVtx_mupt,mueta_dzGenRecoVtx_wpt'  -W '{weight}' --pdir plots/testNanoAOD/vertexStudy/Wboson_noPUorPrefire/{cut}/ -p 'Wmunu_plus_preVFP,Wmunu_plus_postVFP,Wmunu_minus_preVFP,Wmunu_minus_postVFP,Zmumu_preVFP,Zmumu_postVFP' --legendFontSize 0.045 --allProcInLegend --n-column-legend 2 --setLegendCoordinates 0.2,0.76,0.9,0.92  -j 8  --nanoaod-tree   --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_vertexStudy.txt --rdf-alias 'goodMuonsCharge: goodMuonsPlus:.*' --rdf-alias 'goodMuonsOther: goodMuonsMinus:.*' -v 3  --plotmode nostack -U {cut} --skipPlot "
        # --skipPlot

        print(f"### CUT {cut}")
        print(command)
        print()
