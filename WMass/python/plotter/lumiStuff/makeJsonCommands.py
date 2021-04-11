#!/usr/bin/env python

from runPerEra import runsForEra as rfe

## for UL
#eras = list(rfe.keys())
#goldenJson = "pileupStuff/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
#inputDir = "/data/shared/originalNANO/NanoV8Data/"
#tag = ""

## for ReReco
eras = ["B", "C", "D", "E", "F", "G", "H"] 
goldenJson = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt"
inputDir = "/data/shared/NANO02Apr2020_fromReReco/"
tag = "_ReReco"

for era in eras:
    command = f"python lumiStuff/makeJsonFromData.py -i {inputDir}/Run2016{era}/ -o lumiStuff/ -n NanoAOD_json_Run2016{era}_GoldenJsonFilter{tag}.txt -j {goldenJson}"
    print(command)
    print()
