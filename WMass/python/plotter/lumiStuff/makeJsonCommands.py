#!/usr/bin/env python

from runPerEra import runsForEra as rfe

eras = list(rfe.keys())
goldenJson = "pileupStuff/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
inputDir = "/data/shared/originalNANO/NanoV8Data/"

for era in eras:
    command = f"python lumiStuff/makeJsonFromData.py -i {inputDir}/Run2016{era}/ -o lumiStuff/ -n NanoAOD_json_Run2016{era}_GoldenJsonFilter.txt -j {goldenJson}"
    print(command)
    print()
