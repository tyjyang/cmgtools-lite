#!/usr/bin/env python3

import os, re, array, math
import argparse
## safe batch mode
import sys
import json

allEras = {"preVFP"  : ["B", "C", "D", "E", "F"],
           "postVFP" : ["F_postVFP", "G", "H"]
           }

inputJsonBase = "lumiStuff/NanoAOD_json_Run2016XXX_GoldenJsonFilter.txt" # will change XXX with actual eras

for eraVFP in allEras.keys():

    final = None
    for i,era in enumerate(allEras[eraVFP]):
        jsonFile = inputJsonBase.replace("XXX", era)
        if i:
            info = json.load(open(jsonFile))
            final.update(info)
        else:
            final = json.load(open(jsonFile))
    output = inputJsonBase.replace("XXX", eraVFP)
    with open(output, 'w') as outfile:
        # can add indent=0 for json.dumps, but it adds newline after each LS rather than after blocks
        outfile.write(json.dumps(final))
        print(f"Json file for {allEras[eraVFP]} written in {output}")
