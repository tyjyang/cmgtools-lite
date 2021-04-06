#!/usr/bin/env python  

# script to check integrated luminosity per era, summing the values from brilcalc for single runs or even by LS

from shutil import copyfile
import re, sys, os, os.path, subprocess, json, ROOT
import numpy as np
import argparse

from runPerEra import runsForEra as rfe
from makeJsonFromData import compileMacro

sys.path.append(os.getcwd())
from mcAnalysis import initializeJson

parser = argparse.ArgumentParser()
parser.add_argument('--by-lumi', dest='byLumi',  action='store_true', help='Check sum from single lumisections (default is by run). Requires a json')
parser.add_argument('-j','--json',  default='', type=str, help='Json file to be used to filter runs and lumis in dataset')
 
args = parser.parse_args()

inputfile = "lumiStuff/runlumi.csv" # copied from /afs/cern.ch/user/b/bendavid/cmspublic/runlumi.csv
if args.byLumi:
    inputfile = "lumiStuff/bylsoutput.csv" # /afs/cern.ch/user/b/bendavid/cmspublic/bylsoutput.csv 
    if args.json == "":
        print("Error: option --by-lumi requires using a json with -j/--json <jsonfile>")
        quit()
        
filterWithJson = False
if args.json != "":
    #print("ROOT libraries")                                                                                       
    #print(ROOT.gSystem.GetLibraries())                                                                            
    if "/jsonManager_cc.so" not in ROOT.gSystem.GetLibraries():
        compileMacro("../ccFiles/jsonManager.cc")
    if hasattr(ROOT, "jsonMap_all"):
        initializeJson(ROOT.jsonMap_all, args.json)
        print("Initialized json file")
        filterWithJson = True

    
eras = list(rfe.keys())
lumiPerEra = {era : 0.0 for era in eras}

with open(inputfile) as f:
    for line in f:
        if line.startswith('#'):
            continue
        tokens = line.strip().split(',')
        if args.byLumi:
            tokens[0] = tokens[0].split(':')[0]
            tokens[1] = tokens[1].split(':')[0]
            run,LS,lumi = tokens[0],tokens[1],tokens[6]
            LS = int(LS)
        else:
            run,lumi = tokens[0],tokens[1]
        run = int(run)
        lumi = float(lumi)
        #print(f"{run} -> {lumi}")
        if filterWithJson and not ROOT.isGoodRunLS(1, run, LS):
            continue
        for era in eras:
            if run in rfe[era]:
                lumiPerEra[era] += lumi
                break # if run was found, no need to check following eras

if args.json != "":
    print(f"Using JSON: {args.json}")
for era,intlumi in iter(lumiPerEra.items()):
    print(f"{era} : {intlumi}")
