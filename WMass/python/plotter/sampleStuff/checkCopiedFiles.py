#!/usr/bin/env python3

# after using sampleStuff/getSampleFiles.py
# python3 sampleStuff/checkCopiedFiles.py -i sampleStuff/ -o /data/shared/originalNANO/

import os
import re
import argparse
import subprocess
from background_samples import campaign, samples
import copyFileXrdcp

import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True


parser = argparse.ArgumentParser("")
parser.add_argument("-i", "--inputdir",  type=str, default=None, help="Input folder with files containing files that where copied (which need to be checked)")
parser.add_argument("-o", "--outputdir",  type=str, default=None, help="Output folder with copied files")
parser.add_argument("-r", "--regexp",  type=str, default=None, help="Regular expression to filter samples")
parser.add_argument("-x", "--exclude-regexp", dest="excludeRegexp", type=str, default=None, help="Regular expression to exclude samples")
parser.add_argument("-e", "--era",     type=str, choices=["all", "preVFP", "postVPF"], default="all", help="Era to use")
parser.add_argument("-d", "--dry_run", dest="dryRun", action="store_true",    help="Print but do not run")
args = parser.parse_args()

if not args.outputdir.endswith('/'):
    args.outputdir += '/'
if not args.inputdir.endswith('/'):
    args.inputdir += '/'

regexp = None
if args.regexp:
    regexp = re.compile(args.regexp)

antiregexp = None
if args.excludeRegexp:
  antiregexp = re.compile(args.excludeRegexp)
    
for c in campaign.keys():
    if args.era != "all" and args.era != c:
        continue
    for s in samples.keys():
        folder = f"{s}_{c}"
        if regexp:
            if not regexp.match(folder):
                continue
        if antiregexp:
            if antiregexp.match(folder):
                continue
        fname = f"{folder}.txt"
        fullInputFile = args.inputdir + fname
        fullOutputDir = args.outputdir + folder + '/'
        missingFileList = []
        #print(fullInputFile)
        f = open(fullInputFile)
        if f:
            lines = [x.strip() for x in f.readlines()]
            for line in lines:
                bname = os.path.basename(line)
                #print(line)
                if not os.path.exists(fullOutputDir + bname):
                    #print(fullOutputDir + bname)
                    missingFileList.append(line)
                else:
                    rf = ROOT.TFile.Open(fullOutputDir + bname, "READ")
                    if not rf or rf.IsZombie():
                        print("ERROR! Unable to open file %s" % fullOutputDir + bname)
                        missingFileList.append(line)
                    elif rf.TestBit(ROOT.TFile.kRecovered):
                        print("WARNING! Attempt was made to recover file %s" % fullOutputDir + bname)
                        missingFileList.append(line)
        else:
            print(f"Error in opening file {fullInputFile}")
            quit()
        print('-'*30)
        print(f"{folder}: missing {len(missingFileList)} files")
        print('-'*30)
        tmpfilename = fullInputFile.replace(".txt", "_TMPTOREMOVE.txt")
        tmpfile = open(tmpfilename, "w")
        if not tmpfile:
            print(f"Error opening file {tmpfilename}")
            quit()
        for line in missingFileList:
            if line.startswith("/store/"):
                infname = f"root://xrootd-cms.infn.it//{line}"
            else:
                infname = line
            tmpfile.write(infname + "\n")
        tmpfile.close()
        copyFileXrdcp.copyFile(tmpfilename, fullOutputDir, dryRun=args.dryRun)
        os.system(f"rm {tmpfilename}")
