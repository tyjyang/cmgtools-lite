#!/usr/bin/env python
import os, sys
import argparse
import subprocess

from background_samples import campaign, samples

parser = argparse.ArgumentParser("")
parser.add_argument("-i", "--infile",  type=str, default=None, help="File to read samples from")
parser.add_argument("-e", "--era",     type=str, choices=["all", "preVFP", "postVPF"], default="all", help="Era to use")
parser.add_argument("-t", "--tier",    type=str, choices=["NANOAODSIM"], default="NANOAODSIM", help="Dataset tier")
parser.add_argument("-d", "--dry_run", dest="dryRun", action="store_true",    help="Print but do not run")
parser.add_argument("-s", "--save",    type=str, default=None,   help="Save list of files from DAS to txt files, named after the sample and campaign, in this folder")
parser.add_argument("-c", "--copy",    type=str, default=None,   help="Call copyFileXrdcp.py to copy files to this output folder (sample and campaign name is used as subfolder). Only works in combination with option -s")
args = parser.parse_args()

if args.copy and not args.save:
    print("Error: option -c is only used in combination with option -c")
    quit()
    
outdir = "./"
if args.save:
    outdir = args.save
    if not outdir.endswith("/"):
        outdir += "/"
    os.system(f"mkdir -p {outdir}") 

if args.copy and not os.path.exists(args.copy):
    if not args.copy.endswith("/"):
        args.copy += "/"
    os.system(f"mkdir -p {args.copy}") 

filesToCopy = []
    
for c in campaign.keys():
    if args.era != "all" and args.era != c:
        continue
    for s in samples.keys():
        dataset = f"/{samples[s]}/{campaign[c]}/NANOAODSIM" 
        command = "/cvmfs/cms.cern.ch/common/dasgoclient -query='file dataset={d}'".format(d=dataset)
        files = []
        if args.dryRun:
            print(command)
        else:
            files = subprocess.check_output(command, shell=True, universal_newlines=True).strip().split()
        if args.save:
            fname = f"{outdir}{s}_{c}.txt"
            if len(files):
                if args.dryRun:
                    outf = open(fname, "w")
                    for f in files:
                        outf.write(f"{f}\n")
                    outf.close()
                    print(f"Files written in {fname}")
                filesToCopy.append(fname)
            else:
                print(">>> Warning: no files for this campaign. No output file written")
                print(f">>> Campaign: {c}   Sample: {s}")
                print('-'*30)
        else:
            print('='*30)
            print(f"Campaign: {c}   Sample: {s}")
            print('-'*30)
            for f in files:
                print(f)

if args.copy:
    for fname in filesToCopy:
        outfolder = args.copy + os.path.basename(fname).split(".txt")[0] + "/"
        command = f"python3 copyFileXrdcp.py {fname} {outfolder}"
        if args.dryRun:
            print(command)
        else:
            os.system(command)
