#!/usr/bin/env python3

from shutil import copyfile
import re, sys, os, os.path
import subprocess
from concurrent.futures import ThreadPoolExecutor
import argparse

def xrdcp(files):
    subprocess.run(["xrdcp", "-f", files[0], files[1]])

def copyFile(inputfile, outdir, jobs=32, dryRun=False):

    outname = outdir
    if not outname.endswith("/"):
        outname += "/"
    os.system("mkdir -p {outdir}".format(outdir=outname))

    infiles = []
    outfiles = []

    with open(inputfile) as f:
        for line in f:
            if line.startswith("/store/"):
                fname = f"root://xrootd-cms.infn.it//{line}"
            else:
                fname = line
            fname = fname.strip()
            infiles.append(fname)
            outfiles.append(outname+os.path.basename(fname))

    if dryRun:
        for infile,outfile in zip(infiles,outfiles):
            print(f"xrdcp {infile} {outfile}")
    else:
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            executor.map(xrdcp, zip(infiles,outfiles))


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", type=str, nargs=1, help="Input file with files");
    parser.add_argument("outdir",    type=str, nargs=1, help='output directory to save things (should end with /)')
    parser.add_argument('-d', '--dry-run', dest='dryRun', action='store_true', help='If true, print commands but do not execute them')
    #parser.add_argument('--options', type=str, default="", help='Options to be passed to xrdcp')
    parser.add_argument("-j","--jobs", type=int, default=32)
    args = parser.parse_args()

    copyFile(args.inputfile[0], args.outdir[0], jobs=args.jobs, dryRun=args.dryRun)
    
