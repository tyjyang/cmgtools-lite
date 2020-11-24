#!/usr/bin/env python  

# script to check which ntuples are missing for VBF Higgs analysis
# they are usually named as XXX_N.root, where N goes from 0 to a given number
# this script check the maximum N that is available and assess all n with 0<=n<=N that are missins
# one should filter names using a regular expression passed to option -r

from shutil import copyfile
import re, sys, os, os.path, subprocess, json, ROOT
import numpy as np

from os import listdir
from os.path import isfile, join

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] folder")
#parser.add_option('-x', '--exclude-regexp', dest='exclude_regexp', default='', type='string', help='Like option -r, but will reject object if name matches')
parser.add_option('-s', '--silent'  ,       dest='silent',         default=False, action='store_true', help='Silent mode: just print summary, not all entries')  
parser.add_option('-r', '--regexp',         dest='regexp',         default='', type='string', help='Filter output passing a regular expression to be matched in the object name')
(options, args) = parser.parse_args()

if len(sys.argv) < 1:
    parser.print_usage()
    quit()

indir = args[0]
if not indir.endswith("/"):
    indir += "/"

files = [f for f in listdir(indir) if isfile(join(indir, f)) and f.endswith(".root")]
#print "="*40
#for f in files:
#    print f
if len(options.regexp):
    prunedFiles = [f for f in files if re.match(options.regexp,f)]
    files = prunedFiles
files = sorted(files, key= lambda x: int(x.split('_')[-1].rstrip(".root")))
if not options.silent:
    print "="*40
    for f in files:
        print f
maxN = int(files[-1].split('_')[-1].rstrip(".root"))
filePrefix = "_".join(files[-1].split('_')[:-1])
missingN = []
for n in range(maxN+1):
    checkFile = indir + filePrefix + "_" + str(n) + ".root"
    if not isfile(checkFile):
        missingN.append(n)

print "-"*40
print "File pattern: %s" % filePrefix
if len(missingN):
    print "Following %d/%d elements are missing" % (len(missingN), maxN)
    print missingN
else:
    print "No missing elements (assuming they were really %d in total)" % maxN


