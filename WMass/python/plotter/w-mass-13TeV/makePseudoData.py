#!/usr/bin/env python3

# script to modify the root files with the histograms for the fit, creating a pseudodata histogram to be used in the fit with -t 0
# when preparing the fit, text2hdf5.py has an option -D imported from HiggsAnalysis/CombinedLimit/python/DatacardParser.py to define the name of the dataset
# therefore, the original file can be safely modified since one is only adding a new histogram, but a copy is kept as backup

# example
# python w-mass-13TeV/makePseudoData.py plots/testNanoAOD/testNtuplesAltPDF/testCode_TH3_chargePlus/plotsPDF/shapes_Wmunu_plus_postVFP.root cards/wmass_fixMassWeights_splitW/Wmunu_plus_shapes.root -p pseudodataWmunuFromNNPDF30 -n pseudodata --add-modified-shapes etaPt_Wmunu_plus_postVFP_nomi_nnpdf30_2d --add-original-shapes x_Wmunu_minus,x_Wtaunu_plus,x_Wtaunu_minus,x_Zmumu,x_Ztautau,x_Top,x_Diboson,x_data_fakes

import os, re
import argparse
from array import array

import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

import utilities
utilities = utilities.util()

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(1)
ROOT.TH1.SetDefaultSumw2()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input file with histogram to take to form new dataset")
    parser.add_argument("originalshapefile", type=str, nargs=1, help="Root file where to add the new pseudodata histogram")
    parser.add_argument('-p','--postfix', dest='postfix', default='addPseudodata', type=str, help='Postfix for output file')
    parser.add_argument('-n','--name', default='x_pseudodata', type=str, help='name for pseudodata histogram to be created')
    parser.add_argument(     '--add-original-shapes', dest='originalShapes', default='x_Wmunu_minus,x_Wtaunu_plus,x_Wtaunu_minus,x_Zmumu,x_Ztautau,x_Top,x_Diboson,x_data_fakes', type=str, help='Histograms from original root file that will be added to form the pseudodata (comma-separated list of histogram names)')
    parser.add_argument(     '--add-modified-shapes', dest='modifiedShapes', default='x_Wmunu_plus', type=str, help='Histograms from new root file that will be added to form the pseudodata (comma-separated list of histogram names, if empty nothing is used)')
    parser.add_argument(     '--overwrite', action='store_true', help="If the new created histogram already exists in the input it will be overwritten if using this option")
    args = parser.parse_args()

    if not args.name.startswith("x_"):
        args.name = "x_" + args.name
    
    pseudodata = None    
    infile = safeOpenFile(args.originalshapefile[0], mode="READ")
    # check that the pseudodata histogram does not exist, unless one wants to update it
    htmp = safeGetObject(infile, args.name, quitOnFail=False, silent=True)
    if htmp and not args.overwrite:
        print("Warning: histogram {args.name} already exists in {args.originalshapefile[0]}. Abort without overwriting")
        quit()

    for ih,hname in enumerate(list(args.originalShapes.split(","))):
        if ih == 0:
            pseudodata = safeGetObject(infile, hname)
        else:
            pseudodata.Add(safeGetObject(infile, hname))
    infile.Close()

    infile = safeOpenFile(args.rootfile[0], mode="READ")
    for ih,hname in enumerate(list(args.modifiedShapes.split(","))):
        pseudodata.Add(safeGetObject(infile, hname))
    infile.Close()

    pseudodata.SetName(args.name)
    pseudodata.SetTitle("Pseudodata")
    
    outfilename = args.originalshapefile[0].replace(".root", f"_{args.postfix}.root")
    print(f"Creating updated root file {outfilename}")
    print(f"Adding histogram {args.name}")
    os.system(f"cp {args.originalshapefile[0]} {outfilename}")
    outfile = safeOpenFile(outfilename, mode="UPDATE")
    pseudodata.Write()
    outfile.Close()

