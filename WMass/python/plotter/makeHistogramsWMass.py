#!/usr/bin/env python

#from mcPlots import *
import re, os, os.path
import math
import copy
import argparse
import logging

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True


from tree2yield import setLogging
from cropNegativeTemplateBins import cropNegativeContent

def getQCDScaleIndices():
    # first number is for renormalization scale, the other is for factorization scale
    # [("05","05"), ("05","1"), ("05", "2"), ("1", "05"), \
    #  ("1","1"), ("1","2"), ("2", "05"), ("2", "1"), ("2", "2")]):
    # and use even indices in the LHEScaleWeight array for NNPDF3.1 (odd are for NNPDF3.0)
    # so e.g. muRmuFUp is (2,2), corresponding to index 8*2 (index for first pair is 0)
    # we don't use the case where one scale goes up and the other down, so exclude pairs like (0.5,2)
    ret = {  0 : "muRmuFDown",
             2 : "muRDown",
             6 : "muFDown",
            10 : "muFUp",
            14 : "muRUp",
            16 : "muRmuFUp"
    }
    return ret

def getTH2fromTH3(hist3D, name, binStart, binEnd=None):

    if binEnd == None:
        binEnd = binStart
    hist3D.GetZaxis().SetRange(binStart,binEnd)
    # Order yx matters to have consistent axes!
    hist2D = hist3D.Project3D("yxe") # make TH2 with y axis versus x axis
    hist2D.SetName(name)
    return hist2D
    
def mirrorShape(nominal,alternate,mirror):
    # assumes any regularization (e.g. cropping negative bin content to make it 0) already  happened outside
    # same for normalization

    # FIXME: need to add possible conditions when either y0 or yA are negative
    # while this choice is only relevant for some bins with low stat in MC, it can induce non trivial 
    # effects on constraints of some particular nuisance parameters

    for bx in range(1,nominal.GetNbinsX()+1):
        for by in range(1,nominal.GetNbinsY()+1):
            y0 = nominal.GetBinContent(bx, by)
            yA = alternate.GetBinContent(bx, by)
            yM = max(0,2*y0-yA)
            mirror.SetBinContent(bx, by, yM)
            mirror.SetBinError(bx, by, alternate.GetBinError(bx, by))
    return mirror


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile",  type=str, default=None, help="File to read histos from")
parser.add_argument("-o", "--outfile", type=str, default=None, help="output file name") 
#parser.add_argument("--od", "--outdir", type=str, default=None, help="output folder") 
parser.add_argument("--crop-negative-bin", dest="cropNegativeBin", action="store_true", help="Set negative bins to 0")
parser.add_argument("-v", "--verbose", type=int, default=3, choices=[0,1,2,3,4], help="Set verbosity level with logging, the larger the more verbose");

args = parser.parse_args()

setLogging(args.verbose)

infilename = args.infile 
outfilename = args.outfile
if os.path.abspath(infilename) == os.path.abspath(outfilename):
    logging.warning(" input and output file names are the same. Abort")
    quit()

hnomi = {} # {process name : histo}
hsyst = {} # {syst name : {process name : histo}}

# read input files to get all histograms, but do not append Up/Down in names yet
fin = ROOT.TFile.Open(infilename)
if not fin or not fin.IsOpen():
    raise(RuntimeError('Unable to open file {fn}'.format(fn=infilename)))
for ikey,e in enumerate(fin.GetListOfKeys()):
    name = e.GetName()
    obj  = e.ReadObj()
    if not obj:
        raise(RuntimeError('Unable to read object {n}'.format(n=name)))
    syst,proc = name.split("__")
    if proc == "data":
        proc = "data_obs"
    newname = "x_" + proc
    if name.startswith("nominal"):
        hnomi[proc] = obj.Clone(newname)
        hnomi[proc].SetDirectory(0)
    else:
        if syst not in hsyst:
            hsyst[syst] = {}
        hsyst[syst][proc] = obj.Clone(newname+"_"+syst)
        hsyst[syst][proc].SetDirectory(0)
fin.Close()

print('-'*30)
print("Processes: %s" % ", ".join(str(x) for x in list(hnomi.keys())))
print('-'*30)
print("Systematics and relevant processes")
systs = sorted(list(hsyst.keys()))
for syst in systs:
    proclist = ", ".join(str(x) for x in sorted(hsyst[syst].keys()))
    print("{: <20}: {p}".format(syst,p=proclist))
print('-'*30)

qcdscales = getQCDScaleIndices()

outf = ROOT.TFile.Open(outfilename,'RECREATE')
if not outf or not outf.IsOpen():
    raise(RuntimeError('Unable to open file {fn}'.format(fn=outfilename)))
outf.cd()
for proc in list(hnomi.keys()):
    hnomi[proc].Write()

for syst in systs:
    procs = list(hsyst[syst].keys())
    for proc in procs:
        h3D = hsyst[syst][proc]
        if "qcdScale" in syst:
            indices = sorted(list(qcdscales.keys()))
            if "qcdScaleVptBin" in syst:
                ptbin = syst.split("VptBin")[1] # value starts from 1, so can use 0 to signal its absence
            else:
                ptbin = 0 # told you ;)
            for i in indices:
                systname = qcdscales[i]
                if ptbin:
                    tmp = systname.replace("Up","").replace("Down","") + str(ptbin) + ("Up" if systname.endswith("Up") else "Down")
                    systname = tmp 
                name = "x_" + proc + "_" + systname
                h2D = getTH2fromTH3(h3D, name, i+1, i+1) # root histogram bin number starts from 1
                h2D.Write()
        if "pdf" in syst:
            # this includes actual pdf hessians (bins 1 to 100) and alphaSUp and alphaSDown (bin 101 and 102)
            # pdfxx needs mirroring, alphaS already has Up and Down
            for i in range(1,103):
                if i <= 100:
                    name = "x_" + proc + "_pdf%dUp" % i # define this as Up variation 
                    h2D = getTH2fromTH3(h3D, name, i, i)
                    h2D_mirror = h2D.Clone(name.replace("Up", "Down"))
                    h2D_mirror = mirrorShape(hnomi[proc], h2D, h2D_mirror)
                    h2D.Write()
                    h2D_mirror.Write()
                else:
                    name = "x_" + proc + "_alphaS%s" % ("Up" if i == 101 else "Down")
                    h2D = getTH2fromTH3(h3D, name, i, i)
                    h2D.Write()
                                           
outf.Close()
print(f"Histograms saved in file {outfilename}")


