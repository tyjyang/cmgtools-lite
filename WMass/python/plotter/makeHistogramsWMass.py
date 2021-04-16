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
parser.add_argument("-c", "--charge", type=str, default=None, choices=["plus", "minus"], help="Charge for this channel")
parser.add_argument("--decorrelate-by-charge", dest="decorrByCharge", type=str, default=None, help="Matching regular expression for nuisances to be decorrelated by charge (or comma-separate list of expressions). The corresponding histograms have to be named according to the charge") 
parser.add_argument("--wlike", dest="isWlike", action="store_true", help="Flag for W-like analysis (have to change histogram name for signal to include the charge)")

args = parser.parse_args()

setLogging(args.verbose)

infilename = args.infile 
outfilename = args.outfile
if os.path.abspath(infilename) == os.path.abspath(outfilename):
    logging.warning(" input and output file names are the same. Abort")
    quit()

if args.decorrByCharge:
    regexp = args.decorrByCharge.replace(',','|')
    matchDecorr = re.compile(regexp)
    chargeKey = "Plus" if args.charge == "plus" else "Minus"
    
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
    if proc == "Zmumu" and args.isWlike:
        proc = proc + "_" + args.charge 
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
print("Systematics (as named in original file) and relevant processes")
systs = sorted(list(hsyst.keys()))
for syst in systs:
    proclist = ", ".join(str(x) for x in sorted(hsyst[syst].keys()))
    print("{: <20}: {p}".format(syst,p=proclist))
print('-'*30)

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
        if "massWeight" in syst:
            # check https://github.com/WMass/nanoAOD-tools/blob/master/python/postprocessing/wmass/lheWeightsFlattener.py
            cenMassWgt = 11 # element 11 of LHEReweightingWeight is the 12th of the array, so it is bin 12 of the histogram (TO BE CHECKED AGAIN)
            maxMassShift = 100
            massGrid = 10
            for i in range(1, 1 + int(maxMassShift/massGrid)):
                massShift = i * massGrid
                ibinUp = cenMassWgt + 1 + i 
                name = "x_" + proc + "_massShift%dMeVUp" % massShift
                h2DUp = getTH2fromTH3(h3D, name, ibinUp, ibinUp)
                ibinDown = cenMassWgt + 1 - i 
                name = "x_" + proc + "_massShift%dMeVDown" % massShift
                h2DDown = getTH2fromTH3(h3D, name, ibinDown, ibinDown)
                h2DUp.Write()
                h2DDown.Write()           
        if "effStatTnP" in syst:
            for ieff in range(1, 576+1): # need to be kept manually consistent until we save these numbers somewhere
                systname = "effStatTnP%d" % ieff
                if matchDecorr.match(systname):
                    systname = systname + chargeKey
                name = "x_{p}_{s}Up".format(p=proc, s=systname)  # define this as Up variation 
                h2D = getTH2fromTH3(h3D, name, ieff, ieff)
                h2D_mirror = h2D.Clone(name.replace("Up", "Down"))
                h2D_mirror = mirrorShape(hnomi[proc], h2D, h2D_mirror)
                h2D.Write()
                h2D_mirror.Write()
        if "muonL1Prefire" in syst:
            for ieff in range(1, 16+1):
                systname = "muonL1Prefire%d" % ieff
                if matchDecorr.match(systname):
                    systname = systname + chargeKey
                name = "x_{p}_{s}Up".format(p=proc, s=systname)  # define this as Up variation 
                h2D = getTH2fromTH3(h3D, name, ieff, ieff)
                h2D_mirror = h2D.Clone(name.replace("Up", "Down"))
                h2D_mirror = mirrorShape(hnomi[proc], h2D, h2D_mirror)
                h2D.Write()
                h2D_mirror.Write()
        if "qcdScale" in syst:
            qcdscales = getQCDScaleIndices()
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
                if matchDecorr.match(systname):
                    systname = systname.replace("Up", chargeKey + "Up") if systname.endswith("Up") else systname.replace("Down", chargeKey + "Down")
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

nKeys = outf.GetNkeys()
outf.Close()
print(f"{nKeys} histograms saved in file {outfilename}")


