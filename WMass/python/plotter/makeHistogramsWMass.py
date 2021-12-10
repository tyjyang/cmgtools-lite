#!/usr/bin/env python

# unpack histograms produce by mcPlots.py (or plotFakesTemplate.py in the case of Wmass analysis)
# The usually have TH3 with eta-pt-charge, or a 4th dimension with systematics, but we eventually save TH2 eta-pt for each charge separately, because it is easier to deal with the fit

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

sys.path.append(os.getcwd() + "/plotUtils/")
import utility as util

def getQCDScaleIndicesOldNtuples():
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

def getQCDScaleIndices():
    # branch LHEScaleWeight with 9 elements:                                                   *
    # [0] is mur=0.5 muf=0.5; [1] is mur=0.5 muf=1; [2] is mur=0.5 muf=2; [3] is mur=1 muf=0.5 ; [4] is mur=1 muf=1; [5] is mur=1 muf=2; [6] is mur=2 muf=0.5; [7] is mur=2 muf=1 ; [8] is mur=2 muf=2)*
    # we don't use the case where one scale goes up and the other down, so exclude pairs like (0.5,2)    
    ret = {  0 : "muRmuFDown",
             1 : "muRDown",
             3 : "muFDown",
             5 : "muFUp",
             7 : "muRUp",
             8 : "muRmuFUp"
    }
    return ret


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

def writeAndRemove(hists):
    for h in hists:
        h.Write()
        h.Delete()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile",  type=str, default=None, help="File to read histos from")
# changing this part: select output folder, rather than output file: the file will be created and named automatically for W or Z analysis based on --wlike, adding also the charge as appropriate: usually Wmunu_{charge}_shapes_{postfx}.root
#parser.add_argument("-o", "--outfile", type=str, default=None, help="output file name") 
parser.add_argument("--outdir", required=True, type=str, default=None, help="output folder")
parser.add_argument("-p",   "--postfix", type=str, default="", help="Postfix added to output file name, to avoid overwriting an existing one in case of tests");
parser.add_argument("--crop-negative-bin", dest="cropNegativeBin", action="store_true", help="Set negative bins to 0")
parser.add_argument("-v", "--verbose", type=int, default=3, choices=[0,1,2,3,4], help="Set verbosity level with logging, the larger the more verbose");
parser.add_argument("-c", "--charge", type=str, default=None, choices=["plus", "minus"], help="Charge for this channel")
parser.add_argument("--decorrelate-by-charge", dest="decorrByCharge", type=str, default=None, help="Matching regular expression for nuisances to be decorrelated by charge (or comma-separate list of expressions). The corresponding histograms have to be named according to the charge") 
parser.add_argument("--wlike", dest="isWlike", action="store_true", help="Flag for W-like analysis (have to change histogram name for signal to include the charge)")
parser.add_argument("--alpha-from-pdf-histogram", dest="alphaFromPdfHisto", action="store_true", help="Alpha is usually an independent histogram with respect to PDFs, but for NNPDF3.1, if this option is true, it is actually made from the TH3 containing PDFs (last two bins are alpha with 0.002 variation)")

args = parser.parse_args()

setLogging(args.verbose)

infilename = args.infile 
#outfilename = args.outfile

outdir = args.outdir + "/"
util.createPlotDirAndCopyPhp(outdir)

channelKeyword = "Zmumu" if args.isWlike else "Wmunu"
postfix = "" if not args.postfix else f"_{args.postfix}"
outfilename = outdir + f"{channelKeyword}_{args.charge}_shapes{postfix}.root"

if os.path.abspath(infilename) == os.path.abspath(outfilename):
    logging.warning(" input and output file names are the same. Abort")
    quit()

chargeBin = 1 if args.charge == "minus" else 2
    
if args.decorrByCharge:
    regexp = args.decorrByCharge.replace(',','|')
    matchDecorr = re.compile(regexp)
    chargeKey = "Plus" if args.charge == "plus" else "Minus"

hnomi = {} # {process name : histo}
hsyst = {} # {syst name : {process name : histo}}

# read input files to get all histograms, but do not append Up/Down in names yet
# project them already to select the desired charge
fin = ROOT.TFile.Open(infilename)
if not fin or not fin.IsOpen():
    raise(RuntimeError('Unable to open file {fn}'.format(fn=infilename)))
for ikey,e in enumerate(fin.GetListOfKeys()):
    name = e.GetName()
    obj  = e.ReadObj()
    if not obj:
        raise(RuntimeError('Unable to read object {n}'.format(n=name)))
    syst,proc = name.split("__")
    if args.alphaFromPdfHisto and any(x in syst for x in ["alphaS0117NNPDF31", "alphaS0119NNPDF31"]):
        continue
    if proc == "data":
        proc = "data_obs"
    elif proc == "Zmumu" and args.isWlike:
        proc = proc + "_" + args.charge 

    newname = "x_" + proc
    if name.startswith("nominal"):
        tracker = hnomi
    else:
        newname += f"_{syst}"
        if syst not in hsyst:
            hsyst[syst] = {}
        tracker = hsyst[syst]
    tracker[proc] = util.projectChargeHist(obj, args.charge, newname)

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

hasNNPDF31alphaS0001var = False
if all(x in systs for x in ["alphaS0117NNPDF31", "alphaS0119NNPDF31"]) and not args.alphaFromPdfHisto:
    hasNNPDF31alphaS0001var = True
    
outf = ROOT.TFile.Open(outfilename,'RECREATE')
if not outf or not outf.IsOpen():
    raise(RuntimeError('Unable to open file {fn}'.format(fn=outfilename)))
outf.cd()
for proc in list(hnomi.keys()):
    hnomi[proc].Write()

# histograms are TH3 eta-pt-charge or THn with dimension 4 embedding systematic variations
# we need tosplit charges
    
for syst in systs:
    if args.alphaFromPdfHisto and any(x in syst for x in ["alphaS0117NNPDF31", "alphaS0119NNPDF31"]):
        continue
    clabel = "" if not matchDecorr.match(syst) else chargeKey
    procs = list(hsyst[syst].keys())
    for proc in procs:
        h3D = hsyst[syst][proc] # this might actually be a THn with n=4
        #isTHn = False
        #if "THn" in h3D.ClassName():
        #    isTHn = True
        if "massWeight" in syst:
            # note, we are now using weights from LHEReweightingWeightCorrectMass
            # it has 21 elements for W and 23 for Z, where the last two values for Z should not be mass shift weights (maybe they are for width)
            # the sorting is based on increasing weight value, where the nominal==1 and is the 11th elements of the array 
            cenMassWgt = 10 # element xx=10 of LHEReweightingWeight[xx] is the 11th bin of the hist and the 10th entry of the projected array
            # 10 steps to 100 MeV
            groups = [[cenMassWgt-i, cenMassWgt+i] for i in range(1,11)]
            # This is a somewhat hacky way to get it to display the names properly (taking advantage of the indexing)
            writeAndRemove(util.makeVariationHistsForCharge(h3D, f"x_{proc}_massShift"+"{i}0MeV", groups))
        if "luminosity" in syst:
            writeAndRemove(util.makeVariationHistsForCharge(h3D, f"x_{proc}_lumi", [(0,1)]))
        if "effStatTnP" in syst:
            writeAndRemove(util.makeVariationHistsForCharge(h3D, f"x_{proc}_effStatTnP{clabel}", util.mirrorGroups(h3D), hnomi[proc], addMirror=True))
        if "effSystTnP" in syst:
            # here h3D is actually a TH2
            vars2d = [h3D, ROOT.mirrorHist(h3D, hnomi[proc], h3D.Clone("tmp"+syst))]
            writeAndRemove(util.buildVariationHistsForCharge(vars2d, f"x_{proc}_effSystTnP", [(0,1)]))
        if "CMS_scale_m" in syst:
            # Already has the mirrored hists
            writeAndRemove(util.makeVariationHistsForCharge(h3D, f"x_{proc}_muonScale", util.mirrorGroups(h3D, addMirror=False)))
            # As a cross check
            #writeAndRemove(util.makeVariationHistsForCharge(h3D, f"x_{proc}_muonScaleMirror", [(i,i+576) for i in range(288)], hnomi[proc], addMirror=True))
        if "muonL1PrefireStat" in syst:
            writeAndRemove(util.makeVariationHistsForCharge(h3D, f"x_{proc}_muonL1PrefireStat{clabel}", util.mirrorGroups(h3D), hnomi[proc], addMirror=True))
        if "muonL1PrefireSyst" in syst:
            # In this case the nominal is in bin 1
            writeAndRemove(util.makeVariationHistsForCharge(h3D, f"x_{proc}_muonL1PrefireSyst", [(1,2)]))
        if "muonPtScaleTest" in syst:
            for iptsyst in range(1, 96+1): # 48 *2, for up and down variations
                systname = "muonPtScaleTest%d" % (iptsyst if iptsyst <= 48 else (iptsyst-48))
                if matchDecorr.match(systname):
                    systname = systname + chargeKey
                name = "x_{p}_{s}{updown}".format(p=proc, s=systname, updown="Up" if iptsyst <=48 else "Down")
                h2D = ROOT.projectTH2FromTH3(h3D, name, iptsyst, iptsyst)
                h2D.Write()
        if "qcdScale" in syst:
            if "qcdScaleVptBin" in syst:
                ptbin = syst.split("VptBin")[1] # value starts from 1, so can use 0 to signal its absence
            else:
                ptbin = 0 # told you ;)
                # new ntuples have 9 elements, older had 18
            qcdscales = getQCDScaleIndices() if h3D.GetNbinsZ() == 9 else getQCDScaleIndicesOldNtuples()
            indices = sorted(list(qcdscales.keys()))
            for i in indices:
                systname = qcdscales[i]
                if ptbin:
                    tmp = systname.replace("Up","").replace("Down","") + str(ptbin) + ("Up" if systname.endswith("Up") else "Down")
                    systname = tmp 
                if matchDecorr.match(systname):
                    systname = systname.replace("Up", chargeKey + "Up") if systname.endswith("Up") else systname.replace("Down", chargeKey + "Down")
                name = "x_" + proc + "_" + systname
                h2D = ROOT.projectTH2FromTH3(h3D, name, i+1, i+1) # root histogram bin number starts from 1
                h2D.Write()
        #if "pdf" in syst: # now we have more PDF sets, so names are more complicated
        # some temporary hacks to force using alpha for 0.001 variations if present (those in the pdf histogram are the 0.002 variations)
        if "pdfNNPDF31" in syst:
            # this includes actual pdf hessians (bins 1 to 100) and alphaSDown and alphaSUp by 0.002 (bin 101 and 102)
            writeAndRemove(util.makeVariationHistsForCharge(h3D, f"x_{proc}_pdf"+"{i}NNPDF31", util.mirrorGroups(h3D)[:100], hnomi[proc], True))
        elif "pdfCT18" in syst:
            writeAndRemove(util.makeVariationHistsForCharge(h3D, f"x_{proc}_"+"pdf{i}CT18", util.pairGroups(58)))
        if "alphaS" in syst:
            # At some point need to scale these, I guess in the cardmaker
            writeAndRemove(util.makeVariationHistsForCharge(h3D, f"x_{proc}_pdfAlphaS002", [(0,1)]))
                    
nKeys = outf.GetNkeys()
outf.Close()
print(f"{nKeys} histograms saved in file {outfilename}")

