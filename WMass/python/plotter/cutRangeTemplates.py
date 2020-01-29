#!/usr/bin/env python                                                       

#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np
import root_numpy

from w_helicity_13TeV.utilities import util
utilities = util()

from w_helicity_13TeV.make_diff_xsec_cards import getXYBinsFromGlobalBin
from w_helicity_13TeV.make_diff_xsec_cards import getArrayParsingString
from w_helicity_13TeV.make_diff_xsec_cards import getDiffXsecBinning
from w_helicity_13TeV.make_diff_xsec_cards import templateBinning
from w_helicity_13TeV.make_diff_xsec_cards import get_ieta_ipt_from_process_name
from w_helicity_13TeV.make_diff_xsec_cards import get_ieta_from_process_name
from w_helicity_13TeV.make_diff_xsec_cards import get_ipt_from_process_name

from w_helicity_13TeV.templateRolling import roll1Dto2D, dressed2D
from w_helicity_13TeV.rollingFunctions import unroll2Dto1D

def getHistoNarrowPtRange(h,inbinning,outbinning):
    # exploit the fact that eta binning is the same
    # and the unrolling is made with consecutive pt bins (i.e. last bins are those at hig pt)
    nEta = inbinning[0]
    nPtLarge = inbinning[2]
    nPtNarrow = outbinning[2]
    nBinsOut = nEta * nPtNarrow
    hret = ROOT.TH1D(h.GetName()+"_cutPt", h.GetTitle(), nBinsOut, 0.5, 0.5 + float(nBinsOut))
    for i in range(1, 1+ nBinsOut):
        hret.SetBinContent(i, h.GetBinContent(i))
        hret.SetBinError(i, h.GetBinError(i))
    return hret

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options]")
    #parser.add_option("-i", "--indir",    dest="indir",  type="string", default="", help="Input folder");
    parser.add_option("-i", "--indir",   dest="indir", type="string", default="", help="Input folder name");
    parser.add_option("-o", "--outdir",   dest="outdir", type="string", default="./", help="Output folder (current one as default)");
    #parser.add_option("-p", "--processes",   dest="processes", type="string", default=".*", help="Regular expression for processes");
    #parser.add_option("--xp", "--exclude-processes",   dest="excludeProcesses", type="string", default="x_W.*", help="Regular expression to exclude processes");
    parser.add_option('-b','--etaPtbinning', dest='etaPtbinning', default='', type='string', help='eta-pt binning for templates, with reduced range')
    parser.add_option("-f", "--flavour", dest="flavour", type="string", default='el', help="Channel: either 'el' or 'mu'");
    #parser.add_option("-c", "--charge",  dest="charge",  type="string", default='', help="Charge: either 'plus' or 'minus'");
    (options, args) = parser.parse_args()

    if not len(options.indir):
        print "Warning: you must specify an input folder name with option --indir"
        quit()

    # manage output folder
    outdir = options.outdir
    if not outdir.endswith('/'): outdir += "/"
    if outdir != "./":
        if not os.path.exists(outdir):
            print "Creating folder", outdir
            os.system("mkdir -p " + outdir + "/part0/")

    if options.flavour not in ["el", "mu"]:
       print "Warning: you must specify a lepton flavour with option -f el|mu"
       quit()
    #if options.charge not in ["plus", "minus"]:
    #   print "Warning: you must specify a charge with option -c plus|minus"
    #   quit()

    #charge = options.charge
    flavour = options.flavour
    if ("_"+flavour+"_") not in options.indir:
        print "Warning: flavour not appearing in input folder. Did you set it correctly?"
        quit()
    isMu = True if flavour == "mu" else False


    inbinfile = options.indir + "/binningPtEta.txt"
    inetaPtBinningVec = getDiffXsecBinning(inbinfile, "reco")
    inrecoBins = templateBinning(inetaPtBinningVec[0],inetaPtBinningVec[1])
    inbinning = [inrecoBins.Neta, inrecoBins.etaBins, inrecoBins.Npt, inrecoBins.ptBins]

    outbinfile = options.etaPtbinning
    outetaPtBinningVec = getDiffXsecBinning(outbinfile, "reco")
    outrecoBins = templateBinning(outetaPtBinningVec[0],outetaPtBinningVec[1])
    outbinning = [outrecoBins.Neta, outrecoBins.etaBins, outrecoBins.Npt, outrecoBins.ptBins]

    rootfiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(options.indir,followlinks=True) for f in fn if f.endswith('.input.root')]
    for f in rootfiles:
        print "-"*30
        print f
        # now open file, and open new one which will have the templates with different range
        tfno = ROOT.TFile.Open(f,"READ")
        if not tfno or not tfno.IsOpen():
            raise RuntimeError('Unable to open file {fn}'.format(fn=f))
        fname = os.path.basename(f)
        fout = outdir + "/part0/" + fname
        of = ROOT.TFile(fout,'recreate')
        if not of or not of.IsOpen():
            raise RuntimeError('Unable to open file {fn}'.format(fn=fout))

        nKeys = tfno.GetNkeys()
        nCopiedKeys = 0
        for ikey,e in enumerate(tfno.GetListOfKeys()):
            name = e.GetName()
            obj  = e.ReadObj()
            objnew = getHistoNarrowPtRange(obj,inbinning,outbinning)
            objnew.Write(name)
            nCopiedKeys += 1
        print "Overwrote {n}/{tot} keys to {of}".format(n=str(nCopiedKeys),tot=str(nKeys), of=fout)
        of.Close()
        tfno.Close()

print ""
print "Created {n} files".format(n=len(rootfiles))
print ""

print ""
print "="*30
print "DONE"
print "="*30
print ""
