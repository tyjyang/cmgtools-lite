#!/usr/bin/env python                                                       

# script to scale content of histograms in shapes.root.
 
#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np

from w_helicity_13TeV.make_diff_xsec_cards import getXYBinsFromGlobalBin
from w_helicity_13TeV.make_diff_xsec_cards import getArrayParsingString
from w_helicity_13TeV.make_diff_xsec_cards import getDiffXsecBinning
from w_helicity_13TeV.make_diff_xsec_cards import getArrayBinNumberFromValue
from w_helicity_13TeV.make_diff_xsec_cards import templateBinning

from w_helicity_13TeV.templateRolling import roll1Dto2D, dressed2D
from w_helicity_13TeV.rollingFunctions import unroll2Dto1D

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
#parser.add_option("-i", "--indir",    dest="indir",  type="string", default="", help="Input folder");
parser.add_option("-f", "--infile",   dest="infile", type="string", default="", help="Input file name");
parser.add_option("-o", "--outdir",   dest="outdir", type="string", default="./", help="Output folder (current one as default)");
parser.add_option("-m", "--match",   dest="match",   type="string", default="x_W.*_lep_.*_ipt_.*",    help="Regular expression for match");
parser.add_option("-e", "--exclue-match",   dest="excludeMatch",   type="string", default="data.*",    help="Regular expression for anti-match");
parser.add_option("-n", "--name",   dest="name",   type="string", default="",    help="Name for output file (if not given use default one)");
parser.add_option("-s", "--scale",   dest="scale",   type="float", default=1.0, help="Scale factor");
# other options to scale some bins based on a histogram
parser.add_option("--scaleFromHist",   dest="scaleFromHist",   type="string", default="",    help="pass name of root file and 2D histogram inside it, comma separated. The binning of it should be consistent with the shapes when rolled back to a 2D");
parser.add_option(      '--binfile'  , dest='binfile', default='binningPtEta.txt', type='string', help='eta-pt binning for templates, by default it is expected to be in input folder (needed when scaling only some bins)')
parser.add_option(       '--eta-range', dest='eta_range', action="append", type="float", nargs=2, default=[], help='Will only scale reco eta bins in this range. Takes two float as arguments (increasing order) and can specify multiple times. They should match bin edges and a bin is not scaled if at least one edge is outside this range (so can set left edge as 1.99 to scale from eta=2.0')
parser.add_option(       '--pt-range', dest='pt_range', action="append", type="float", nargs=2, default=[], help='Will only scale reco pt bins in this range. Takes two float as arguments (increasing order) and can specify multiple times. They should match bin edges and a bin is not scaled if at least one edge is outside this range (so can set left edge as 29.9 to scale from pt=30.0')
parser.add_option(     '--suppress-root-warnings', dest='suppressRootWarnings' , default=False , action='store_true',   help='Suppress root warnings (it happens with TH1::Multiply that two identical histograms are seen as having different binning, it depends on how they were created, so you can safely suppress the warnings)')
#parser.add_option("-f", "--flavour", dest="flavour", type="string", default='el', help="Channel: either 'el' or 'mu'");
# charge is only used for the file name, not for the match
#parser.add_option("-c", "--charge",  dest="charge",  type="string", default='', help="Charge: either 'plus' or 'minus'");
(options, args) = parser.parse_args()

if options.suppressRootWarnings:
    #ROOT.gSystem.RedirectOutput("mylogfile_rootwarnings.txt", "w" if ikey == 0 else "a")
    savErrorLevel = ROOT.gErrorIgnoreLevel; 
    ROOT.gErrorIgnoreLevel = ROOT.kError;

# if not len(options.indir):
#     print "Warning: you must specify a folder with signal root files with option --indir"
#     quit()
if not len(options.infile):
    print "Warning: you must specify an input file name with option --infile"
    quit()

# manage output folder
outdir = options.outdir
if not outdir.endswith('/'): outdir += "/"
if outdir != "./":
    if not os.path.exists(outdir):
        print "Creating folder", outdir
        os.system("mkdir -p " + outdir)

#if options.flavour not in ["el", "mu"]:
#    print "Warning: you must specify a lepton flavour with option -f el|mu"
#    quit()
#if options.charge not in ["plus", "minus"]:
#    print "Warning: you must specify a charge with option -c plus|minus"
#    quit()

#charge = options.charge
#flavour = options.flavour

h2 = None
h2_unrolledTo1D = None
if options.scaleFromHist:
    fn = options.scaleFromHist.split(',')[0]
    hn = options.scaleFromHist.split(',')[1]
    rf = ROOT.TFile.Open(fn,"READ")
    if not rf or not rf.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=fn))
    h2 = rf.Get(hn)
    if not h2:
        raise RuntimeError('Unable to get histogram {hn}'.format(hn=hn))
    h2.SetDirectory(0)
    rf.Close()
    print "h2 has %d * %d bins" % (h2.GetNbinsX(),h2.GetNbinsY())
    h2_unrolledTo1D = unroll2Dto1D(h2, newname=h2.GetName()+"_unrolledTo1D")
    print "h2_unrolledTo1D has %d bins from %.1f to %.1f" % (h2_unrolledTo1D.GetNbinsX(),
                                                             h2_unrolledTo1D.GetXaxis().GetBinLowEdge(1),
                                                             h2_unrolledTo1D.GetXaxis().GetBinUpEdge(h2_unrolledTo1D.GetNbinsX()))
    # set error to 0 to allow multiplication without affecting uncertainty of original histogram (not fully ok, but simpler)
    for ib in range(0,2+h2_unrolledTo1D.GetNbinsX()):
        h2_unrolledTo1D.SetBinError(ib, 0.0) 


#infile = options.indir + "/" + options.infile
infile = options.infile

etaPtBinningFile = os.path.dirname(infile) + "/" + options.binfile
# get eta-pt binning for both reco and gen
etaPtBinningVec = getDiffXsecBinning(etaPtBinningFile, "reco")
recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])
nEtabins = recoBins.Neta
nPtbins  = recoBins.Npt

# consider only some reco bins to be scaled
###
# pt
ptBinToScale = []  # will store a bool to assess whether the given ipt index has to be scaled
hasPtRange = False
for bin in range(recoBins.Npt):
    ptBinToScale.append(False)
ptRanges = options.pt_range
if len(ptRanges):
    hasPtRange = True
    print "Only bins with reco pt in the following ranges will be scaled"
    print options.pt_range            
    for index in range(recoBins.Npt):
        for pair in ptRanges:
        #print pair
            if recoBins.ptBins[index] >= pair[0] and recoBins.ptBins[index+1] <= pair[1]:
                ptBinToScale[index] = True
else:
    hasPtRange = False
###
# eta
etaBinToScale = []  # will store a bool to assess whether the given ieta index has to be scaled
hasEtaRange = False
for bin in range(recoBins.Neta):
    etaBinToScale.append(False)
etaRanges = options.eta_range
if len(etaRanges):
    hasEtaRange = True
    print "Only bins with reco eta in the following ranges will be scaled"
    print options.eta_range            
    for index in range(recoBins.Neta):
        for pair in etaRanges:
        #print pair
            if recoBins.etaBins[index] >= pair[0] and recoBins.etaBins[index+1] <= pair[1]:
                etaBinToScale[index] = True
else:
    hasEtaRange = False
###
###

outfile = outdir + os.path.basename(infile).replace(".root","_SCALEXSEC.root")
if options.name:
    outfile = outdir + options.name

match = re.compile(options.match)
antimatch = re.compile(options.excludeMatch)

tfno = ROOT.TFile.Open(infile,"READ")
if not tfno or not tfno.IsOpen():
    raise RuntimeError('Unable to open file {fn}'.format(fn=infile))

of = ROOT.TFile(outfile,'recreate')
if not of or not of.IsOpen():
    raise RuntimeError('Unable to open file {fn}'.format(fn=outfile))
nKeys = tfno.GetNkeys()
nCopiedKeys = 0
nScaledKeys = 0
for ikey,e in enumerate(tfno.GetListOfKeys()):
    name = e.GetName()
    obj  = e.ReadObj()
    if not obj.InheritsFrom("TH1"): continue
    if not obj:
        raise RuntimeError('Unable to read object {n}'.format(n=name))
    newname = name
    newobj = obj.Clone(newname)
    # if re.match(options.match,name):
    if match.match(name) and not antimatch.match(name):
        if hasEtaRange or hasPtRange:
            # here we will scale only some bins
            for ib in range(newobj.GetNbinsX()):
                etabin = (ib % nEtabins)   # from 0 to nEtabins-1
                ptbin = int(ib / nEtabins) # from 0 to nPtBins-1
                if etaBinToScale[etabin] or ptBinToScale[ptbin]:
                    if options.scaleFromHist: 
                        scale = h2.GetBinContent(etabin+1,ptbin+1)
                        scale *= options.scale  # global scale, default is one so no harm
                    else: 
                        scale = option.scale
                    newobj.SetBinContent(ib+1,scale*newobj.GetBinContent(ib+1))
                    newobj.SetBinError(ib+1,scale*newobj.GetBinError(ib+1))
        elif options.scaleFromHist:
            newobj.Scale(options.scale)
            # for ib in range(1,newobj.GetNbinsX()+2):
            #     diff = abs(newobj.GetXaxis().GetBinLowEdge(ib)-newobj.GetXaxis().GetBinLowEdge(ib))
            #     if diff > 0.0:
            #         print "bin %d: lower edge difference = %.5f" % (ib,diff)
            #print "newobj has %d bins from %.1f to %.1f" % (newobj.GetNbinsX(),
            #                                                newobj.GetXaxis().GetBinLowEdge(1),
            #                                                newobj.GetXaxis().GetBinUpEdge(newobj.GetNbinsX())) 
            #if ikey > 1: quit()
            newobj.Multiply(h2_unrolledTo1D)            
            # this would be slower than Multiply I think
            # for ib in range(newobj.GetNbinsX()):
            #     etabin = (ib % nEtabins)   # from 0 to nEtabins-1
            #     ptbin = int(ib / nEtabins) # from 0 to nPtBins-1
            #     scale = h2.GetBinContent(etabin+1,ptbin+1)
            #     scale *= options.scale  # global scale, default is one so no harm
            #     newobj.SetBinContent(ib+1,scale*newobj.GetBinContent(ib+1))
            #     newobj.SetBinError(ib+1,scale*newobj.GetBinError(ib+1))
        else:
            newobj.Scale(options.scale)
        nScaledKeys += 1
    newobj.Write(newname)        
    nCopiedKeys += 1
    #sys.stdout.write('Key {num}/{tot}   \r'.format(num=ikey+1,tot=nKeys))
    #sys.stdout.flush()
    if (ikey+1) % 100 == 0:
        sys.stdout.write('Key {0:.2%}     \r'.format(float(ikey+1)/nKeys))
        sys.stdout.flush()
print "Copied {n}/{tot} from {fn} to {fo}".format(n=str(nCopiedKeys),tot=str(nKeys),fn=infile,fo=outfile)
print "Scaled {n}/{tot} from {fn} to {fo}".format(n=str(nScaledKeys),tot=str(nKeys),fn=infile,fo=outfile)
#of.Write()
of.Close()
tfno.Close()

if options.suppressRootWarnings:
    ROOT.gErrorIgnoreLevel = savErrorLevel
