#!/usr/bin/env python                                                                                                                                                        
#from shutil import copyfile
import re, sys, os, os.path, ROOT
from array import array
#import numpy as np

from w_helicity_13TeV.make_diff_xsec_cards import getXYBinsFromGlobalBin
from w_helicity_13TeV.make_diff_xsec_cards import getArrayParsingString
from w_helicity_13TeV.make_diff_xsec_cards import getDiffXsecBinning
from w_helicity_13TeV.make_diff_xsec_cards import templateBinning
from w_helicity_13TeV.make_diff_xsec_cards import getArrayBinNumberFromValue

from w_helicity_13TeV.mergeCardComponentsAbsY import mirrorShape

# originally developed to make cards for diff xsec in lepton pt and |eta|
# transform TH3 into TH1 usable by combine
# this is supposed to be done only on signal histograms, then they will be merged to data and background ones by hand

## python makeTH1FromTH3.py cards/diffXsec_mu_2018_11_24_group10_onlyBkg/wmass_varhists_mu_withScale.root -o cards/diffXsec_mu_2018_11_24_group10_onlyBkg/ -f mu -c plus --binfile cards/diffXsec_mu_2018_11_24_group10_onlyBkg/binningPtEta.txt

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] shapes.root")
parser.add_option("-o", "--outdir",    dest="outdir", type="string", default="./", help="Output folder (current one as default)");
parser.add_option("-n", "--name",      dest="name",   type="string", default="", help="Name for output root file (if not given, name is W<flavour>_<charge>_shapes_signal.root ).");
parser.add_option("-s", "--suffix",    dest="suffix", type="string", default="", help="Suffix to add to output file before extension. Ineffective if using option -n < name>");
parser.add_option("-f", "--flavour",   dest="flavour", type="string", default='', help="Channel: either 'el' or 'mu'");
#parser.add_option("-g", "--group",     dest="group", type="int", default='0', help="In case signal bins for diff.xsec were grouped, the number of grouped bins appeara in the name of the process as '<name>_group<n>'");
parser.add_option("-c", "--charge",    dest="charge", type="string", default='', help="Charge: either 'plus' or 'minus'");
parser.add_option(      '--binfile'  , dest='binfile', default='binningPtEta.txt', type='string', help='eta-pt binning for templates.')
parser.add_option(      "--effStat-all", dest="effStatAll",   action="store_true", default=False, help="If True, assign any EffStat syst to any eta bin: otherwise, it is associated only to the corresponding eta bin");
parser.add_option(      "--symmetrize-syst",  dest="symSyst", type="string", default='.*ErfPar.*EffStat.*|.*pdf.*|.*fsr.*', help="Regular expression matching systematics whose histograms should be symmetrized with respect to nominal (Up and Down variations not already present)");
parser.add_option(      "--shape-only-symmetrized-syst",  dest="shapeOnlySymSyst", type="string", default='.*fsr.*', help="Regular expression matching systematics that are mirrored and which shouldbe shap-only");
(options, args) = parser.parse_args()

if len(sys.argv) < 1:
    parser.print_usage()
    quit()
    
# manage output folder
outdir = options.outdir
if not outdir.endswith('/'): outdir += "/"
if outdir != "./":
    if not os.path.exists(outdir):
        print "Creating folder", outdir
        os.system("mkdir -p " + outdir)

# if not options.group:
#     print "Warning: you must specify a number of bins grouped together with option -g <ngroup>"
#     quit()
if options.flavour not in ["el", "mu"]:
    print "Warning: you must specify a lepton flavour with option -f el|mu"
    quit()
if options.charge not in ["plus", "minus"]:
    print "Warning: you must specify a charge with option -c plus|minus"
    quit()


charge = options.charge
flavour = options.flavour
#binname = options.bin if len(options.bin) else "W%s" % flavour
signalMatch = "W%s" % charge
symSystMatch = options.symSyst

# input histograms named as h3_minus_eta_pt_globalBin_Erf0EffStat41
# eta - pt - global gen bin on x-y-z axis
# in principle we only need the reco bins for the unrolling of each TH2 slice into a TH1
etaPtBinningFile = options.binfile
# get eta-pt binning for both reco and gen
etaPtBinningVec = getDiffXsecBinning(etaPtBinningFile, "reco")
recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])
etaPtGenBinningVec = getDiffXsecBinning(etaPtBinningFile, "gen")
genBins  = templateBinning(etaPtGenBinningVec[0],etaPtGenBinningVec[1])
#
print ""
recoBins.printBinAll()
genBins.printBinAll()
print ""

# creating binning to define TH1D later
thd1binning = [(0.5 + float(i)) for i in range(1+recoBins.NTotBins)]
#print thd1binning

# signal processes are like
# h3_minus_eta_pt_globalBin_Erf1EffStat29
# must become either (examples)
# x_Wplus_mu_outliers_Wplus_mu_group_46_pdf36Up
# or
# x_Wplus_mu_ieta_1_ipt_4_Wplus_mu_group_9_pdf57Up
# group_<n> is actually not needed, it is a legacy of the time when signal cards where made with makeShapeCards.py grouping signal bins to reduce the number of jobs

shapename = ""
if options.name == "":
    shapename = "{od}W{fl}_{ch}_shapes_signal.root".format(od=outdir, fl=flavour, ch=charge)
    if len(options.suffix): shapename = shapename.replace(".root","_{sf}.root".format(sf=options.suffix))
else:
    shapename = outdir + options.name

# open input file
tf = ROOT.TFile.Open(args[0],"READ")
if not tf or not tf.IsOpen():
    raise RuntimeError('Unable to open file {fn}'.format(fn=args[0]))

# open output file
of = ROOT.TFile(shapename,'recreate')
if not of or not of.IsOpen():
    raise RuntimeError('Unable to open file {fn}'.format(fn=shapename))

newname = ""
name = ""
obj = None
nominalTH1 = {}
symmetrizedTH1 = {}
# efficiencies are made with 0.1 eta bins between -2.5 and 2.5, even though the template might have a coarser binning
# if genEta ends at 2.4, only the nuisances EffStatXX with X from 2 to 49 are used
effstatOffset = 25 if flavour == "mu" else 26  
#effstatOffset = int(1 + 10*(genBins.etaBins[-1] + 0.001))

print "Going to unroll TH3 into TH1 ..."
nKeys = tf.GetNkeys()

pattEffStat = re.compile('(EffStat)(\d+)')
ietaTemplate = -1

# in order to prepare for flavor combination, use lep instead of mu|el in name
#signalRoot = "W{ch}_{fl}".format(ch=charge,fl=flavour)
signalRoot = "W{ch}_lep".format(ch=charge)

for ikey,e in enumerate(tf.GetListOfKeys()):

    # print status
    sys.stdout.write('Key {num}/{tot}   \r'.format(num=ikey+1,tot=nKeys))
    sys.stdout.flush()

    name = e.GetName()
    obj  = e.ReadObj()
    if not obj:
        raise RuntimeError('Unable to read object {n}'.format(n=name))
    if not obj.InheritsFrom("TH3"): continue
    if charge not in name: continue
    if "_globalBin" not in name: continue

    lasttoken = name.split("_")[-1]
    toBeMirrored = False
    isNominal = False
    isEffStat = False
    etaEffStat = -1
    if lasttoken == "globalBin": isNominal = True
    elif re.match(symSystMatch,lasttoken): toBeMirrored = True

    # if needed, select the gen eta bin that will be affected by effstat
    if not options.effStatAll and re.match('.*EffStat.*',lasttoken) :
        isEffStat = True
        tkn = pattEffStat.findall(lasttoken)  # if lasttoken == "EffStatYY", tkn is ('EffStat', 'YY')
        etaEffStat = int(tkn[0][1])
        etaEffStat = etaEffStat - effstatOffset # assess whether it is in the first half, corresponding to negative eta
        # since ieta goes from 0 to XX, make etaEffStat goes from 0 to XX
        if etaEffStat < 0: 
            etaEffStat = abs(etaEffStat) - 1
        # we always have 48 or 50 efficiency bins, but the reco binning might be coarser: get the eta for the efficiency bin
        etaBinCenter = etaEffStat * 0.1 + 0.05  
        ietaTemplate = getArrayBinNumberFromValue(genBins.etaBins,etaBinCenter)        
        # if the reco binning along eta is narrower than the EffStat histogram used for the reweighting, skip this etaEffStat
        # this happens for example for bin 1 and 50 in the electron channel if the reco binning (and therefore also the gen) stops at |eta|=2.4
        if etaBinCenter < recoBins.etaBins[0] or etaBinCenter > recoBins.etaBins[-1]:
            continue

    # now loop on bins along Z, each is a signal template, bin 0 (underflow) has the outliers (this convention might change, please check)
    for iz in range(obj.GetNbinsZ()+1):         
        globalbin = iz
        if globalbin > 0:
            ieta,ipt = getXYBinsFromGlobalBin(globalbin-1,genBins.Neta,binFrom0=True)
            if isEffStat and ietaTemplate != ieta: continue  # if this is EffStat and the ieta bin does not correspond to the effstat continue
            newname = "x_" + signalRoot + "_ieta_{ie}_ipt_{ip}".format(ie=str(ieta),ip=str(ipt))
        else:
            newname = "x_" + signalRoot + "_outliers"
        if not isNominal:            
            if "lepScale" in lasttoken: 
                newname += "_" + lasttoken.replace("lepScale","CMS_Wmu_muscale" if flavour == "mu" else "CMS_We_elescale")
            elif "lepUncorrScale" in lasttoken:
                # distinct from previous one because there is an integer number after lepUncorrScale
                newname += "_" + lasttoken.replace("lepUncorrScale","CMS_Wmu_muscale" if flavour == "mu" else "CMS_We_elescale")
            elif "lepEff" in lasttoken:
                newname += "_" + lasttoken.replace("lepEff","CMS_Wmu_sig_lepeff" if flavour == "mu" else "CMS_We_sig_lepeff")
            else:
                newname += "_" + lasttoken
            #if isEffStat and "ErfPar" not in lasttoken: newname = newname.replace("Erf","ErfPar")  # patch for bad names inside input file (ErfXX instead of ErfParXX) # no longer needed
            if newname.endswith('Dn'): newname = newname[:-2] + "Down"
            # for EffStat, add <flavour><charge> at the end, because we use independent systs for charge and flavour, so the names must be different
            # note that at this point the histograms don't have Up/Down in their name
            if isEffStat:
                suffixToAdd = "{fl}{ch}".format(fl=flavour, ch=charge)
                if "BinUncEffStat" in newname:
                    # add suffix before Up/Down
                    if newname.endswith("Up"): newname = newname[:-2] + suffixToAdd + "Up"
                    elif newname.endswith("Down"): newname = newname[:-4] + suffixToAdd + "Down"
                else:
                    newname += suffixToAdd
        hist = ROOT.TH1D(newname,"", len(thd1binning)-1, array('d', thd1binning))
        # unroll a slice of TH3 at fixed z (this is a TH2) into a 1D histogram
        for ix in range(1,1+obj.GetNbinsX()):
            for iy in range(1,1+obj.GetNbinsY()):
                bin = ix + obj.GetNbinsX() * (iy-1)
                hist.SetBinContent(bin, obj.GetBinContent(ix,iy,iz))
                hist.SetBinError(bin, obj.GetBinError(ix,iy,iz))
        if toBeMirrored:
            # save but not write it to file yet
            if newname not in symmetrizedTH1: symmetrizedTH1[newname] = hist.Clone()
        else:
            if isNominal:
                # save nominal to allow to retrieve it easily later: it will be used to symmetrize PDF and EffStat histograms
                if newname not in nominalTH1: nominalTH1[newname] = hist.Clone()
            hist.Write(newname)  # write in output file

# print "Printing key of symmetrizedTH1"
# for key in symmetrizedTH1: print key,
# print ""
# print "Printing key of nominalTH1"
# for key in nominalTH1: print key,
# print ""

print "Creating mirror histogram for PDFs and EffStat"     
# now symmetrize what need to be
nKeys = len(symmetrizedTH1)
for i,key in enumerate(symmetrizedTH1):

    sys.stdout.write('Key {num}/{tot}   \r'.format(num=i+1,tot=nKeys))
    sys.stdout.flush()

    # get name prefix to retrieve nominal
    pfx = "_".join(key.split('_')[:-1])  # remove last piece (might be 'pdf10')
    if re.match(options.shapeOnlySymSyst,key):        
        (alternate,mirror) = mirrorShape(nominalTH1[pfx],symmetrizedTH1[key],key,alternateShapeOnly=True,use2xNomiIfAltIsZero=True)
    else:
        (alternate,mirror) = mirrorShape(nominalTH1[pfx],symmetrizedTH1[key],key,use2xNomiIfAltIsZero=True)
    for alt in [alternate,mirror]:
        alt.Write()    

of.Close()

tf.Close()

    
print ""
print "Wrote root file in %s" % shapename
print ""
