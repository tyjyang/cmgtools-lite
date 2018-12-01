#!/usr/bin/env python                                                                                                                                                        
#from shutil import copyfile
import re, sys, os, os.path, ROOT
from array import array
#import numpy as np

from w_helicity_13TeV.make_diff_xsec_cards import getXYBinsFromGlobalBin
from w_helicity_13TeV.make_diff_xsec_cards import getArrayParsingString
from w_helicity_13TeV.make_diff_xsec_cards import getDiffXsecBinning
from w_helicity_13TeV.make_diff_xsec_cards import templateBinning

from w_helicity_13TeV.mergeCardComponentsAbsY import mirrorShape

# originally developed to make cards for diff xsec in lepton pt and |eta|
# transform TH3 into TH1 usable by combine
# this is supposed to be done only on signal histograms, then they will be merged to data and background ones by hand

## python makeTH1FromTH3.py cards/diffXsec_mu_2018_11_24_group10_onlyBkg/wmass_varhists_mu_withScale.root -o cards/diffXsec_mu_2018_11_24_group10_onlyBkg/ -f mu -c plus --binfile cards/diffXsec_mu_2018_11_24_group10_onlyBkg/binningPtEta.txt

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] shapes.root")
parser.add_option("-o", "--outdir",    dest="outdir", type="string", default="./", help="Output folder (current one as default)");
parser.add_option("-n", "--name",      dest="name",   type="string", default="", help="Name for output root file (if not given, name is W<flavou>_<charge>_shapes_signal.root ).");
parser.add_option("-s", "--suffix",    dest="suffix", type="string", default="", help="Suffix to add to output file before extension. Ineffective if using option -n < name>");
parser.add_option("-f", "--flavour",   dest="flavour", type="string", default='', help="Channel: either 'el' or 'mu'");
#parser.add_option("-g", "--group",     dest="group", type="int", default='0', help="In case signal bins for diff.xsec were grouped, the number of grouped bins appeara in the name of the process as '<name>_group<n>'");
parser.add_option("-c", "--charge",    dest="charge", type="string", default='', help="Charge: either 'plus' or 'minus'");
parser.add_option(      '--binfile'  , dest='binfile', default='binningPtEta.txt', type='string', help='eta-pt binning for templates.')
parser.add_option(      "--effStat-all", dest="effStatAll",   action="store_true", default=False, help="If True, assign any EffStat syst to any eta bin: otherwise, it is associated only to the corresponding eta bin");
parser.add_option(      "--symmetrize-syst",  dest="symSyst", type="string", default='.*EffStat.*|.*pdf.*', help="Regular expression matching systematics whose histograms should be symmetrized with respect to nominal (Up and Down variations not already present)");
parser.add_option(      "--fixmuRmuFNUp", dest="fixmuRmuFNUp",   action="store_true", default=False, help="Make a fix: muRmuXXFUp should become muRmuFXXUp");
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
etaPtBinningVec = getDiffXsecBinning(etaPtBinningFile, "gen")
genBins  = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])
#
print ""
recoBins.printBinAll()
genBins.printBinAll()
print ""

# creating binning to define TH1D lated
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
effstatOffset = 25 if (genBins.etaBins[-1] < 2.45) else 26  # if gen eta goes up to 2.4, set offset to 25, it is used below

print "Going to unroll TH3 into TH1 ..."
nKeys = tf.GetNkeys()

pattBug = re.compile('(muRmu)(\d+)(FUp)')
pattEffStat = re.compile('(EffStat)(\d+)')

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
    signalRoot = "W{ch}_{fl}".format(ch=charge,fl=flavour)

    lasttoken = name.split("_")[-1]
    toBeMirrored = False
    isNominal = False
    isEffStat = False
    etaEffStat = -1
    if lasttoken == "globalBin": isNominal = True
    elif re.match(symSystMatch,lasttoken): toBeMirrored = True

    if options.fixmuRmuFNUp and "muRmu" in lasttoken:
        tkn = pattBug.findall(lasttoken)
        if len(tkn):
            wrong = "muRmu" + str(tkn[0][1]) + "FUp"
            correct = "muRmuF" + str(tkn[0][1]) + "Up"
            lasttoken = lasttoken.replace(wrong,correct)

    # if needed, select the gen eta bin that will be affected by effstat
    if not options.effStatAll and re.match('.*EffStat.*',lasttoken) :
        isEffStat = True
        tkn = pattEffStat.findall(lasttoken)  # if lasttoken == "EffStatYY", tkn is ('EffStat', 'YY')
        etaEffStat = int(tkn[0][1])
        etaEffStat = etaEffStat - effstatOffset # assess whether it is in the first half, corresponding to negative eta
        # since ieta goes from 0 to XX, make etaEffStat goes from 0 to XX
        if etaEffStat < 0: 
            etaEffStat = abs(etaEffStat) - 1
        #else:
        #    etaEffStat = etaEffStat        

    # now loop on bins along Z, each is a signal template, bin 0 (underflow) has the outliers (this convention might change, please check)
    for iz in range(obj.GetNbinsZ()+1):         
        globalbin = iz
        if globalbin > 0:
            ieta,ipt = getXYBinsFromGlobalBin(globalbin-1,genBins.Neta,binFrom0=True)
            if isEffStat and etaEffStat != ieta: continue  # if this is EffStat and the ieta bin does not correspond to the effstat continue
            newname = "x_" + signalRoot + "_ieta_{ie}_ipt_{ip}_".format(ie=str(ieta),ip=str(ipt)) + signalRoot
        else:
            newname = "x_" + signalRoot + "_outliers_" + signalRoot
        if not isNominal:            
            if "lepScale" in lasttoken: 
                newname += "_" + lasttoken.replace("lepScale","CMS_Wmu_muscale" if flavour == "mu" else "CMS_We_elescale")
            elif "lepEff" in lasttoken:
                newname += "_" + lasttoken.replace("lepEff","CMS_Wmu_sig_lepeff" if flavour == "mu" else "CMS_We_sig_lepeff")
            else:
                newname += "_" + lasttoken
            if isEffStat and "ErfPar" not in lasttoken: newname = newname.replace("Erf","ErfPar")  # patch for bad names inside input file (ErfXX instead of ErfParXX)
            if newname.endswith('Dn'): newname = newname[:-2] + "Down"
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
    (alternate,mirror) = mirrorShape(nominalTH1[pfx],symmetrizedTH1[key],key)
    for alt in [alternate,mirror]:
        alt.Write()    

of.Close()

tf.Close()

    
print ""
print "Wrote root file in %s" % shapename
print ""
