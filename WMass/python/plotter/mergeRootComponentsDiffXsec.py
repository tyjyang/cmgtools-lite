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
# merge signal shapes with Z and data+backgrounds
# first, merge shapes with Z and remove data from there (data is just the MC)
# then, merge signal, data and merged Z

# Z and data are assumed to be in folder given with --indir-bkg
# signal is assumed to be named W<flavour>_<charge>_shapes_signal.root

## python mergeRootComponentsDiffXsec.py -f mu -c minus --indir-bkg  cards/diffXsec_mu_2018_11_24_group10_onlyBkg/part0/  --indir-sig cards/diffXsec_mu_2018_11_24_group10_onlyBkg/ -o cards/diffXsec_mu_2018_11_24_group10_onlyBkg/ -d

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
parser.add_option("-d", "--dry-run", dest="dryrun",   action="store_true", default=False, help="Dry run: print commands but do not merge");
parser.add_option(      "--indir-bkg", dest="indirBkg", type="string", default="", help="Input folder with root files for Z, data and other backgrounds");
parser.add_option(      "--indir-sig", dest="indirSig", type="string", default="", help="Input folder with root files for signal");
parser.add_option("-o", "--outdir",    dest="outdir", type="string", default="./", help="Output folder (current one as default)");
parser.add_option("-n", "--name",      dest="name",   type="string", default="", help="Name for output root file (if not given, name is W<flavou>_<charge>_shapes.root ).");
parser.add_option("-s", "--suffix",    dest="suffix", type="string", default="", help="Suffix to add to output file before extension. Ineffective if using option -n < name>");
parser.add_option("-f", "--flavour",   dest="flavour", type="string", default='', help="Channel: either 'el' or 'mu'");
parser.add_option("-c", "--charge",    dest="charge", type="string", default='', help="Charge: either 'plus' or 'minus'");
(options, args) = parser.parse_args()
    
# manage output folder
outdir = options.outdir
if not outdir.endswith('/'): outdir += "/"
if outdir != "./":
    if not os.path.exists(outdir):
        print "Creating folder", outdir
        os.system("mkdir -p " + outdir)

if options.flavour not in ["el", "mu"]:
    print "Warning: you must specify a lepton flavour with option -f el|mu"
    quit()
if options.charge not in ["plus", "minus"]:
    print "Warning: you must specify a charge with option -c plus|minus"
    quit()

charge = options.charge
flavour = options.flavour
#binname = options.bin if len(options.bin) else "W%s" % flavour

if not len(options.indirBkg):
    print "Warning: you must specify a folder with data and background root files with option --indir-bkg"
    quit()
if not len(options.indirSig):
    print "Warning: you must specify a folder with signal root files with option --indir-sig"
    quit()

# define ultimate output file
shapename = ""
if options.name == "":
    shapename = "{od}W{fl}_{ch}_shapes.root".format(od=outdir, fl=flavour, ch=charge)
    if len(options.suffix): shapename = shapename.replace(".root","_{sf}.root".format(sf=options.suffix))
else:
    shapename = outdir + options.name

print ""
print "-"*20
print ""

## prepare the relevant files. First merge Z with correct charge
zMatch = "^Z_{fl}_{ch}.*".format(fl=flavour,ch=charge)
zfiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(options.indirBkg) for f in fn if (f.endswith('.input.root') and re.match(zMatch,f))]

#for f in zfiles:
#    print f

# merge Z files (will have to remove x_data_obs later)
tmpZfile="{od}Z_{fl}_{ch}_mergeTMP.root".format(od=outdir, fl=flavour, ch=charge)
cmdMerge = "hadd -f -k -O {tmp} {zfs}".format(tmp=tmpZfile, zfs=" ".join([f for f in zfiles]))
#print cmdMerge
Zfile="{od}Z_{fl}_{ch}_merge.root".format(od=outdir, fl=flavour, ch=charge)

print "Merging Z"
if not options.dryrun: os.system(cmdMerge)

print "Remove x_data_obs from Z, and replace 'x_Z_dy_' with 'x_Z_'"
print "Also changing Dn to Down"
nZcopied = 0
# open new Z file to remove data from input file
#----------------------------------
newname = ""
zpdf = []
znominal = None
nMirroredPDF = 0
if not options.dryrun:
    tf = ROOT.TFile.Open(tmpZfile,"READ")
    if not tf or not tf.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=tmpZfile))

    # open output file
    of = ROOT.TFile(Zfile,'recreate')
    if not of or not of.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=Zfile))

    for ikey,e in enumerate(tf.GetListOfKeys()):
        name = e.GetName()
        obj  = e.ReadObj()
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        if "data_obs" in name: continue
        if "x_Z_dy_" in name: newname = name.replace("x_Z_dy_","x_Z_")
        else: newname = name
        if name == "x_Z": znominal = obj.Clone(name)
        if "_pdf" in name:
            zpdf.append(obj.Clone(newname))
            continue
        if newname.endswith("Dn"): newname = newname.replace("Dn","Down")
        newobj = obj.Clone(newname)
        newobj.Write(newname)
        nZcopied += 1

    for h in zpdf:
        nMirroredPDF += 1
        (alternate,mirror) = mirrorShape(znominal,h,h.GetName())
        for alt in [alternate,mirror]:
            alt.Write()

    of.Close()
    tf.Close()
#----------------------------------

print "Copied {n} histograms in {zf} (x_data_obs removed)".format(n=str(nZcopied),zf=Zfile)
print "Created {n} mirrored histograms for PDFs in {zf} ".format(n=str(nMirroredPDF),zf=Zfile)
print "Removing temporary file {tmp}".format(tmp=tmpZfile)
if not options.dryrun: os.system("rm {tmp}".format(tmp=tmpZfile))    
print "-"*20
print ""

dataAndBkgFile = "{obkg}bkg_and_data_{fl}_{ch}.input.root".format(obkg=options.indirBkg, fl=flavour, ch=charge)
dataAndBkgFileTmp = dataAndBkgFile.replace(".input.root","TMP.input.root")
print "Creating temporary file {bkg} to remove 'x_data' histogram".format(bkg=dataAndBkgFileTmp)
print "Also changing Dn to Down"
# now remove x_data from bkg_and_data_*.input.root (we only need x_data_obs, this x_data should be removed)
#----------------------------------
if not options.dryrun:
    tf = ROOT.TFile.Open(dataAndBkgFile,"READ")
    if not tf or not tf.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=dataAndBkgFile))

    # open output file
    of = ROOT.TFile(dataAndBkgFileTmp,'recreate')
    if not of or not of.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=dataAndBkgFileTmp))
    
    nKeys = tf.GetNkeys()
    nCopiedKeys = 0
    for ikey,e in enumerate(tf.GetListOfKeys()):
        name = e.GetName()
        obj  = e.ReadObj()
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        if name == "x_data": continue
        newname = name
        if newname.endswith("Dn"): 
            newname = newname.replace("Dn","Down")        
        newobj = obj.Clone(newname)
        newobj.Write(newname)
        nCopiedKeys += 1

    print "Copied {n}/{tot} from {bkg}".format(n=str(nCopiedKeys),tot=str(nKeys),bkg=dataAndBkgFile)
    of.Close()
    tf.Close()
#----------------------------------
print "-"*20
print ""


print "Now merging signal + Z + data + other backgrounds"
sigfile = "{osig}W{fl}_{ch}_shapes_signal.root".format(osig=options.indirSig, fl=flavour, ch=charge)

cmdFinalMerge="hadd -f -k -O {of} {sig} {zf} {bkg}".format(of=shapename, sig=sigfile, zf=Zfile, bkg=dataAndBkgFileTmp)
print "Final merging ..."
print cmdFinalMerge
if not options.dryrun: os.system(cmdFinalMerge)

print "-"*20
print ""
print "Removing temporary file {bkg}".format(bkg=dataAndBkgFileTmp)
if not options.dryrun: os.system("rm {bkg}".format(bkg=dataAndBkgFileTmp))

print "-"*20
print ""
print ""
print "Wrote root file in %s" % shapename
print ""
