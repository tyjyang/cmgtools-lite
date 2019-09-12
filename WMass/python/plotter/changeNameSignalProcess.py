#!/usr/bin/env python                                                       

# script to change name of histograms inside Wxx_shapes.root, for some pt bins. 
# This is needed to prepare for combination, and is used on the electron channel
 
#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np

from w_helicity_13TeV.make_diff_xsec_cards import get_ipt_from_process_name

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
parser.add_option("-i", "--indir",    dest="indir",  type="string", default="", help="Input folder");
parser.add_option("-o", "--outdir",   dest="outdir", type="string", default="./", help="Output folder (current one as default)");
parser.add_option("--pt-bins-name-change", dest="ptBinsNameChange", type="string", default='0,1', help="Comma separated list of pt bins for which the process name should be changed");
parser.add_option("-r", "--replace",   dest="replace", type="string", nargs=2,    help="pass two strings: first is match, second is replace");
parser.add_option("-f", "--flavour",   dest="flavour", type="string", default='el', help="Channel: either 'el' or 'mu'");
parser.add_option("-c", "--charge",    dest="charge",  type="string", default='', help="Charge: either 'plus' or 'minus'");
(options, args) = parser.parse_args()

if not len(options.indir):
    print "Warning: you must specify a folder with signal root files with option --indir"
    quit()

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
ptBinsMatch = [int(x) for x in options.ptBinsNameChange.split(",")]

infile = options.indir + "/W{fl}_{ch}_shapes.root".format(fl=flavour,ch=charge)
outfile = outdir + "W{fl}_{ch}_shapes_NEWPROCESSNAME.root".format(fl=flavour,ch=charge)

match = options.replace[0]
replace = options.replace[1]

tfno = ROOT.TFile.Open(infile,"READ")
if not tfno or not tfno.IsOpen():
    raise RuntimeError('Unable to open file {fn}'.format(fn=infile))

of = ROOT.TFile(outfile,'recreate')
if not of or not of.IsOpen():
    raise RuntimeError('Unable to open file {fn}'.format(fn=outfile))
nKeys = tfno.GetNkeys()
nCopiedKeys = 0
for ikey,e in enumerate(tfno.GetListOfKeys()):
    name = e.GetName()
    obj  = e.ReadObj()
    if not obj:
        raise RuntimeError('Unable to read object {n}'.format(n=name))
    newname = name
    if re.match("x_W{ch}.*_ipt_".format(ch=charge),name):
        ipt = get_ipt_from_process_name(name)
        if ipt in ptBinsMatch: 
            newname = name.replace(match,replace)
    newobj = obj.Clone(newname)
    newobj.Write(newname)        
    nCopiedKeys += 1
    sys.stdout.write('Key {num}/{tot}   \r'.format(num=ikey+1,tot=nKeys))
    sys.stdout.flush()
print "Copied {n}/{tot} from {fn}".format(n=str(nCopiedKeys),tot=str(nKeys),fn=infile)
of.Close()
tfno.Close()

