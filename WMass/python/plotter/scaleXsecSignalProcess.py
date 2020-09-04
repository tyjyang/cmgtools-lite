#!/usr/bin/env python                                                       

# script to scale content of histograms for some pt bins in xsec_shapes.root. Needed on the combined root fileswhen in the electron channel some bins had different names than in the muon channel
 
#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np

from w_helicity_13TeV.make_diff_xsec_cards import get_ipt_from_process_name

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
#parser.add_option("-i", "--indir",    dest="indir",  type="string", default="", help="Input folder");
parser.add_option("-f", "--infile",   dest="infile", type="string", default="", help="Input file name");
parser.add_option("-o", "--outdir",   dest="outdir", type="string", default="./", help="Output folder (current one as default)");
parser.add_option("--pt-bins-scale-xsec", dest="ptBinsScaleXsec", type="string", default='0,1', help="Comma separated list of pt bins for which the content should be scaled");
parser.add_option("-m", "--match",   dest="match",   type="string", default="x_W.*_lep_.*_ipt_.*",    help="Regular expression for match");
parser.add_option("-s", "--scale",   dest="scale",   type="float", default=1.0, help="Scale factor");
#parser.add_option("-f", "--flavour", dest="flavour", type="string", default='el', help="Channel: either 'el' or 'mu'");
# charge is only used for the file name, not for the match
#parser.add_option("-c", "--charge",  dest="charge",  type="string", default='', help="Charge: either 'plus' or 'minus'");
(options, args) = parser.parse_args()

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
ptBinsMatch = [int(x) for x in options.ptBinsScaleXsec.split(",")]

#infile = options.indir + "/" + options.infile
infile = options.infile
outfile = outdir + os.path.basename(infile).replace(".root","_SCALEXSEC.root")

match = re.compile(options.match)

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
    newobj = obj.Clone(newname)
    # if re.match(options.match,name):
    if match.match(name):
        ipt = get_ipt_from_process_name(name)
        if ipt in ptBinsMatch: 
            newobj.Scale(options.scale)
    newobj.Write(newname)        
    nCopiedKeys += 1
    #sys.stdout.write('Key {num}/{tot}   \r'.format(num=ikey+1,tot=nKeys))
    #sys.stdout.flush()
    if (ikey+1) % 100 == 0:
        sys.stdout.write('Key {0:.2%}     \r'.format(float(ikey+1)/nKeys))
        sys.stdout.flush()
print "Copied {n}/{tot} from {fn} to {fo}".format(n=str(nCopiedKeys),tot=str(nKeys),fn=infile,fo=outfile)
#of.Write()
of.Close()
tfno.Close()

