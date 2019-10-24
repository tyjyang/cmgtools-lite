#!/usr/bin/env python                                                       

# script to change name of histograms inside Wxx_shapes.root
 
#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np

from w_helicity_13TeV.make_diff_xsec_cards import get_ipt_from_process_name

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
parser.add_option("-f", "--infile",   dest="infile",  type="string", default="", help="Input folder");
parser.add_option("-o", "--outdir",   dest="outdir", type="string", default="./", help="Output folder (current one as default)");
parser.add_option("-r", "--replace",   dest="replace", action="append", type="string", nargs=2, default=[],   help="pass two strings: first is match, second is replace. Can specify multiple times");
parser.add_option("-x", "--regexp",   dest="regexp", type="string", default=".*", help="Regexp to match only some keys");
parser.add_option("-v", "--verbose", dest="verbose" , default=False, action="store_true", help="Print what is being changed")
parser.add_option("--save-only-changed-histo", dest="saveOnlyChangedHisto" , default=False, action="store_true", help="If true, just save histo whose name was changed (output file can then be merged to input using hadd). This is faster than looping and copying all keys")
(options, args) = parser.parse_args()

if not len(options.infile):
    print "Warning: you must specify a folder with signal root files with option -f/--infile [arg]"
    quit()

# manage output folder
outdir = options.outdir
if not outdir.endswith('/'): outdir += "/"
if outdir != "./":
    if not os.path.exists(outdir):
        print "Creating folder", outdir
        os.system("mkdir -p " + outdir)

infile = options.infile
if not infile.endswith(".root"):
    print "Warning: input file does not end with '.root'. Is it a root file? Abort"
    quit()

outfile = outdir + os.path.basename(infile).replace(".root","_NEWPROCESSNAME.root")
#"W{fl}_{ch}_shapes_NEWPROCESSNAME.root".format(fl=flavour,ch=charge)

regexp = options.regexp
match_regexp = ""
if len(regexp):
    match_regexp = re.compile(regexp)

matchReplace = {}
if len(options.replace):
    for i in range(len(options.replace)):
        matchReplace[options.replace[i][0]] = options.replace[i][1]
#match = options.replace[0]
#replace = options.replace[1]

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
    # does it read to memory of file?    
    #obj  = e.ReadObj()
    of.cd()
    obj = None
    newname = name
    nameHasChanged = False
    if len(regexp) and match_regexp.match(name):        
        for key in matchReplace:
            if key in newname:
                newname = newname.replace(key,matchReplace[key])
                nameHasChanged = True
        if options.verbose and nameHasChanged:
            print "Changing {i} to {f}".format(i=name, f=newname)
    if not options.saveOnlyChangedHisto or nameHasChanged:
        obj = tfno.Get(name)    
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        newobj = obj.Clone(newname)
        newobj.Write(newname)        
        #obj.Write(newname)        
        nCopiedKeys += 1
    if (ikey+1) % 100 == 0:
        sys.stdout.write('Key {num}/{tot}   \r'.format(num=ikey+1,tot=nKeys))
        sys.stdout.flush()
print "Copied {n}/{tot} \nfrom {fn} \nto {of}".format(n=str(nCopiedKeys),tot=str(nKeys),
                                                     fn=infile,
                                                     of=outfile)
of.Close()
tfno.Close()

