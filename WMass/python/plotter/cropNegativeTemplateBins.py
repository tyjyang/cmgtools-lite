#!/usr/bin/env python                                                       
 
#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np
import root_numpy
import argparse

# in the fit, bins with 0 content in the nominal are ignored.
# so, after we made the shape files, we can just set to 0 the negative content of bins in the
# nominal templates (no need to touch the alternate as well)
# note that, by default, this function modify the original input file

def cropNegativeContent(h, silent=True, cropError=False):

    dim = h.GetDimension()
    nbins = 0
    if   dim == 1: nbins = h.GetNbinsX() + 2
    elif dim == 2: nbins = (h.GetNbinsX() + 2) * (h.GetNbinsY() + 2)
    elif dim == 3: nbins = (h.GetNbinsX() + 2) * (h.GetNbinsY() + 2) * (h.GetNbinsZ() + 2)

    integral = h.Integral()
    for i in range(nbins):
        nom  = h.GetBinContent(i)
        if nom<0.0: 
            h.SetBinContent(i, 0)
            if cropError:
                h.SetBinError(i, 0)
    integralNonNeg = h.Integral()
    if not silent:
        print("{n}: original integral = {i} changed by {r}".format(n=h.GetName(),
                                                                   i=str(integral),
                                                                   r=str(integralNonNeg/integral)))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--infile", type=str, default="", help="Input file name");
    parser.add_argument("-o", "--outdir", type=str, default="./", help="Output folder (current one as default)");
    parser.add_argument("-p", "--processes", type=str, default="x_Z|x_W.*|x_Top.*|.*DiBosons.*|x_data_fakes.*", help="Regular expression for histograms to be affected");
    parser.add_argument("-a", "--all", dest="allHists",   action="store_true", help="Crop all templates, not just those for systematics");
    parser.add_argument("-s", "--silent",  action="store_true", help="Print info when acting on histograms");
    args = parser.parse_args()
    options = args
    
    if not len(args.infile):
        print("Warning: you must specify an input file name with option --infile")
        quit()

    # manage output folder
    outdir = args.outdir
    if not outdir.endswith('/'): outdir += "/"
    if outdir != "./":
        if not os.path.exists(outdir):
            print("Creating folder", outdir)
            os.system("mkdir -p " + outdir)

    basenameFile = os.path.basename(args.infile)

    charge = ""
    flavour = ""
    if not any(ch in basenameFile for ch in ["plus", "minus"]):
        print("Warning: I could not understand charge plus|minus from file name. Abort")
        quit()
    else:
        charge = "plus" if "plus" in basenameFile else "minus"

    if not any(fl in basenameFile for fl in ["el_", "mu_"]):
        print("Warning: I could not understand flavour el|mu from file name. Abort")
        quit()
    else:
        flavour = "mu" if "mu_" in basenameFile else "el"

    isMu = True if flavour == "mu" else False

    infileCopy = outdir + basenameFile.replace(".root","_ORIGINAL.root")
    print("Copying original file with shapes into a backup")
    cpcmd = "cp {f} {fnew}".format(f=args.infile,fnew=infileCopy)
    print(cpcmd)
    os.system(cpcmd)
    outfile = outdir + basenameFile
    if os.path.abspath(outfile) != os.path.abspath(args.infile):
        print("Copying original file with shapes again in output folder")
        cpcmd = "cp {f} {fnew}".format(f=args.infile,fnew=outfile)
        print(cpcmd)
        os.system(cpcmd)
    else:
        print("Output folder is the same as the input file. I will update the input file")

    tfno = ROOT.TFile(args.infileCopy,'READ')
    if not tfno or not tfno.IsOpen():
        raise(RuntimeError('Unable to open file {fn}'.format(fn=args.infile)))
    of = ROOT.TFile(outfile,'UPDATE')
    if not of or not of.IsOpen():
        raise(RuntimeError('Unable to open file {fn}'.format(fn=outfile)))
    nKeys = tfno.GetNkeys()
    nCopiedKeys = 0
    for ikey,e in enumerate(tfno.GetListOfKeys()):
        name = e.GetName()
        if not re.match(args.processes,name): continue
        if not args.allHists and re.match(".*(Up|Down)",name): continue
        obj  = e.ReadObj()
        if not obj:
            raise(RuntimeError('Unable to read object {n}'.format(n=name)))
        cropNegativeContent(obj, silent=args.silent)
        obj.Write(name,ROOT.TObject.kOverwrite)
        nCopiedKeys += 1
        if (ikey+1) % 100 == 0:
            sys.stdout.write('Key {0:.2%}     \r'.format(float(ikey+1)/nKeys))
            sys.stdout.flush()
   
    print("Overwrote {n}/{tot} from {fn}".format(n=str(nCopiedKeys),tot=str(nKeys),fn=args.infile))
    print("Updated file saved in {of}".format(of=outfile))
    of.Close()
    tfno.Close()

