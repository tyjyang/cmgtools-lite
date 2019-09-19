#!/bin/env python

#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np

from w_helicity_13TeV.make_diff_xsec_cards import get_ipt_from_process_name
from w_helicity_13TeV.make_diff_xsec_cards import get_ieta_from_process_name

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
#parser.add_option("-i", "--indir",    dest="indir",  type="string", default="", help="Input folder");
parser.add_option("-f", "--infile",   dest="infile", type="string", default="", help="Input file name");
parser.add_option("-o", "--outdir",   dest="outdir", type="string", default="./", help="Output folder (current one as default)");
parser.add_option("--pt-bins-add-group", dest="ptBinsAddGroup", type="string", default='0,1', help="Comma separated list of pt bins to add in group");
(options, args) = parser.parse_args()

# manage output folder
outdir = options.outdir
if not outdir.endswith('/'): outdir += "/"
if outdir != "./":
    if not os.path.exists(outdir):
        print "Creating folder", outdir
        os.system("mkdir -p " + outdir)

ptBinsMatch = [int(x) for x in options.ptBinsAddGroup.split(",")]

datacard = options.infile
newdatacard = outdir + os.path.basename(datacard).replace(".txt","_NEWGROUP.txt")

# now create new datacard with modified lines according to input options
#newdatacard = outdir + dcname.replace(dctag,"Wlep")
newdcfile = open(newdatacard,'w')

with open(datacard) as thisfile:
    for l in thisfile.readlines():        
        if re.match('.* sumGroup = .*',l):
            tokens = l.split(" sumGroup = ")
            groupname = tokens[0]
            pois = tokens[1]
            # look for charge in group name
            charge = ""
            if "plus" in groupname: 
                charge = "plus"
            elif "minus" in groupname: 
                charge = "minus"
            else:
                raise RuntimeError('Error: charge not found in group with name {gn}'.format(gn=groupname))
            # look for var in group name
            var = ""
            ivar = -1
            if "_ieta_" in groupname:
                var = "eta"
                ivar = get_ieta_from_process_name(groupname)
            elif "_ipt_" in groupname:
                var = "pt"
                ivar = get_ipt_from_process_name(groupname)
            else:
                raise RuntimeError('Error: ieta or ipt not found in group with name {gn}'.format(gn=groupname))
                
            addpois = ""
            if var == "eta":
                for ipt in ptBinsMatch:
                    poi = "W{ch}_el_ieta_{ieta}_ipt_{ipt} ".format(ch=charge, ieta=ivar, ipt=ipt)
                    addpois += poi
            elif var == "pt":
                if ivar in ptBinsMatch:
                    neta = len(pois.split(" "))
                    for ieta in range(neta):
                        poi = "W{ch}_el_ieta_{ieta}_ipt_{ipt} ".format(ch=charge, ieta=ieta, ipt=ivar)
                        addpois += poi
                else:
                    addpois = ""

            newdcfile.write(groupname + " sumGroup = " + addpois + " " + pois)        
        else:
            newdcfile.write(l)
newdcfile.write('\n') 
newdcfile.close()
print "-------------------------------------------"
print "Created new datacard: %s" % newdatacard
print "-------------------------------------------"


