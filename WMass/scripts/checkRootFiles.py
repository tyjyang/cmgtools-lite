#!/usr/bin/env python

import sys,os,re
import ROOT

if __name__ == "__main__":
    
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] -d dir_with_rootfiles -f file_with_badfiles.txt")
    parser.add_option("-d", "--dir",      dest="fdir",    type="string", default='./', help="Directory with the ROOT files to be checked");
    parser.add_option("-f", "--filetxt",  dest="filetxt", type="string", default=None, help="If given, output the list of bad files in a txt file");

    (options, args) = parser.parse_args()
    
    fdir = options.fdir
    
    allfiles = os.listdir(fdir)
    rootfiles = filter(lambda x: x.endswith('.root'), allfiles)
    print "ROOT files to be checked in dir ",fdir," :"
    print rootfiles

    badrootfiles = []
    for f in rootfiles:
        tf = ROOT.TFile(f)
        good=True
        if not tf: good=False
        else:
            if tf.IsZombie() or tf.TestBit(ROOT.TFile.kRecovered): good=False
        print '{f} # {qual}'.format(f=f,qual="OK" if good else "zombie")
        if not good:
            badrootfiles.append(f)

    if options.filetxt:
        txt = open(options.filetxt,'w')
        for b in badrootfiles:
            txt.write("rm {b}\n".format(b=b))
        txt.close()
        print "list of bad files in ",options.filetxt

    print "DONE."
