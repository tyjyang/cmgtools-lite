#!/usr/bin/env python  

from shutil import copyfile
import re, sys, os, os.path, subprocess, json, ROOT
import numpy as np

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] infile.root")
parser.add_option('-c', '--class',          dest='className',      default='', type='string', help='Filter output based on the object class (e.g. -c TH1D)')
parser.add_option('-x', '--exclude-regexp', dest='exclude_regexp', default='', type='string', help='Like option -r, but will reject object if name matches')
parser.add_option('-f', '--output-format',  dest='output_format',  default='all', type='string', help='Print type, name and (for histograms only) integral (all, default); just name (name)')
parser.add_option('-s', '--silent'  ,       dest='silent',         default=False, action='store_true', help='Silent mode: just print summary, not all entries')
parser.add_option('-i', '--inherith',       dest='inherit',        default='', type='string', help='Filter output based on whether the object class inheriths from this (e.g. -i TH1 will match all root histogram objects)')  
parser.add_option('-r', '--regexp',         dest='regexp',         default='', type='string', help='Filter output passing a regular expression to be matched in the object name')
parser.add_option(      '--print-min-negative', dest='printMinNegative', default=False, action='store_true', help='If true, print also minimum bin content when negative (only for histograms)')
parser.add_option(      '--print-integral-cap-negative', dest='printIntegralCapNegative', default=False, action='store_true', help='If true, print integral when capping negative bin to 0 (print also ratio with actual integral')
parser.add_option(      '--print-axis-label-histogram', dest='printAxisLabelHisto', default='', type='string', help='If given, print axis labels of histogram with this name (everything else in the file will be skipped, and other options might be ignored). For a TH2, the x axis is chosen')
parser.add_option(      '--print-branch-tree', dest='printBranchTree', default='', type='string', help='If given, print names of branches in tree with this name (everything else in the file will be skipped, and other options might be ignored).')

(options, args) = parser.parse_args()

if len(sys.argv) < 1:
    parser.print_usage()
    quit()

tf = ROOT.TFile.Open(args[0],"READ")
nRead = 0
nNullIntegral = 0
nNegativeBin = 0

if not options.silent:
    print "-"*50
    if options.output_format == "all":    
        print "{: <10} {: <20}    {y}    {txt} ".format("Class","name",
                                                        y="integral(TH only)",
                                                        txt="integral(nonNeg)/integral" if options.printIntegralCapNegative else "")
    elif options.output_format == "name":
        print "name"
    print "-"*50
        
for k in tf.GetListOfKeys() :
    name=k.GetName()
    if options.printAxisLabelHisto != "" and name != options.printAxisLabelHisto:
        continue
    if options.printBranchTree != "" and name != options.printBranchTree:
        continue

    if len(options.regexp) and not re.match(options.regexp,name): continue
    if len(options.exclude_regexp) and re.match(options.exclude_regexp,name): continue
    obj=k.ReadObj()
    if len(options.className) and obj.ClassName() != options.className: continue
    if len(options.inherit) and not obj.InheritsFrom(options.inherit): continue
    nRead += 1
    integral = -1
    integralOnlyNonNegativeBin = 0.0
    minBinContent = 1
    if (obj.ClassName().startswith("TH") and obj.InheritsFrom("TH1")):

        if options.printAxisLabelHisto != "" and name == options.printAxisLabelHisto:
            print "-"*30
            print "Labels of histogram %s" % name
            print "-"*30
            for ib in range(1,1+obj.GetNbinsX()):
                print "{: <10} {l}".format(ib,l=obj.GetXaxis().GetBinLabel(ib))
            print "-"*30

        integral = obj.Integral() 
        minBinContent = obj.GetBinContent(obj.GetMinimumBin())
        if options.printIntegralCapNegative:
            for ibin in range(1,obj.GetNbinsX()+1):
                integralOnlyNonNegativeBin += max(0.0,obj.GetBinContent(ibin))

    if integral == 0.0: nNullIntegral += 1
    if minBinContent < 0: nNegativeBin += 1
    if not options.silent:
        if options.output_format == "all":
            #print "%s   %s   %s" % (obj.ClassName(), name, str(integral) if integral >= 0 else "")
            print "{: <10} {: <20}    {y}    {r}    {m} ".format(obj.ClassName(), name, 
                                                                 y=str(integral) if integral >= 0 else "",
                                                                 r=str(integralOnlyNonNegativeBin/integral) if options.printIntegralCapNegative else "",
                                                                 m=str(minBinContent) if (options.printMinNegative and minBinContent < 0) else "")
        elif options.output_format == "name":
            print "%s" % name


    if (obj.ClassName().startswith("TTree") or obj.InheritsFrom("TTree")):
        print "-"*30
        print "Branches of tree %s" % name
        print "-"*30
        lok = obj.GetListOfLeaves()
        for ip,p in enumerate(lok):
            print "{: <10} {l}".format(ip,l=p.GetName())
        print "-"*30
        

if options.printAxisLabelHisto != "" or options.printBranchTree != "":
    quit()

print "########################################"
print "########################################"
print " Summary"
print "----------------------------------------"
print "There were %d keys in the file" % tf.GetNkeys()
print "There were %d keys passing filters" % nRead
print "There were %d histograms with integral = 0" % nNullIntegral
print "There were %d histograms with one or more negative bins" % nNegativeBin

    


