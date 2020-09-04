#!/usr/bin/env python
import os, sys
import re

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] inputDir outputDir")
    parser.add_option("-p", "--pretend", dest="pretend",   action="store_true", default=False, help="Don't run anything");
    parser.add_option("-v", "--verbose", dest="verbose",   action="store_true", default=False, help="Use option -v for cp");
    parser.add_option("-m", "--match",   dest="match",     type="string", default=None, help="Select only folders matching this regular expression")
    parser.add_option("-r", "--replace", dest="replace",   type="string", nargs=2, default=None, help="Replace first argument with second argument when copying folder or file")
    (options, args) = parser.parse_args()

if len(args)<2:
    parser.print_help()
    sys.exit(1)

inputdir = args[0]
outputDir = args[1]
if inputdir == outputDir:
    print "Error: input folder is the same as the output folder."
    quit()

print "-"*30
print "Copying folders from inputdir %s to outputDir %s" % (inputdir, outputDir)
if options.match != None:
    print "Selecting folders matching '{reg}'".format(reg=options.match)
if options.replace != None:
    print "Changing name replacing '%s' with '%s'" % (options.replace[0], options.replace[1])
print "-"*30

if not os.path.exists(outputDir):
    print "Warning: folder {out} does not exist".format(out=outputDir)
    print "It will be created now"
    os.makedirs(outputDir)

folders = os.listdir(inputdir)
for fld in folders:
    if options.match == None or re.match(options.match,fld):
        cmd = "cp -r{v} {ind}/{f} {out}/{newf}".format(v="v" if options.verbose else "", 
                                                       ind=inputdir,f=fld, out=outputDir, 
                                                       newf=fld if options.replace == None else fld.replace(options.replace[0], options.replace[1]))
        if options.pretend:
            print cmd
        else:
            os.system(cmd)
            

print "-"*30
print "DONE"
print "-"*30
