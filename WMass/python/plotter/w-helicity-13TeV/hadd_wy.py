# usage:  python w-helicity-13TeV/hadd_wy.py <input_folder> [options]
#
# first argument is list of folders with distributions for no or nominal selection
# See examples below

# script developed to work passing a folder containing other folders like the following:
#
# [mciprian@pccmsrm29 plotter]$ ls plots/gen/
# wgen_fullsel_minus         wgen_fullsel_pfmt30_plus   wgen_fullsel_pfmt40_plus   wgen_fullsel_pfmt50_plus  wgen_nosel_minus
# wgen_fullsel_pfmt30_minus  wgen_fullsel_pfmt40_minus  wgen_fullsel_pfmt50_minus  wgen_fullsel_plus         wgen_nosel_plus
# [mciprian@pccmsrm29 plotter]$ python w-helicity-13TeV/hadd_wy.py plots/gen/ -c pfmt
#
# it will create root files containing histograms with no selection and with any of the other selections (no pfmt or pfmtXX in this case).
# For the example above, it will create 4 files, one for nominal and one for each of the 3 pfmt selections
#
# It works also in you have only nosel and nominal selection (4 folders), which is the default.
#
# In general we use NLO samples, but we could also use LO (for example, for systematics). 
# In this case the convention is that the last folder must end with '_LO'
# In that case, the histograms will also be named with '_LO' at the end of their original names, and will be put in the same output file as the others.
# This requires passing option -l, otherwise the script complains becasue it sees more than 4 folders 
# (if you have exactly 4 folders and all are named with _LO the script still works).
# The rest works in the same way (if you have additional selections and so on)
#
# Note: the script should still work if all the fullsel folders have some cut tag in their name, like in the following:
#
# [mciprian@pccmsrm29 plotter]$ ls plots/gen_eff_tightCharge/     
# wgen_fullsel_pfmt40_minus     wgen_fullsel_pfmt40_plus     wgen_nosel_minus     wgen_nosel_plus
# wgen_fullsel_pfmt40_minus_LO  wgen_fullsel_pfmt40_plus_LO  wgen_nosel_minus_LO  wgen_nosel_plus_LO
#
# In this case, you may not specify the cut with -c <cut_name_ID> option, and the files will be treated as usual, but you need -l

import sys, os
import ROOT, datetime, array
import re

from optparse import OptionParser
parser = OptionParser(usage='python %prog <input_folder> [options] ')
parser.add_option('-c','--cut-name', dest='cutName', default='', type='string', help='name of cut. It is a tag that should be present in some folders')
parser.add_option('-l','--LO', dest="hasLO", action="store_true", default=False, help="Specify if there are folders for LO samples (they must end with '_LO')")
(options, args) = parser.parse_args()

if len(args) == 0:
    parser.print_usage()
    quit()

inputdir = args[0]

# varcut = "pfmt"
# if varcut != "":
#     noAdditionalCut = False
# else:
#     noAdditionalCut = True

varcut = ""
noAdditionalCut = True
if options.cutName != '':
    varcut = options.cutName
    noAdditionalCut = False

#files = [ f for f in os.listdir(inputdir) if f.endswith('.root') ]
#files = list( [os.path.join(inputdir, f) for f in files] )

files = list()
for root, dirs, tmpfiles in os.walk(inputdir):
    for f in tmpfiles:
        if f.endswith(".root"):
            thisfile = os.path.join(root, f)
            print("Getting file --> %s " % str(thisfile))
            files.append(str(thisfile))

if options.hasLO:
    expectedFolders = 8
else:
    expectedFolders = 4

if len(files) > expectedFolders and noAdditionalCut:
    print "==================================="
    print "WARNING: you have not specified any folder name tag, but I see more than %d folders (%d)" % (int(expectedFolders),len(files))
    if int(expectedFolders) == 4 and len(files) == 8:
        print "Did you want to include LO samples but forgot option --LO ?"
    print "I cannot manage this situation correctly, please check. Exit"
    print "==================================="
    quit()

# for f in files:
#     print f
# quit()

helicities = ["right", "left", "long"]

tmpplots = []
varcut_thr_list = set([])
for f in files:
    print "Opening file: ",f
    tf = ROOT.TFile.Open(f)
    for k in tf.GetListOfKeys() :
        name=k.GetName()
        obj=k.ReadObj()
        #if '_wy_' in name and 'background' not in name and obj.InheritsFrom("TH1"):
        if '_wy_' in name and obj.InheritsFrom("TH1") and any(h in name for h in helicities):
            if 'fullsel' in f:
                tokens = name.split('_')
                if varcut != "" and varcut in f:
                    regex = re.compile(varcut+'([0-9]*)')
                    varcut_thr = regex.findall(f)
                    if len(varcut_thr):
                        #print int(varcut_thr[0])
                        varcut_thr_list.add(int(varcut_thr[0]))
                        newname = '_'.join( tokens[:2]+['reco_%s%d' % (varcut, int(varcut_thr[0]))]+tokens[2:] )                        
                else: 
                    newname = '_'.join( tokens[:2]+['reco']+tokens[2:] )
            else:
                newname = name
            lastFolder = os.path.basename(os.path.dirname(f))
            if lastFolder.endswith('_LO'):
                newname = newname + '_LO'
            newh = obj.Clone(newname)
            newh.SetDirectory(None)
            tmpplots.append(newh)
    #tf.Close()

outputfile = 'mc_reco_eff.root'
mergedFile = ROOT.TFile.Open(outputfile,'recreate')
mergedFile.cd()

print "#####################"
print "File ", outputfile
for p in tmpplots:
    if noAdditionalCut or varcut not in p.GetName():
        print "Writing histo: ",p.GetName()
        p.Write()
mergedFile.Close()


if not noAdditionalCut:

    print ""    
    print "thresholds for " + varcut + " cut: ",varcut_thr_list
    for thr in varcut_thr_list:
        outputfile = 'mc_reco_%s%d_eff.root' % (varcut, thr) 
        mergedFile = ROOT.TFile.Open(outputfile,'recreate')
        mergedFile.cd()
        print "#####################"
        print "File ", outputfile
        for p in tmpplots:
            if 'reco' not in p.GetName() or (str(thr) in p.GetName() and varcut in p.GetName()):
                tmpNameNoCut = p.GetName().replace("%s%d_" % (varcut, thr),"")
                print "%s%d: Writing histo: %s (original name --> %s)" % (varcut, thr, tmpNameNoCut, p.GetName())
                p.Write(tmpNameNoCut)
        mergedFile.Close()
    