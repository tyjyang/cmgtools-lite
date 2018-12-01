import re, sys, os, os.path

#folder = "/eos/cms/store/group/dpg_ecal/comm_ecal/localreco/TREES_1LEP_80X_V3/"
#folder = "/eos/cms/store/group/dpg_ecal/comm_ecal/localreco/TREES_1LEP_80X_V3_WENUSKIM_V5_TINY/"
#folder = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/plotter/TREE_4_XSEC_AFS/"
#folder = "/afs/cern.ch/work/e/emanuele/TREES/TREES_electrons_1l_V6_TINY/"
folder = "/eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_SIGSKIM_WMUNU_FULLSEL_NOMT/"
print "Folder: %s " % folder


#command = "ls %s | grep DYJetsToLL_M50_part" % folder
#print str(command)
#os.system(str(command))
#print "os.listdir(%s)" % folder
#quit()

#excludeDirMatch = [ 'NoSkim_WJetsToLNu_NLO' ]
excludeDirMatch = []

samples = os.listdir(folder)
print "Content:" 
for i,sample in enumerate(samples):
    print "%d) %s " % (i,sample)
print "--------------------------------------"
print ""

# regular expression pattern
# ^ means 'startswith'
# sampleRoots = ['^DYJetsToLL_M50.*', 
#                '^WW.*',
#                '^WZ.*',
#                '^ZZ.*',
#                '^TBar_tWch_ext.*',
#                '^T_tWch_ext.*',
#                '^TTJets_SingleLeptonFromT_.*',
#                '^TTJets_SingleLeptonFromTbar_.*',
#                '^T_tch_powheg.*',
#                '^TToLeptons_sch_amcatnl.*',
#                '^TBar_tch_powheg_.*',
#                '^WJetsToLNu_NLO.*',
#                '^WJetsToLNu_LO.*',
#                '^NoSkim_WJetsToLNu_NLO.*',
# #               '.*_Mu5.*',
# #               '.*_Mu15.*',
#                '.*_bcTo.*',
#                '.*_EMEnriched.*'
#                ]
sampleRoots = ['^WJetsToLNu_.*',               
               ]

sumgenwgt = {}

for sampleRoot in sampleRoots:

    genwgt = {}
    for sample in samples:
        if any(x in str(sample) for x in excludeDirMatch): continue

        #if str(sampleRoot) in str(sample):
        if re.match(sampleRoot,str(sample)):
            #print str(sample)
            filename= folder + sample + "/skimAnalyzerCount/SkimReport.txt"
            with open(str(filename), 'r') as thisfile:
                for line in thisfile:
                    if "Sum Weights" in str(line):
                        parts = line.split()
                        #print parts
                        #print str(parts[2].strip())
                        genwgt[str(sample)] = float(parts[2].strip())

    isFirst = True
    sortedKeys = sorted(genwgt.keys())
    for key in sortedKeys:
        print "%s : %s" % (key, genwgt[key])
        if isFirst:
            sumgenwgt[str(sampleRoot)] = float(genwgt[key])
        else:
            sumgenwgt[str(sampleRoot)] += float(genwgt[key])
        isFirst = False

    print "####"
    print "N(keys) = %d" % len(sortedKeys)
    print "sumgenwgt[%s] : %s" % (str(sampleRoot), sumgenwgt[str(sampleRoot)])
    print "####"
    print ""

# for key in sumgenwgt:
#     print "sumgenwgt[%s] : %s" % (key, sumgenwgt[key])


map_sampleName_sumGenWgt = {}
samplesNotMatchedToKey = []

for sample in samples:
    if any(x in str(sample) for x in excludeDirMatch): continue
    sampleMatchedToKey = False
    for key in sumgenwgt:
        if re.match(key,str(sample)):
            map_sampleName_sumGenWgt[str(sample)] = sumgenwgt[key]
            print 'map_sampleName_sumGenWgt["%s"] = %s;' % (sample, sumgenwgt[key])
            sampleMatchedToKey = True
            break
    if not sampleMatchedToKey:
        samplesNotMatchedToKey.append(str(sample))

print ""
print "List of samples not matched to any key"
for i,sample in enumerate(samplesNotMatchedToKey):
    print "%s : %s" % (i,sample)


print ""
print "List of matches used to exclude folders"
for i,dirmatch in enumerate(excludeDirMatch):
    print "%s : %s" % (i,dirmatch)

# i = 0
# for folder in folders:
#     print "%s --> %s " % (i, folder)
