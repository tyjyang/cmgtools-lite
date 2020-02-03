#!/usr/bin/env python

## WLIKE (2mu) skim:
##============================================ 
#       python  w-mass-13TeV/skims.py w-mass-13TeV/wlike_mu/skimming/mca-wlike.txt w-mass-13TeV/wlike_mu/skimming/skimCuts2mu.txt /afs/cern.ch/work/e/emanuele/TREES/TREES_wlike_mu_V1 /eos/cms/store/cmst3/group/wmass/w-mass-13TeV/ntuples/TREES_2LEP_wlike_mu_V1 -f w-mass-13TeV/wlike_mu/skimming/varsToKeep.txt --log <logdir> --mo --pretend -q condor
# to just check, remove the "-q condor" (but leave --pretend not to run interactively)
#
# Then, do the friends. Same command as above, but remove the "-f <varsToKeep> --mo" and add "--fo". And then:
# condor_submit <logdir>/condor_ftrees.condor
##============================================ 

## MUONS
##============================================
# DATA:
#       python  w-helicity-13TeV/skims.py w-helicity-13TeV/wmass_mu/skimming/mca-wmu-skim-bkg-data.txt w-helicity-13TeV/wmass_mu/skimming/skimCuts.txt /eos/cms/store/cmst3/user/mdunser/wMassTrees/2018-02-03-legacySingleMu/ /eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_2018_07_16_DATA_legacy/ -f w-helicity-13TeV/wmass_mu/skimming/varsToKeep.txt --mo -q 8nh --log skim_logs/
# MC:  
#       python  w-helicity-13TeV/skims.py w-helicity-13TeV/wmass_mu/skimming/mca-wmu-skim-bkg-data.txt w-helicity-13TeV/wmass_mu/skimming/skimCuts.txt /eos/cms/store/group/dpg_ecal/comm_ecal/localreco/TREES_1LEP_80X_V3/ /eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_2018-03-21_MC_1muskim/ -f w-helicity-13TeV/wmass_mu/skimming/varsToKeep.txt --mo

# SIGNAL:  
#       python  w-helicity-13TeV/skims.py w-helicity-13TeV/wmass_mu/skimming/mca-signal.txt w-helicity-13TeV/wmass_mu/skimming/signalCuts.txt /eos/cms/store/group/dpg_ecal/comm_ecal/localreco/TREES_1LEP_80X_V3_NoSkim5/ /eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_2018-03-21_SIGNAL_1muskim/ -f w-helicity-13TeV/wmass_mu/skimming/varsToKeep.txt

# SIGNAL FRIENDS:  
#       python  w-helicity-13TeV/skims.py w-helicity-13TeV/wmass_mu/skimming/mca-signal.txt w-helicity-13TeV/wmass_mu/skimming/signalCuts.txt /eos/cms/store/group/dpg_ecal/comm_ecal/localreco/TREES_1LEP_80X_V3_NoSkim5/ /eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_2018-03-21_SIGNAL_1muskim/ -f w-helicity-13TeV/wmass_mu/skimming/varsToKeep.txt

# add -q 8nh --log logs to run in batch 1 job/component (and --pretend to just check the command that will be run)


import os, subprocess, re

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] mc.txt cuts.txt treeDir outputDirSkims ")
    parser.add_option("-f", "--varfile",  dest="varfile", type="string", default=None, action="store",  help="File with the list of Branches to drop, as per TTree::SetBranchStatus")
    parser.add_option("--fo", "--friend-only",  dest="friendOnly", action="store_true", default=False,  help="Do not redo skim of the main trees, only of the friends")
    parser.add_option("--mo", "--main-only",  dest="mainOnly", action="store_true", default=False,  help="Do not make skim of the friend trees, only of the main")
    parser.add_option("--max-entries",     dest="maxEntries", default=1000000000, type="int", help="Max entries to process in each tree") 
    #parser.add_option("-n", "--job-name", dest="jobName",   type="string", default=None, help="Name assigned to jobs (if not given, use the default of skimTrees.py and skimFTrees.py"); # already in skimTrees,  cannot be duplicated
    from CMGTools.WMass.plotter.skimTrees import addSkimTreesOptions
    addSkimTreesOptions(parser)
    (options, args) = parser.parse_args() 

    mcargs = args[:2]
    treeDir = args[2]
    outputDirSkims = args[3]

    outputDirFSkims = outputDirSkims+"/friends"

    if not options.friendOnly:
        if not os.path.exists(outputDirSkims):
            os.makedirs(outputDirSkims)
            os.makedirs(outputDirFSkims)
        else:
            if options.pretend:
                print "The skim output dir ",outputDirSkims," exists, but you passed option --pretend"
                print "I guess you just want to print commands to test behaviour or resubmit few jobs"
                print "I will add suffix _v1 to the output folder to let you avoid making mistakes and overwriting things"            
                outputDirSkims = outputDirSkims.rstrip('/') + "_v1/"
                outputDirFSkims = outputDirSkims+"/friends"
                print "New output dir is " + outputDirSkims
            else:
                print "The skim output dir ",outputDirSkims," exists. Will remove it and substitute with new one. \nDo you agree?[y/N]\n"
                if raw_input()!='y':
                    print 'Aborting'
                    exit()
                os.system("rm -rf "+outputDirSkims)
                os.makedirs(outputDirSkims)
                os.makedirs(outputDirFSkims)
        os.system('cp {vf} {od}'.format(od=outputDirSkims,vf=options.varfile))
        os.system('cp {sf} {od}'.format(od=outputDirSkims,sf=args[1])) ## this should work??
    else: 
        print "Make only the friend trees in dir ",outputDirFSkims
        os.system('cp {vf} {od}'.format(od=outputDirSkims,vf=options.varfile)) ## this should work??

    OPTS = ' --obj tree -P '+treeDir+' --s2v -j 4 -F Friends "{P}/friends/tree_Friend_{cname}.root" '
    OPTS += ' --max-entries %d ' % options.maxEntries 
    BATCH_OPTS = ''
    if options.pretend: BATCH_OPTS += ' --pretend '
    if options.queue: BATCH_OPTS += ' -q %s ' % options.queue
    if options.logdir: BATCH_OPTS += ' --log %s ' % options.logdir
    if options.jobName != None: BATCH_OPTS += ' -n %s ' % options.jobName

    varsToKeep = []; DROPVARS = ''
    if options.varfile!=None:
        with open(options.varfile) as f:
            varsToKeep = f.read().splitlines()
        DROPVARS = " --dropall --keep "+" --keep ".join(varsToKeep)
        OPTS += DROPVARS
    
    cmdSkim = "python skimTrees.py "+" ".join(mcargs)+" " + outputDirSkims + OPTS + BATCH_OPTS
    cmdFSkimEv = " python skimFTrees.py "+outputDirSkims+" "+treeDir+"/friends/ "+outputDirFSkims+' -f tree_Friend -t "Friends" ' + DROPVARS + BATCH_OPTS
    if options.matchComponent:
        cmdSkim    += " --match-component '{reg}'".format(reg=options.matchComponent)
        cmdFSkimEv += " --match-component '{reg}'".format(reg=options.matchComponent)

    if not options.friendOnly:
        print "Now skimming the main trees, keeping the following vars:\n",varsToKeep
        print "This step may take time...\n"
        print cmdSkim
        os.system(cmdSkim)
    if not options.mainOnly:
        print "Now skimming the event variables friend trees:\n"
        os.system(cmdFSkimEv)
        # print "Now skimming the fake rate friend trees:\n"
        # os.system(cmdFSkimFr)
        # print "Now skimming the trigger friend trees:\n"
        # os.system(cmdFSkimTg)

    print "VERY DONE\n"


