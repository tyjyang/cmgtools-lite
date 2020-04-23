#!/usr/bin/env python

import os, subprocess, re
import ROOT

# this script is meant to run only the stat variations, the nominal can be run locally, it is fast

def checkClosureFile(indir):

    fname = indir+"/plot_dm_diff.root"
    if not os.path.isfile(fname): 
        return False
    tf = ROOT.TFile.Open(fname)
    if not tf:
        print "ERROR! Unable to open file %s" % (fname)
        return False
    if tf.TestBit(ROOT.TFile.kRecovered):
        print "WARNING! Attempt was made to recover file %s" % (fname)
        return False
    h = tf.Get("plot_dm_diff")
    if not h:
        print "ERROR! plot_dm_diff is missing in %s" % (fname)
        return False
    return True

def checkHistogramFile(indir):

    fname = indir+"/histo3D_mass_pt_eta.root"
    tf = ROOT.TFile.Open(fname)
    if not tf:
        print "ERROR! Unable to open file %s" % (fname)
        return False
    if tf.TestBit(ROOT.TFile.kRecovered):
        print "WARNING! Attempt was made to recover file %s" % (fname)
        return False
    hdata = tf.Get("data_histo_mass_pt_eta")
    hmc = tf.Get("mc_histo_mass_pt_eta")
    if not hdata or not hmc:
        print "ERROR! At least one TH3D histogram is missing in %s" % (fname)
        return False
    return True


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options]")
    parser.add_option("-d", "--dry-run", dest="dryRun", default=False, action="store_true",  help="Just print commands") 
    parser.add_option("-f", "--only-fit", dest="onlyFit", default=False, action="store_true",  help="Just fit, do not run on ntuples (only works if histograms exist already)") 
    parser.add_option("--log", "--log-dir", dest="logdir", type="string", default=None, help="Directory of stdout and stderr");
    parser.add_option("-o", "--outdir", dest="outdir", type="string", default=None, help="Output folder with the nominal closure: this contains the two charge subfolders, and a folder named 'RoccorStatVar' will be created inside each of them");
    parser.add_option("-n", "--job-name", dest="jobName",   type="string", default="zFitterNew", help="Name assigned to jobs");
    parser.add_option("--run-failed", dest="runFailed", default=False, action="store_true",  help="Just run failed jobs (assumes you already ran once but some histograms are missing)") 
    parser.add_option("--only-check-failed", dest="onlyCheckFailed", default=False, action="store_true",  help="Just check failed jobs") 
    parser.add_option("--only-delete-fits", dest="onlyDeleteFits", default=False, action="store_true",  help="Utility option to delete fits when no longer needed, to save space") 
    parser.add_option("--draw-fits",   dest="drawFits", default=False, action='store_true' , help="Save plots for fits (takes quite a lot of space for all the stat variations)");
    parser.add_option("--ptbins",   dest="ptbins", type="string", default="23,30,35,40,45,50,55,60", help="Comma separated list of pt bin edges")
    parser.add_option("--etabins",   dest="etabins", type="string", default="-2.4,-1.6,-0.8,0,0.8,1.6,2.4", help="Comma separated list of eta bin edges")
    parser.add_option("--ntuples-path",   dest="ntuplesPath", type="string", default="/eos/cms/store/cmst3/group/wmass/mciprian/heppyNtuples/skim_Zmumu/", help="Path to ntuples")
    (options, args) = parser.parse_args() 

    ntuplesPath = options.ntuplesPath


    if options.onlyCheckFailed and options.runFailed:
        print "Warning: --run-failed is not compatible with --only-check-failed"
        print "Perhaps you only wanted to check missing files but not run any job"        
        quit()

    if options.onlyCheckFailed:
        print "Will just check bad folders ..."
    if options.runFailed:
        print "Will only run failed jobs (it assumes you already tried once)"
    if options.onlyDeleteFits:        
        print "Will only delete fits to save space"

    ## constructing the command and arguments to run in condor submit file
    runner = "%s/src/CMGTools/WMass/python/postprocessing/lxbatch_runner.sh" % os.environ['CMSSW_BASE']
 
   ## first make the log directory to put all the condor files
    writelog = ""
    logdir   = ""
    if not options.logdir: 
        print 'ERROR: must give a log directory to make the condor files!\nexiting...'
        exit()
    else:
        logdir = options.logdir.rstrip("/")
        if not os.path.exists(logdir):
            os.system("mkdir -p "+logdir)

    failedDir = []

    for charge in ["plus", "minus"]:

        if not options.outdir: 
            print 'ERROR: must give a directory where to store output!\nexiting...'
            exit()
        else:
            outdir = options.outdir.rstrip("/") + charge + "/RoccorStatVar/"
            if not os.path.exists(outdir):
                os.system("mkdir -p "+outdir)


        condor_fn = options.logdir+'/condor_zFitter_{c}.condor'.format(c=charge)
        condor_f = open(condor_fn,'w')
        condor_f.write('''Universe = vanilla
Executable = {runner}
use_x509userproxy = true
Log        = {ld}/{c}_$(ProcId).log
Output     = {ld}/{c}_$(ProcId).out
Error      = {ld}/{c}_$(ProcId).error
getenv      = True
request_memory = 2000
+MaxRuntime = 7200
+JobBatchName = "{name}"\n
'''.format(runner=runner,c=charge,ld=options.logdir,name=options.jobName))
        if os.environ['USER'] in ['mdunser', 'psilva']:
            condor_f.write('+AccountingGroup = "group_u_CMST3.all"\n')
        elif os.environ['USER'] in ['mciprian']:
            condor_f.write('+AccountingGroup = "group_u_CMS.CAF.ALCA"\n')
        condor_f.write('\n')

        for istat in range(100):

            ptvar = "LepGood_rocPt" + ("_stat%d" % istat)
            pdir = "{o}/charge{c}_stat{d}/".format(o=outdir,c="Plus" if charge=="plus" else "Minus",d=istat)

            if options.onlyDeleteFits:
                rmcmd = "rm {pd}plot_pt_*_eta_*.p*".format(pd=pdir)
                if options.dryRun: 
                    print rmcmd
                else:
                    os.system(rmcmd)
                continue

            if options.onlyCheckFailed:
                ret = checkClosureFile(pdir)
                if ret == False:
                    failedDir.append(pdir)
                continue

            runOnlyFit = False
            if options.runFailed:
                ret1 = checkHistogramFile(pdir)
                if ret1:
                    ret = checkClosureFile(pdir)
                    if ret:
                        # all is ok, skip this job
                        continue
                    else:
                        # to be redone, but histogram is present, can just run fit
                        runOnlyFit = True
                else:
                    # to be redone completely
                    pass

            cmd = "python w-mass-13TeV/zFitterNew.py "
            cmd += " --data {p} --refmc {p}".format(p=ntuplesPath)
            cmd += " -c Zmm-full --pdir {pd}".format(pd=pdir)
            cmd += " --ptvar {pt} ".format(pt=ptvar)
            cmd += " -x  mass_2({pt}[0],LepGood_eta[0],LepGood_phi[0],0.1057,{pt}[1],LepGood_eta[1],LepGood_phi[1],0.1057) 80,70,110 ".format(pt=ptvar)
            cmd += " --select-charge {c} --useAllEvents".format(c=charge)
            cmd += " --setRangeClosure 0.1  -s Z-DCB --roofitPrintLevel -1 --plot-extension png "
            cmd += " --ptbins {pt} --etabins {eta}".format(pt=options.ptbins,eta=options.etabins)
            if options.onlyFit or runOnlyFit:
                cmd += " --loadHistoFromFile {pd}histo3D_mass_pt_eta.root ".format(pd=pdir)
            if not options.drawFits:
                cmd += " --no-draw-fits "

            cmdargs = [x.strip() for x in cmd.split()]
            strargs=''
            for a in cmdargs: # join do not preserve " or '
                ## escaping not needed for condor if '{P}' in a or '*' in a: 
                ## escaping not needed for condor     a = '\''+a+'\''
                strargs += ' '+a+' '


            condor_f.write('\narguments = {d} {cmssw} {cmdargs}\n'.format(d=os.getcwd(), cmssw = os.environ['CMSSW_BASE'], cmdargs=strargs))
            condor_f.write('queue 1\n\n')

        if options.onlyDeleteFits:
            continue

        condor_f.close()
        cmd = 'condor_submit {sf}'.format(sf=condor_fn)

        if not options.onlyCheckFailed:
            if options.dryRun: 
                print cmd
            else:
                os.system(cmd)

    if options.onlyDeleteFits:
        quit()

    if options.onlyCheckFailed:
        if len(failedDir):
            print "Following dirs miss the final file with closure"
            print "-"*30
            for i,d in enumerate(failedDir):
                print "%d) %s" % (i,d)
            print "-"*30
            print "Need to resubmit %d jobs :(" % len(failedDir)
        else:
            print "All files are present :)\n\n"
