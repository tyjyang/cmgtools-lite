#!/usr/bin/env python

import os, subprocess, re
import ROOT

# this script is meant to run only the stat variations, the nominal can be run locally, it is fast

def checkClosureFile(indir,closurefile="plot_dm_diff.root"):

    fname = indir+"/"+closurefile
    if not os.path.isfile(fname): 
        return False
    tf = ROOT.TFile.Open(fname)
    if not tf:
        print "ERROR! Unable to open file %s" % (fname)
        return False
    if tf.TestBit(ROOT.TFile.kRecovered):
        print "WARNING! Attempt was made to recover file %s" % (fname)
        return False
    h = tf.Get("plot_dm_diff" if "_diff" in closurefile else "plot_dm")
    if not h:
        print "ERROR! %s is missing in %s" % (closurefile,fname)
        return False
    return True

def checkHistogramFile(indir):

    fname = indir+"/histo3D_mass_pt_eta.root"
    tf = ROOT.TFile.Open(fname)
    if not tf:
        #print "ERROR! Unable to open file %s" % (fname)
        return False
    if tf.TestBit(ROOT.TFile.kRecovered):
        #print "WARNING! Attempt was made to recover file %s" % (fname)
        return False
    hdata = tf.Get("data_histo_mass_pt_eta")
    hmc = tf.Get("mc_histo_mass_pt_eta")
    if not hdata or not hmc:
        #print "ERROR! At least one TH3D histogram is missing in %s" % (fname)
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
    parser.add_option("--closure-filename", dest="closureFilename",   type="string", default="plot_dm.root", help="Name of closure file when checking jobs");
    parser.add_option("--run-failed", dest="runFailed", default=False, action="store_true",  help="Just run failed jobs (assumes you already ran once but some histograms are missing)") 
    parser.add_option("--only-check-failed", dest="onlyCheckFailed", default=False, action="store_true",  help="Just check failed jobs") 
    parser.add_option("--only-check-histogram", dest="onlyCheckHistogram", default=False, action="store_true",  help="Just check failed histograms, not also closure file (this options sets --only-check-failed to True, but only check for histograms)") 
    parser.add_option("--only-delete-fits", dest="onlyDeleteFits", default=False, action="store_true",  help="Utility option to delete fits when no longer needed, to save space") 
    parser.add_option("--draw-fits",   dest="drawFits", default=False, action='store_true' , help="Save plots for fits (takes quite a lot of space for all the stat variations)");
    parser.add_option("--ptbins",   dest="ptbins", type="string", default="23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65", help="Comma separated list of pt bin edges")
    parser.add_option("--etabins",   dest="etabins", type="string", default="-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4", help="Comma separated list of eta bin edges")
    parser.add_option("--rebin-pt",   dest="rebinPt", type="int", default="1", help="Rebin pt by this integer value when fitting")
    parser.add_option("--rebin-eta",   dest="rebinEta", type="int", default="1", help="Rebin eta by this integer value when fitting")
    parser.add_option("--ntuples-path",   dest="ntuplesPath", type="string", default="/eos/cms/store/cmst3/group/wmass/mciprian/heppyNtuples/skim_Zmumu/", help="Path to ntuples")
    parser.add_option("--make-only-histo",    dest="makeOnlyHisto", default=False, action='store_true' , help="only produce histogram and exit. Useful to make a very finely binned histogram to allows making fits later with any desired rebinning")
    parser.add_option("--folder-load-histo",    dest="folderLoadHisto", type="string", default="", help="When loading existing histogram, by default it is taken from same folder as output, but can choose a different one")
    (options, args) = parser.parse_args() 

    ntuplesPath = options.ntuplesPath

    if options.onlyCheckHistogram:
        options.onlyCheckFailed = True

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
            outdir = options.outdir.rstrip("/") + "/" + charge + "/RoccorStatVar/"
            if not os.path.exists(outdir):
                os.system("mkdir -p "+outdir)

        if options.folderLoadHisto:
            folderLoadHisto = options.folderLoadHisto.rstrip("/") + "/" + charge + "/RoccorStatVar/"

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
            pdirFull = pdir + charge + "/"

            if options.folderLoadHisto:
                # when loading, the full path to file must be given, hence the ending charge folder
                thisFolderLoadHisto = "{o}/charge{c}_stat{d}/{ch}/".format(o=folderLoadHisto,c="Plus" if charge=="plus" else "Minus",d=istat,ch=charge)

            if options.onlyDeleteFits:
                rmcmd = "rm {pd}plot_pt_*_eta_*.p*".format(pd=pdirFull)
                if options.dryRun: 
                    print rmcmd
                else:
                    os.system(rmcmd)
                continue

            if options.onlyCheckFailed:
                if options.onlyCheckHistogram:
                    ret = checkHistogramFile(pdirFull)                    
                else:
                    ret = checkClosureFile(pdirFull,closurefile=options.closureFilename)
                if ret == False:
                    failedDir.append(pdirFull)
                continue

            runOnlyFit = False
            if options.runFailed:
                ret1 = checkHistogramFile(pdirFull)
                if ret1:
                    ret = checkClosureFile(pdirFull,closurefile=options.closureFilename)
                    if ret:
                        # all is ok, skip this job
                        continue
                    else:
                        # to be redone, but histogram is present, can just run fit
                        runOnlyFit = True
                else:
                    # to be redone completely
                    pass

            # default is to use template fit with SCALE
            cmd = "python w-mass-13TeV/zFitterNew.py "
            cmd += " --data {p} --refmc {p}".format(p=ntuplesPath)
            cmd += " -c Zmm-full --pdir {pd}".format(pd=pdir)
            cmd += " --ptvar {pt} ".format(pt=ptvar)
            cmd += " -x  mass_2({pt}[0],LepGood_eta[0],LepGood_phi[0],0.1057,{pt}[1],LepGood_eta[1],LepGood_phi[1],0.1057) 800,70,110 ".format(pt=ptvar)
            cmd += " --select-charge {c} --useAllEvents".format(c=charge)
            cmd += " --setRangeClosure 0.001  -s MC-SCALE -b None -p dm "
            cmd += " --fit-strategy 2 --roofitPrintLevel -1 "
            cmd += " --plot-extension png --drawTH2option 'COLZ0TEXT45E' "
            cmd += " --ptbins {pt} --etabins {eta}".format(pt=options.ptbins,eta=options.etabins)
            # hardcoded rebinning
            cmd += " --rebin 5 --rebinTemplateModel 2 --templateSmoothOrder 3 "
            if options.rebinPt > 1:
                cmd += " --rebin-pt {d}".format(d=options.rebinPt) 
            if options.rebinEta > 1:
                cmd += " --rebin-eta {d}".format(d=options.rebinEta) 
            if options.onlyFit or runOnlyFit:
                cmd += " --loadHistoFromFile {pd}histo3D_mass_pt_eta.root ".format(pd=pdir if not options.folderLoadHisto else thisFolderLoadHisto)
            if options.makeOnlyHisto:
                cmd += " --makeOnlyHisto "
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
            if options.onlyCheckHistogram:
                print "Following dirs miss the histogram file (didn't check for the closure)"
            else:
                print "Following dirs miss the final file with closure"
            print "-"*30
            for i,d in enumerate(failedDir):
                print "%d) %s" % (i,d)
            print "-"*30
            print "Need to resubmit %d jobs :(" % len(failedDir)
        else:
            print "All files are present :)\n\n"
