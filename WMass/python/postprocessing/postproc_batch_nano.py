#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from glob import glob
import re, pickle, math
from CMGTools.WMass.postprocessing.framework.postprocessor import PostProcessor

## USAGE:

## to run locally:
# python postproc_batch.py -N 250000  <directoryWithTrees> <targetDirectoryForFriends>  --friend  
# optional: -d <datasetToRun> -c <chunkToRun>

## if you want to submit to condor, add:
## --log friends_log/ --submit  --runtime 240 (optional, default 480 minutes) --memory 2000 (or more if needed)

## if you insist on running on LSF, instead add:
## --log friends_log/ --submit --env lxbatch --queue 1nd (optional, default 8nh)

DEFAULT_MODULES = [("CMGTools.WMass.postprocessing.examples.puWeightProducer", "puWeight"),
                   ("CMGTools.WMass.postprocessing.examples.triggerMatchProducer", "muTrigMatch")
                   #("CMGTools.WMass.postprocessing.examples.lepVarProducer","muCalibratedWithStatVar"),
                   #("CMGTools.WMass.postprocessing.examples.jetReCleaner","jetReCleaner,jetAllReCleaner"),
                   #("CMGTools.WMass.postprocessing.examples.genFriendProducer","genQEDJets"),
                   ]

def writeCondorCfg(logdir, name, flavour=None, maxRunTime=None,maxMemory=4000):

    #if not os.path.isfile(logdir+'/dummy_exec.sh'):
    dummy_exec = open(logdir+'/dummy_exec.sh','w') 
    dummy_exec.write('#!/bin/bash\n')
    dummy_exec.write('cd {here}\n'.format(here=os.environ['PWD']))
    dummy_exec.write('eval `scramv1 runtime -sh` \n')
    dummy_exec.write('python $*\n')
    dummy_exec.close()

    maxruntime="-t" # time in minutes
    if maxRunTime:
        if flavour:
            print "Can't set both flavour and maxruntime"
            sys.exit(1)
        maxruntime = str(60 * int(maxRunTime))
    job_desc = """Universe = vanilla
Executable = {ld}/dummy_exec.sh
use_x509userproxy = true
Log        = {ld}/{name}_$(ProcId).log
Output     = {ld}/{name}_$(ProcId).out
Error      = {ld}/{name}_$(ProcId).error
getenv      = True
environment = "LS_SUBCWD={here}"
request_memory = {mem}
""".format(ld=logdir,
           name=name,
           mem=maxMemory,
           here=os.environ['PWD'])
    if flavour:
        job_desc += '+JobFlavour = "%s"\n' % flavour
    if maxruntime!="":
        job_desc += '+MaxRuntime = %s\n' % maxruntime
    if os.environ['USER'] in ['mdunser', 'psilva']:
        job_desc += '+AccountingGroup = "group_u_CMST3.all"\n'
    if os.environ['USER'] in ['mciprian']:
        job_desc += '+AccountingGroup = "group_u_CMS.CAF.ALCA"\n'
    ##job_desc += 'queue 1\n'

    jobdesc_file = logdir+'/'+name+'.condor'
    with open(jobdesc_file,'w') as outputfile:
        outputfile.write(job_desc)
        outputfile.write('\n\n')
        outputfile.close()
    return jobdesc_file

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] inputDir outputDir")
    parser.add_option("-J", "--json",  dest="json", type="string", default=None, help="Select events using this JSON file")
    parser.add_option("-C", "--cut",  dest="cut", type="string", default=None, help="Cut string")
    parser.add_option("-b", "--branch-selection",  dest="branchsel", type="string", default=None, help="Branch selection")
    parser.add_option("--friend",  dest="friend", action="store_true", default=True, help="Produce friend trees in output (current default is to produce full trees)")
    parser.add_option("--full",  dest="friend", action="store_false",  default=False, help="Produce full trees in output (this is the current default)")
    parser.add_option("--noout",  dest="noOut", action="store_true",  default=False, help="Do not produce output, just run modules")
    parser.add_option("--justcount",   dest="justcount", default=False, action="store_true",  help="Just report the number of selected events") 
    parser.add_option("-I", "--import", dest="imports",  type="string", default=[], action="append", nargs=2, help="Import modules (python package, comma-separated list of ");
    parser.add_option("-z", "--compression",  dest="compression", type="string", default=("LZMA:9"), help="Compression: none, or (algo):(level) ")
    parser.add_option("-d", "--dataset", dest="datasets",  type="string", default=[], action="append", help="Process only this dataset (or dataset if specified multiple times)");
    parser.add_option(      "--dataset-rgx", dest="datasetsRgx",  type="string", default=None, help="Process only datasets that match this regexp (or comma-separated list of regexp). Usually it can be used to match a single root file when making tests");
    parser.add_option("-c", "--chunk",   dest="chunks",    type="int",    default=[], action="append", help="Process only these chunks (mainly for local tests)");
    parser.add_option("-N", "--events",  dest="chunkSize", type="int",    default=1000000, help="Default chunk size when splitting trees");
    parser.add_option("-p", "--pretend", dest="pretend",   action="store_true", default=False, help="Don't run anything");
    parser.add_option("-j", "--jobs",    dest="jobs",      type="int",    default=1, help="Use N threads");
    parser.add_option("-q", "--queue",   dest="queue",     type="string", default='8nh', help="Run jobs on lxbatch instead of locally");
    parser.add_option("-r", "--runtime",    type=int, default=480, help="Condor runtime. In minutes.");
    parser.add_option(      "--memory",     type=int, default=2000, help="Condor memory. In MB.");
    parser.add_option("-s", "--submit",    action='store_true', default=False, help="Submit the jobs to lsf/condor.");
    parser.add_option("-t", "--tree",    dest="tree",      default='Events', help="Pattern for tree name");
    parser.add_option("--log", "--log-dir", dest="logdir", type="string", default=None, help="Directory of stdout and stderr");
    parser.add_option("--env",   dest="env", type="string", default="condor", help="Give the environment on which you want to use the batch system (lsf,condor). Default: condor");
    parser.add_option("--run",   dest="runner",  type="string", default="lxbatch_runner.sh", help="Give the runner script (default: lxbatch_runner.sh)");
    parser.add_option("--mconly", dest="mconly",  action="store_true", default=False, help="Run only on MC samples");
    parser.add_option("--signals", dest="signals", default="DYJetsToMuMu,WplusJetsToMuNu,WminusJetsToMuNu", help="declare signals (CSV list) [%default]",type='string');
    parser.add_option("-m", "--modules", dest="modules",  type="string", default=[], action="append", help="Run only these modules among the imported ones");
    parser.add_option(      "--moduleList", dest="moduleList",  type="string", default='DEFAULT_MODULES', help="use this list as a starting point for the modules to run [%default]")

    (options, args) = parser.parse_args()

    if options.friend:
        if options.cut or options.json: raise RuntimeError("Can't apply JSON or cut selection when producing friends")

    if len(args) != 2:
        #print 'this is args:', args
        parser.print_help()
        sys.exit(1)
    #if len(options.chunks) != 0 and len(options.datasets) != 1:
    #    print "must specify a single dataset with -d if using -c to select chunks"
    #    sys.exit(1)

    if len(options.datasets) == 0 and options.friend:
        print "Must pass one or more dataset names when making friends, these will be used to "
        print "create the subfolders containing the friend trees."
        print "This will have to match a part of the path to trees passed as args[0]."
        print "For instance, one has to pass Run2016G for data (not SingleMuon_Run2016G), "
        print "or WplusJetsToMuNu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos for MC."
        print "The match could be smaller (e.g. WplusJetsToMuNu), but better to have the full field "
        print "that appears between two '/' in the path"
        sys.exit(1)

    treedir = args[0]; outdir=args[1]; args = args[2:]

    jobs = []

    # treedir path must be given in a form acceptable for glob
    # e.g. data for Run2016G is in /eos/cms/store/data/Run2016G/SingleMuon/NANOAOD/Nano02Dec2019-v1/
    # data for Run2016H is in /eos/cms/store/data/Run2016H/SingleMuon/NANOAOD/Nano02Dec2019-v1/
    # so the path change in the middle, and then each folder has a set of subfolders (e.g. /270000/) 
    # containing the trees, which are named like 60D17F68-FE0D-E849-BF40-49E01AC21621.root
    # so one has to pass a path containing some regular expression to use all data, unless one run 
    # this script independently for each path
    # for MC we currently have trees in another location, with a completely different logic in the path name
    #
    # for data one can pass 
    # /eos/cms/store/data/Run2016*/SingleMuon/NANOAOD/Nano02Dec2019-v1/*/*
    # so glob(treedir will return the root files with full path)
    for fnameWithPath in glob(treedir):
        fname = fnameWithPath
        if not fname.endswith(".root"): continue
        
        if os.path.exists(fname):
            shortfname = os.path.basename(fnameWithPath)
            if options.datasets != []:
                # dataset can be like Run2016G, Run2016H for data (i.e. no longer SingleMuon_Run2016G_part3 as it was for CMGTools ntuples)
                # and WplusJetsToMuNu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos for signal MC
                # this appears explicitly in the full path to relevant root files (i.e. in fnameWithPath)
                # one can pass more datasets, but it is enough that only one matches.
                # This is probably no longer needed, before we had one tree for each dataset (where 'dataset' in this context was a DD_partXX of the actual dataset DD), 
                # while now we have a dataset (or few for data, named as above) and all trees inside
                foundDataset = False
                for x in options.datasets:
                    if x in fnameWithPath:
                        dataset = x
                        foundDataset = True
                if not foundDataset: continue

            if options.datasetsRgx != None:
                regexps = options.datasetsRgx.split(',')
                if not any(re.match(rgx,fnameWithPath) for rgx in regexps): continue

            primaryDatasetsData = "DoubleMu DoubleEG MuEG MuonEG SingleMuon SingleElectron".split()
            data = False
            for x in primaryDatasetsData:
                if x in fnameWithPath:
                    data = True
                    #dataset = x + "_" + dataset

            if data and options.mconly: continue
            ## construct proper xrootd file path
            prepath = ''
            if not 'root:/' in fname:
                if   '/eos/user/'      in fname: prepath = 'root://eosuser.cern.ch//'
                elif '/eos/cms/store/' in fname: prepath = 'root://eoscms.cern.ch//'
            fname = prepath+fname
            ## fname is now xrootd compatible
            ROOT.gEnv.SetValue("TFile.AsyncReading", 1);
            
            ## open the file
            #f = ROOT.TFile.Open(fname);
            f = ROOT.TFile.Open(fname+"?readaheadsz=65535");

            ## get the counters from the trees (even after a skim, so just get number of entries)
            tree_events = f.Get('Events')
            entries     = tree_events.GetEntries()
            tree_runs   = f.Get('Runs')            
            sample_nevt = entries
            if not data: 
                for i,event in enumerate(tree_runs):
                    sample_nevt = event.genEventSumw
                    if i: 
                        break
            ## done. but as of now sample_nevt is not used...
            chunk = options.chunkSize
            if entries < chunk:
                print "  ",fnameWithPath,("  DATA" if data else "  MC")," single chunk"
                jobs.append((dataset,
                             fname,
                             sample_nevt,
                             "_Friend_%s"%shortfname.replace(".root",""),
                             data,
                             xrange(entries),
                             -1))
            else:
                nchunk = int(math.ceil(entries/float(chunk)))
                print "  ",fnameWithPath,("  DATA" if data else "  MC")," %d chunks" % nchunk
                for i in xrange(nchunk):
                    if options.chunks != []:
                        if i not in options.chunks: continue
                    r = xrange(int(i*chunk),min(int((i+1)*chunk),entries))
                    jobs.append((dataset,
                                 fname,
                                 sample_nevt,
                                 "_Friend_%s.chunk%d" % (shortfname.replace(".root",""),i),
                                 data,
                                 r,
                                 i))

    print "\n"
    print "I have %d taks to process" % len(jobs)

    print 'I\'m using the following list of modules',options.moduleList
    imports = globals()[options.moduleList] + options.imports
    if options.submit:
        import os, sys

        writelog = ""
        logdir   = ""
        if options.logdir: 
            logdir = os.environ['PWD']+'/'+options.logdir
            logdir = logdir.rstrip("/")
            if not os.path.exists(logdir):
                os.system("mkdir -p "+logdir)

        runner = ""
        super = ""
        if options.env == "lxbatch":
            runner = options.runner
            super  = "bsub -q {queue}".format(queue = options.queue)
        elif  options.env == "condor":
            runner = options.runner
        else:
            print "ERROR. Scheduler ",options.env," not implemented. Choose either 'lsf' or 'condor'."
            sys.exit(1)

        basecmd_tmp = " {self} -N {chunkSize} -t {tree} --signals {signals} --moduleList {moduleList} ".format(self=sys.argv[0], chunkSize=options.chunkSize, tree=options.tree, signals=options.signals, moduleList=options.moduleList)

        friendPost = ""
        if options.friend: 
            friendPost += " --friend " 
        cmds = []
        # name identifies the sample, it is the same for all files, at least at the moment
        # one wants a single condor submit file for each sample, managing all jobs
        firstname = ""
        nameChanged = False
        for (name,fin,sample_nevt,fout,data,_range,chunk) in jobs:
            
            if firstname == "":
                firstname = name
                nameChanged = True
            elif firstname != name:
                firstname = name
                nameChanged = True
            else:
                nameChanged = False

            finNoPrepath = fin.replace("root://eosuser.cern.ch//","").replace("root://eoscms.cern.ch//","")
            basecmd = basecmd_tmp + " {data} {output}".format(data=finNoPrepath, output=outdir+"/"+dataset+"/")

            if not chunk or chunk == -1:
                if nameChanged:
                    condorSubFile = writeCondorCfg(logdir,name,maxRunTime=options.runtime,maxMemory=options.memory)
            #if chunk != -1:
            if options.env == 'condor':
                tmp_f = open(condorSubFile, 'a')
                tmp_f.write('arguments = {base} -d {data} {chunk} {post} \nqueue 1 \n\n'.format(base=basecmd, data=name, chunk='' if chunk == -1 else '-c '+str(chunk), post=friendPost))
                tmp_f.close()
                if chunk <= 0:
                    if nameChanged:
                        # the command is here, it just tracks the condor submit file name
                        # we might still be adding arguments in the condor submit file itself
                        cmds.append('condor_submit ' + condorSubFile)
            else:
                print 'option --env MUST be condor (which it is by default)'
        for cmd in cmds:
            print 'i will now run this if not --pretend'
            print cmd
            if not options.pretend:
                os.system(cmd)

        exit()

    maintimer = ROOT.TStopwatch()
    def _runIt(myargs):
        (dataset,fin,sample_nevt,fout,data,_range,chunk) = myargs
        modules = []
        for mod, names in imports: 
            import_module(mod)
            obj = sys.modules[mod]
            selnames = names.split(",")
            for name in dir(obj):
                if name[0] == "_": continue
                if name in selnames:
                    print "Modules to run: ",options.modules,"  name = ",name
                    if len(options.modules) and name not in options.modules: continue
                    print "Loading %s from %s " % (name, mod)
                    print "Running on dataset = ",dataset
                    signal = any(x in dataset for x in options.signals.split(','))
                    if name=='genQEDJets' and not signal: continue
                    ## the kamuca
                    if name=='kamucaSyst' and 'SingleMuon' in dataset: continue
                    modules.append(getattr(obj,name)())
        if options.noOut:
            if len(modules) == 0: 
                raise RuntimeError("Running with --noout and no modules does nothing!")
        ppargs=[fin]+args
        p=PostProcessor(outdir+"/"+dataset+"/",ppargs,options.cut,options.branchsel,modules,options.compression,options.friend,fout,options.json,options.noOut,options.justcount,_range,treeName="Events",prefixRootFileName="tree")
        p.run()

    print 'this is jobs', jobs
    if options.jobs > 1:
        from multiprocessing import Pool
        pool = Pool(options.jobs)
        pool.map(_runIt, jobs) if options.jobs > 0 else [_runIt(j) for j in jobs]
    else:
        ret = dict(map(_runIt, jobs))
