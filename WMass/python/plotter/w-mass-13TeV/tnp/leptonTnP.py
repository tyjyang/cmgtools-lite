import ROOT, os
ROOT.gROOT.SetBatch(True)

## USAGE:
## python leptonTnP.py -i /eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_recoLeptons_and_triggerMatch_latest/ -c mu -o outdir -a selection

## option -o indir : saves the output trees in the input directory. if none given it saves it where it's run
## option -c       : has to be mu or el
## option -a       : either "trigger" or "selection"


## new
## python leptonTnP.py -i /eos/cms/store/data/Run2016G/SingleMuon/NANOAOD/Nano02Dec2019-v1/ -o /eos/cms/store/cmst3/group/wmass/w-mass-13TeV/tnp/2020-10-28-firstTest/RunG/ -b


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] mc.txt cuts.txt treeDir outputDirSkims ")
    parser.add_option('-c', '--channel'   , dest='channel'   , type='string'      , default='mu',  help='run tnp ntuples for muons/electrons. default mu')
    parser.add_option('-i', '--indir'     , dest='indir'     , type='string'      , default=''  ,  help='directory with the trees')
    parser.add_option('-m', '--maxentries', type='int'       , default=-1         , help='run only these many entries.')
    parser.add_option('-o', '--outdir'    , type='string'    , default='./'       , help='save the files in this outdir')
    parser.add_option('-a', '--analyzer'  , dest='analyzer'  , type='string'      , default='selection',  help='run tnp ntuples for trigger/selection. default selection')
    parser.add_option('-b', '--batch'     , default=False    , action='store_true', help='use the batch system instead (requires setting also -n)')
    parser.add_option('-n', '--nfiles'    , type='int'       , default=10         , help='set the number of files to be run in each batch job')
    parser.add_option('-f', '--files'     , type='string'    , default=''         , help='the list of files for running on the batch')
    parser.add_option('-s', '--special'   , type='string'    , default=''         , help='some special string to go into the filename and so on')
    (options, args) = parser.parse_args() 


    if options.outdir == 'indir':
        options.outdir = options.indir

    ## make the outdir in the script... easier
    os.system('mkdir -p '+options.outdir)

    xrdindir  = options.indir
    xrdoutdir = options.outdir
    ## use xrootd
    if '/eos/cms/store/' in xrdindir and not 'eoscms' in xrdindir:
        xrdindir = 'root://eoscms.cern.ch/'
    if '/eos/user/' in xrdindir and not 'eosuser' in xrdindir:
        xrdindir = 'root://eosuser.cern.ch/'

    ## use xrootd also for output directory and file
    if '/eos/cms/store/' in xrdoutdir and not 'eoscms' in xrdoutdir:
        xrdoutdir = 'root://eoscms.cern.ch/'+xrdoutdir
    if '/eos/user/' in xrdoutdir and not 'eosuser' in xrdoutdir:
        xrdoutdir = 'root://eosuser.cern.ch/'+xrdoutdir


    if options.batch:
        listoffiles = []
        for root, dirnames, filenames in os.walk(options.indir):
            for filename in filenames:
                if '.root' in filename:
                    listoffiles.append(xrdindir+os.path.join(root, filename))

        listoffilechunks = []

        for ff in range(len(listoffiles)/options.nfiles+1):
            listoffilechunks.append(listoffiles[ff*options.nfiles:ff*options.nfiles+options.nfiles])

        #for i in  listoffilechunks:
        #    print i


    ###for sd in os.listdir(options.indir):
    ##for fi,infile in enumerate(listoffiles):
    ##    # if not '00336_246' in sd: continue
    ##    treefile = infile #'/'.join([options.indir,sd])

    ##    if fi : continue

    ##    #if not os.path.isfile(treefile): 
    ##    #    print 'treefile not found!', treefile
    ##    #    continue
    ##    #print treefile

    ##    if (options.channel == 'el' and 'SingleMu' in treefile): continue
    ##    if (options.channel == 'mu' and 'SingleEl' in treefile): continue
    ##    ##if options.analyzer == 'selection' and 'DYJets' in sd: continue
    ##    ##if 'DYJets' not in sd: continue

    ##    print 'runing on file', treefile

    ##    tmp_file = ROOT.TFile(treefile, 'read')
    ##    tmp_tree = tmp_file.Get('Events')
    ##    
    ##    ## make the instance of the worker
    ##    ##if 'el' in options.channel or 'mu' in options.channel:
    ##    ## friends now also for muons!
    ##    ## ffile = xrdindir+'/friendsNEWSCALE/tree_Friend_'+sd+'.root'
    ##    ## tmp_friend_file = ROOT.TFile(ffile)
    ##    ## tmp_friend_tree = tmp_friend_file.Get('Friends')
    ##    ##else: 
    ##    ##    tmp_friend_tree = None
    ##

    if options.batch:

        dm = 'mc' if 'DYJ' in options.indir else 'data'
        runperiod = ''
        if dm == 'data':
            runperiod = '_'+[i for i in options.indir.split('/') if 'Run' in i][0]

        tmp_condor_filename = 'logs/condor_submit_{an}_{ch}_{dm}{rp}{s}.condor'.format(an=options.analyzer,ch=options.channel,dm=dm,rp=runperiod,s=options.special)
        tmp_condor = open(tmp_condor_filename,'w')
        tmp_condor.write('''Executable = dummy_exec.sh
use_x509userproxy = true
getenv      = True

environment = "LS_SUBCWD={here}"
request_memory = 2000
+MaxRuntime = 15000 \n\n'''.format(here=os.environ['PWD']))
        for il,fs in enumerate(listoffilechunks):
            if not len(fs): continue
            tmp_condor.write('arguments = leptonTnP.py -a {a} -c {ch} -f {files} -o {outdir} \n'.format(files=','.join(fs), outdir=options.outdir, a=options.analyzer,ch=options.channel))
            tmp_condor.write('''
Log        = logs/log_condor_{c}_{dm}{rp}{s}_chunk{ch}.log
Output     = logs/log_condor_{c}_{dm}{rp}{s}_chunk{ch}.out
Error      = logs/log_condor_{c}_{dm}{rp}{s}_chunk{ch}.error\n'''.format(ch=il,c=options.channel,dm=dm,rp=runperiod,s=options.special))
            tmp_condor.write('queue 1\n\n')
        tmp_condor.close()

        print 'wrote condor file', tmp_condor_filename
 
        ##print 'this will result in a total of {n} jobs'.format(n=totalNumberOfJobs)
        ##print 'submitting to condor...'
        ##if not options.pretend:
        ##    os.system('mkdir -p logs')
        ##for sc in condorSubmitCommands:
        ##    print sc
        ##    if not options.pretend:
        ##        os.system(sc)
        ##print '...done'


    else:
        chain = ROOT.TChain('Events')
        for i in options.files.split(','):
            chain.Add(i)

        if options.analyzer=='trigger'     : 
            print 'compiling trigger ntuple producer'
            ROOT.gROOT.ProcessLine(".L TnPNtuplesTriggerEfficiency.C+" )
            tnp_worker = ROOT.TnPNtuplesTriggerEfficiency(chain)#,tmp_friend_tree)
        elif options.analyzer=='selection' : 
            print 'compiling selection ntuple producer'
            ROOT.gROOT.ProcessLine(".L TnPNtuplesSelectionEfficiency.C+" )
            tnp_worker = ROOT.TnPNtuplesSelectionEfficiency(chain) #tmp_tree)#,tmp_friend_tree)
        else: 
            print "analyzer can be either trigger/selection."
            exit(1)

        if 'mu' in options.channel:
            tnp_worker.setFlavor(13)
        else: 
            tnp_worker.setFlavor(11)

        tnp_worker.makeMapFromJson("/afs/cern.ch/work/m/mdunser/public/cmssw/w-mass-13TeV/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/w-mass-13TeV/tnp/formattedJSON.json")

        refstring = options.files.split(',')[0].split('/')[-1].replace('.root','')
        dm = 'mc' if 'DYJ' in options.files.split(',')[0] else 'data'

        tnp_worker.setOutfile(xrdoutdir+'/'+options.analyzer+'TnP_'+refstring+'_{dm}_{ch}.root'.format(ch=options.channel,dm=dm))
        tnp_worker.Loop(options.maxentries)
    
        del tnp_worker ## this segfaults otherwise

