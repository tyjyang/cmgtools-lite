#!/bin/env python

# usage: python submitToys.py cards_el/Wel_card_withXsecMask.hdf5 10000 -n 20 --outdir output -q 8nh
#
# python w-helicity-13TeV/submitToys.py cards/diffXsec_2018_06_29_group10_absGenEta_moreEtaPtBin/Wel_plus_card_withXsecMask.meta 11000 -n 2 --outdir toys/diffXsec_2018_06_29_group10_absGenEta_moreEtaPtBin_newGenXsec/ -t 1 -q cmscaf1nd

jobstring  = '''#!/bin/sh
ulimit -c 0 -S
ulimit -c 0 -H
set -e
cd CMSSWBASE
export SCRAM_ARCH=slc6_amd64_gcc530
eval `scramv1 runtime -sh`
cd OUTDIR
COMBINESTRING

'''


jobstring_tf = '''#!/bin/sh
ulimit -c 0 -S
ulimit -c 0 -H
set -e
export SCRAM_ARCH=slc6_amd64_gcc700
cd CMSSWBASE
eval `scramv1 runtime -sh`
cd OUTDIR
COMBINESTRING

'''


import ROOT, random, array, os, sys

if __name__ == "__main__":
    
    from optparse import OptionParser
    parser = OptionParser(usage='%prog workspace ntoys [prefix] [options] ')
    parser.add_option('-n'  , '--ntoy-per-job'  , dest='nTj'           , type=int           , default=None , help='split jobs with ntoys per batch job')
    parser.add_option('-t'  , '--threads'       , dest='nThreads'      , type=int           , default=1    , help='use nThreads in the fit (suggested 2 for single charge, 1 for combination)')
    parser.add_option(        '--dry-run'       , dest='dryRun'        , action='store_true', default=False, help='Do not run the job, only print the command');
    parser.add_option('-q'  , '--queue'         , dest="queue"         , type="string"      , default="1nd", help="Select the queue to use");
    parser.add_option('--outdir', dest='outdir', type="string", default=None, help='outdirectory');
    (options, args) = parser.parse_args()

    ## for tensorflow the ws has to be the datacard!
    
    workspace = args[0]; wsbase = os.path.basename(workspace).split('.')[0]
    ntoys = int(args[1])
    prefix = args[2] if len(args)>2 else wsbase
    charge = 'plus' if 'plus' in wsbase else 'minus'
    floatPOIs = 'XsecMask' not in workspace

    print "Submitting {nt} toys with workspace {ws} and prefix {pfx}...".format(nt=ntoys,ws=workspace,pfx=prefix)

    absopath  = os.path.abspath(options.outdir)
    if not options.outdir:
        raise RuntimeError, 'ERROR: give at least an output directory. there will be a YUGE number of jobs!'
    else:
        if not os.path.isdir(absopath):
            print 'making a directory and running in it'
            os.system('mkdir -p {od}'.format(od=absopath))

    jobdir = absopath+'/jobs/'
    if not os.path.isdir(jobdir):
        os.system('mkdir {od}'.format(od=jobdir))

    random.seed()

    for j in xrange(int(ntoys/int(options.nTj))):
        ## make new file for evert parameter and point
        job_file_name = jobdir+'/job_{j}_toy{n:.0f}To{nn:.0f}.sh'.format(j=j,n=j*int(options.nTj),nn=(j+1)*int(options.nTj))
        tmp_file = open(job_file_name, 'w')

        tmp_filecont = jobstring_tf
        #cmd = 'text2tf.py -t {n} --seed {j}{jn} {dc}'.format(n=int(options.nTj),dc=os.path.abspath(workspace),j=j*int(options.nTj)+1,jn=(j+1)*int(options.nTj)+1)
        cmd = 'combinetf.py -t {n} --seed {j}{jn} {dc} --nThreads {nthr}'.format(n=int(options.nTj),dc=os.path.abspath(workspace),j=j*int(options.nTj)+1,jn=(j+1)*int(options.nTj)+1,nthr=options.nThreads)
        if floatPOIs: cmd += ' --POIMode none '
        tmp_filecont = tmp_filecont.replace('COMBINESTRING', cmd)
        tmp_filecont = tmp_filecont.replace('CMSSWBASE', os.environ['CMSSW_BASE']+'/src/')
        tmp_filecont = tmp_filecont.replace('OUTDIR', absopath+'/')
        tmp_file.write(tmp_filecont)
        tmp_file.close()
        os.system('chmod u+x {f}'.format(f=job_file_name))
        cmd = 'bsub -o {log} -q {queue} {job}'.format(log=job_file_name.replace('.sh','.log'),queue=options.queue,job=job_file_name)
        if options.dryRun:
            print cmd
        else:
            os.system(cmd)
        
    sys.exit()

