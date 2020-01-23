#!/usr/bin/env python
# USAGE: 
# 1. if you have already the output of friendChunkCheck.sh -z fdir > zombies.txt
# ./scripts/friendChunkResub.py frienddir maintreedir zombies.txt -N 500000
# 2. if you don't have it
#  ./scripts/friendChunkResub.py frienddir maintreedir --run-checker -N 500000

import sys,os,re

def writeCondorCfg(logdir, name, flavour=None, maxRunTime=None):

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
request_memory = 4000
+MaxRuntime = 7200 \n
""".format(ld=logdir,
           name=name,
           here=os.environ['PWD'])
    if os.environ['USER'] in ['mdunser', 'psilva']:
        job_desc += '+AccountingGroup = "group_u_CMST3.all"\n'
    ##job_desc += 'queue 1\n'

    jobdesc_file = logdir+'/'+name+'.condor'
    with open(jobdesc_file,'w') as outputfile:
        outputfile.write(job_desc)
        outputfile.write('\n\n')
        outputfile.close()
    return jobdesc_file

if __name__ == "__main__":
    
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] dir_with_friend_trees dir_with_trees [checkedfile.txt]")
    parser.add_option("-c",      "--run-checker", dest="runChecker",  action='store_true', default=False, help="Run the script friendChunkCheck.sh");
    parser.add_option("-N", "--events",  dest="chunkSize", type="int",    default=1000000, help="Default chunk size when splitting trees");
    (options, args) = parser.parse_args()

    fdir = args[0]
    maindir = args[1]

    if options.runChecker:
        print "Running friendChunkCheck.py now on ",fdir,". Will take time..."
        if len(args)==3:
            print "Can't run the checker if you pass the output file of a previous check (for safety)"
            sys.exit(1)
        pyCheckScript = '{cmssw}/src/CMGTools/WMass/python/postprocessing/scripts/friendChunkCheck.py -z {frienddir} > tmpcheck.txt'.format(cmssw=os.environ['CMSSW_BASE'],
                                                                                                                                              frienddir=fdir)
        if os.path.isfile('tmpcheck.txt'):
            print "File tmpcheck.txt exists. It means you could be running friendChunkCheck.sh. If not, remove it, and run it again."
            sys.exit()
        else:
            os.system(pyCheckScript)
        print "Done. The list of files is in tmpcheck.txt."

    tmpfile = 'tmpcheck.txt' if len(args)<3 else args[2]
    txtfile=open(tmpfile,'r')
    cmds = []
    for line in txtfile:
        l = line.rstrip()
        if l.startswith('#') or l.startswith('DONE'): continue
        if not l.endswith('OK'):
            base = os.path.basename(l)
            tokens = base.split('.')
            dataset = '_'.join(tokens[0].split('_')[2:])
            chunk = tokens[1].split('chunk')[-1]
            #print "# resubmitting dataset = ",dataset," chunk ",chunk
            cmd = "python postproc_batch.py {maintreedir} {frienddir} --friend -N {nevents} -d {dataset} -c {chunk} ".format(maintreedir=maindir,frienddir=fdir,nevents=options.chunkSize,dataset=dataset,chunk=chunk)
            cmds.append(cmd)
            #print cmd
    f_name   = writeCondorCfg('friends_log', 'friendsResubmission')
    f_condor = open(f_name, 'a')
    for cmd in cmds:
        f_condor.write(cmd.replace('python', 'arguments = ')+' \n')
        f_condor.write('queue 1 \n\n')
    f_condor.close()

    print 'run this on now:'
    print 'condor_submit '+f_name


