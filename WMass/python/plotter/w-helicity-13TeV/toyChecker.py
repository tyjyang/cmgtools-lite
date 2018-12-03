#!/usr/bin/env python
# python w-helicity-13TeV/toyChecker.py toys/diffXsec_el_2018_09_20_group10_legacySF/ -c -q cmscaf1nd > testResubToys.txt
# then
# cat testResubToys.txt | grep bsub | bash

import os.path
import sys,ROOT,os
ROOT.gROOT.SetBatch(True)

class CardsChecker:
    def __init__(self, toy_dir, options):
        self.toy_dir = toy_dir
        self.options = options
        self.sourceFiles = {}
        self.pycmd = {}
        retrydirs = [d for d in os.listdir(toy_dir) if 'retry' in d]
        itry = len(retrydirs)+1
        self.resub_card_dir='{cdir}/retry_{i}'.format(cdir=toy_dir,i=itry)
        os.mkdir(self.resub_card_dir)
        for f in os.listdir(toy_dir+'/jobs/'):
            if not f.endswith('.sh'): 
                continue
            if "dummy_exec" in f: continue
            tmp_f = open(toy_dir+'/jobs/'+f, 'r')
            lines = tmp_f.readlines()
            for i in lines:
                if 'combinetf.py' in i: pycmds = i
            seed = pycmds.split(' --seed ')[1].split()[0]
            self.sourceFiles[seed] = f
        print '## Expecting {n} root files'.format(n=len(self.sourceFiles))

    def makeCondorFile(self, srcFiles):
        jobdir = self.resub_card_dir+'/jobs/'
        os.mkdir(jobdir)
        dummy_exec = open(jobdir+'/dummy_exec.sh','w')
        dummy_exec.write('#!/bin/bash\n')
        dummy_exec.write('bash $*\n')
        dummy_exec.close()
     
        condor_file_name = jobdir+'/condor_submit.condor'
        condor_file = open(condor_file_name,'w')
        # creating output directories for logs
        os.system("mkdir -p {jd}/logs/".format(jd=self.resub_card_dir))
        os.system("mkdir -p {jd}/errs/".format(jd=self.resub_card_dir))
        os.system("mkdir -p {jd}/outs/".format(jd=self.resub_card_dir))
        logdir = self.resub_card_dir + "/logs/"
        outdirCondor = self.resub_card_dir + "/outs/"
        errdir = self.resub_card_dir + "/errs/"
        condor_file.write('''Universe = vanilla
Executable = {de}
use_x509userproxy = $ENV(X509_USER_PROXY)
Log        = {ld}/$(ProcId).log   
Output     = {od}/$(ProcId).out   
Error      = {ed}/$(ProcId).error
getenv      = True
environment = "LS_SUBCWD={here}"
next_job_start_delay = 1
request_memory = 4000
+MaxRuntime = {rt}\n
'''.format(de=os.path.abspath(dummy_exec.name), ld=os.path.abspath(logdir), od=os.path.abspath(outdirCondor), ed=os.path.abspath(errdir), 
           rt=int(options.runtime*3600), here=os.environ['PWD'] ) )

        if os.environ['USER'] in ['mdunser', 'psilva']:
            condor_file.write('+AccountingGroup = "group_u_CMST3.all"\n\n\n')
        for sf in srcFiles:
            condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(sf)))
        condor_file.close()
        return condor_file_name


    def checkToys(self):
        resubcmds = {}
        srcfiles = [] 
        for seed,source in self.sourceFiles.iteritems():
            f = "fitresults_" + str(seed) + ".root"
            if not os.path.exists(self.toy_dir+'/'+f) or os.path.getsize(self.toy_dir+'/'+f) < 1000.:
                if self.options.verbose>1: print '# file ',f,'  not present (or size < 1000) in ',self.toy_dir
                tmp_job_path = self.toy_dir+'/jobs/'+source
                srcfiles.append(os.path.abspath(tmp_job_path))                

        if len(srcfiles) == 0: return 0
        else:
            cf = self.makeCondorFile(srcfiles)
            resubcmds[0] = 'condor_submit {rf} '.format(rf = cf)
            return resubcmds

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage='%prog toydir')
    parser.add_option('-c', '--check-toys', dest='checkToys', default=False, action='store_true', help='Check if there are all the datacards and ROOT files');
    parser.add_option('-q', '--queue', dest='queue', type='string', default='1nd', help='choose the queue to submit batch jobs (default is 8nh)');
    parser.add_option('-v', '--verbose', dest='verbose', default=0, type=int, help='Degree of verbosity (0=default prints only the resubmit commands)');
    parser.add_option('-r', '--runtime', default=24, type=int,  help='New runtime for condor resubmission in hours. default None: will take the original one.');
    (options, args) = parser.parse_args()

    if options.checkToys:
        toydir = args[0]
        if len(args)<1: print 'needed inputs: toys_dir '; quit()
        cc = CardsChecker(toydir, options)
        result = cc.checkToys()
        if result == 0: print 'All files are GOOD.'
        else: 
            keys = result.keys()
            keys.sort()
            for k in keys: print result[k]
            print '## in total have to resubmit {n} jobs'.format(n=len(result))

