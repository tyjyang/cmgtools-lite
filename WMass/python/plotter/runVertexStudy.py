#!/usr/bin/env python

import os, subprocess, re
import ROOT

#workingPoints = ["alwaystrue", "genEtaPt", "vertexPresel", "muonInAccept", "muMediumId", "muTightIso", "mtl1pf40", "trigger"]
workingPoints = ["alwaystrue", "genMuNoEtaPt", "vertexPresel", "muonInAccept", "muMediumId", "muTightIso", "mtl1pf40", "trigger"]


path = "{cmssw}/src/CMGTools/WMass/python/plotter/".format(cmssw=os.environ['CMSSW_BASE'])
mcafile  = path + "w-mass-13TeV/wmass_mu/mca_wmu_forTest.txt"
cutfile  = path + "w-mass-13TeV/wmass_mu/cuts_wmu_VertexStudy.txt"
plotfile = path + "w-mass-13TeV/wmass_mu/plots_forTest.txt"
proc = "QCD" # Wnopt
# 
#plots = "dzVertex_gen_primary__Wpt,dzVertex_gen_primary__dressedLepPt"

# this script is meant to run the study on vertex using W MC or other MC

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options]")
    parser.add_option("-d", "--dry-run", dest="dryRun", default=False, action="store_true",  help="Just print commands") 
    parser.add_option("--log", "--log-dir", dest="logdir", type="string", default=None, help="Directory of stdout and stderr");
    parser.add_option("-o", "--outdir", dest="outdir", type="string", default=None, help="Output folder for all plots (each one will be in a specific subfolder)");
    parser.add_option("-n", "--job-name", dest="jobName",   type="string", default="vertexStudy", help="Name assigned to jobs");
    parser.add_option("--ntuples-path",   dest="ntuplesPath", type="string", default="/eos/cms/store/cmst3/group/wmass/mciprian/heppyNtuples/TREES_W_VERTEXSTUDY_94X_V2/", help="Path to ntuples")
    (options, args) = parser.parse_args() 

    ntuplesPath = options.ntuplesPath

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

    if not options.outdir: 
        print 'ERROR: must give a directory where to store output!\nexiting...'
        exit()
    else:
        outdir = options.outdir.rstrip("/") + "/"
        if not os.path.exists(outdir):
            os.system("mkdir -p "+outdir)


        condor_fn = options.logdir+'/condor_vertexStudy.condor'
        condor_f = open(condor_fn,'w')
        condor_f.write('''Universe = vanilla
Executable = {runner}
use_x509userproxy = true
Log        = {ld}/$(ProcId).log
Output     = {ld}/$(ProcId).out
Error      = {ld}/$(ProcId).error
getenv      = True
request_memory = 2000
+MaxRuntime = 7200
+JobBatchName = "{name}"\n
'''.format(runner=runner,ld=options.logdir,name=options.jobName))
        if os.environ['USER'] in ['mdunser', 'psilva']:
            condor_f.write('+AccountingGroup = "group_u_CMST3.all"\n')
        elif os.environ['USER'] in ['mciprian']:
            condor_f.write('+AccountingGroup = "group_u_CMS.CAF.ALCA"\n')
        condor_f.write('\n')


        for wp in workingPoints:
            cmd = "python mcPlots.py "
            cmd += " {mca} {cut} {plot}".format(mca=mcafile,cut=cutfile,plot=plotfile)
            cmd += " -f -l 35.9 --s2v --tree treeProducerWMass --obj tree --noCms -j 8 "
            cmd += " --legendFontSize 0.05 --setLegendCoordinates 0.2,0.77,0.9,0.92 --allProcInLegend --noLegendRatioPlot --n-column-legend 2 "            
            cmd += " --sP {plots} ".format(plots=plots)
            cmd += " -P {ntp} -F Friends {ntp}friends/tree_Friend_{{cname}}.root ".format(ntp=ntuplesPath)
            cmd += " -p {proc} -W 1.0 -U {wp} --pdir {o}/{wp}/ ".format(proc=proc,o=outdir,wp=wp)        
            cmdargs = [x.strip() for x in cmd.split()]
            strargs=''
            for a in cmdargs: # join do not preserve " or '
                ## escaping not needed for condor if '{P}' in a or '*' in a: 
                ## escaping not needed for condor     a = '\''+a+'\''
                strargs += ' '+a+' '

            condor_f.write('\narguments = {d} {cmssw} {cmdargs}\n'.format(d=os.getcwd(), cmssw = os.environ['CMSSW_BASE'], cmdargs=strargs))
            condor_f.write('queue 1\n\n')

        condor_f.close()
        cmd = 'condor_submit {sf}'.format(sf=condor_fn)

        if options.dryRun: 
            print cmd
        else:
            os.system(cmd)
