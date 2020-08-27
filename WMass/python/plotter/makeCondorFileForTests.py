import array
import json
import ROOT
import re 
import os

# I wanted to add python commands directly in the condor submit file, but it doesn't work easily
# following comments are obsolete, as now I put python commands in a bash script and run condor on it
# caution with arguments requiring quotes. They need special formatting to be fed into condor arguments
# see e.g. https://github.com/htcondor/htcondor/blob/master/src/condor_utils/condor_arglist.h
# the best thing is to use ' whenever an argument need to be put inside quotation marks
# for options with multiple arguments one must ...

# ideally one would have a link in python/plotter pointing to a folder on cern web page
plotterPath = "{cms}/src/CMGTools/WMass/python/plotter".format(cms=os.environ['CMSSW_BASE']) 
outfolderPlots = plotterPath + "/plots/Wlike/TREE_4_WLIKE_MU/testIsoQCD/withMtCut_rebin_fillstyle_iso0p15/"  # will add an option
ntuples = "/afs/cern.ch/work/m/mdunser/public/wmassTrees/SKIMS_muons_latest/"
mca   = "w-mass-13TeV/wmass_mu/mca_testIsoQCD.txt"
cuts  = "w-mass-13TeV/wmass_mu/cuts_wmu_test.txt"
plots = "w-mass-13TeV/wmass_mu/plots_forMuonFR.txt"

absetaBinning = [0.0, 1.2, 2.4]
netaBins = len(absetaBinning) - 1

##
## do not add --pdir here, it is done inside the script, just define the name above
##
basecmd = "python {pp}/mcPlots.py  {pp}/{mca} {pp}/{cuts} {pp}/{plots}".format(pp=plotterPath,
                                                                               mca=mca,
                                                                               cuts=cuts,
                                                                               plots=plots)
# foo
basecmd += " -f -l 35.9 --s2v --tree treeProducerWMass --obj tree --lspam '#bf{CMS} #it{Preliminary}' --noCms -j 8  --setTitleXoffset 0.95 "
# processes
#basecmd += " -p 'data_isoMore0p15,data_isoMore0p25,QCD_isoLess0p15,QCD_isoMore0p15,QCD_isoMore0p25' "
basecmd += " -p 'data_isoMore0p15,QCD_isoMore0p15' "
# plotting mode and axis settings 
basecmd += " --plotmode norm --contentAxisTitle 'arbitrary units'  "
# legend
basecmd += " --n-column-legend 2 --setLegendCoordinates 0.2,0.75,0.85,0.92 --legendFontSize 0.038 --allProcInLegend --noLegendRatioPlot "
# plots
basecmd += " --sP 'ptl1coarse,etal1coarse,isol1,mtcoarse' " 
# ntuples
basecmd += " -P {ntp}/ -F Friends {ntp}/friends/tree_Friend_{{cname}}.root ".format(ntp=ntuples)
# event weight
basecmd += " -W 'puw2016_nTrueInt_36fb(nTrueInt)*_get_muonSF_recoToSelection(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta)*_get_muonSF_selectionToTrigger(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta,LepGood1_charge)*prefireJetsWeight(LepGood1_eta)' "
# additional cuts
basecmd += " -X muTightIso "
# ratio options
basecmd += " --showRatio --maxRatioRange 0.6 1.4 --fixRatioRange --ratioDen QCD_isoMore0p15 --ratioNums 'data_isoMore0p15' --ratioYLabel 'data/QCD(iso>0.15)' --ratioYLabelSize 0.11 "
# for tests
#basecmd += " --max-entries 10 "

## other options, such as specific cuts for binning or output folder are managed below


if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    # parser.add_option('-i','--indir',     dest='indir',     default='',   type='string', help='Folder with the condor file (the one with all the logs)')
    parser.add_option('-o','--outdir',     dest='outdir',     default='',   type='string', help='outdput directory to write the new condor file')
    parser.add_option('-f','--fname',     dest='fname',     default='condor_test',   type='string', help='Name of condor file to be run (.condor extension optional, added automatically if not passed)')
    parser.add_option('-n','--job-name',     dest='jobname',     default='testPlots',   type='string', help='Name assigned to condor jobs.')
    # parser.add_option('--ntuples-dir-in',     dest='ntuplesDirIn',     default='',   type='string', help='Path to input unskimmed ntuples')
    # parser.add_option('--ntuples-dir-out',     dest='ntuplesDirOut',     default='',   type='string', help='Path to output skimmed ntuples')
    parser.add_option('-d', '--dry-run', dest='dryrun' , default=False , action='store_true',   help='Print commands, but do not run them')
    parser.add_option('-l', '--run-local', dest='runlocal' , default=False , action='store_true',   help='Run commands locally, not submitting jobs (not suggested, unless condor does not work)')
    # parser.add_option('-f', '--friends', dest='friends' , default=False , action='store_true',   help='Specify if doing friend trees, things need to be handled differently (paths in --ntuples-dir-in and --ntuples-dir-out do no need to end with "/friends/"')
    (options, args) = parser.parse_args()

    # create base folder
    if not os.path.exists("jobsLogCondor/"):
        os.system('mkdir -p jobsLogCondor/')
        dummyfile = open("jobsLogCondor/dummy_exec.sh","w")
        dummyfile.write("#! /bin/bash\n\n")
        dummyfile.write("bash $*\n")
        dummyfile.close()

    if not options.outdir:
        print "Warning: you must specified an output folder with option -o." 
        print "It will be created inside jobsLogCondor/. Abort"
        quit()

    outdir = "jobsLogCondor/" + options.outdir 
    if not outdir.endswith("/"): 
        outdir +=  "/"
    os.system('mkdir -p {od}'.format(od=outdir))
    os.system('mkdir -p {od}src/'.format(od=outdir))
    os.system('mkdir -p {od}logs/'.format(od=outdir))    

    fname = options.fname
    if not fname.endswith(".condor"):
        fname += ".condor"
    outf = outdir + fname
    of = open(outf,"w")

    ###
    of.write("Universe = vanilla\n")
    of.write("Executable = {pp}/jobsLogCondor/dummy_exec.sh\n".format(pp=plotterPath))
    of.write("use_x509userproxy = true\n")
    #of.write("x509userproxy = $ENV(X509_USER_PROXY)\n")
    of.write("Log        = {od}/logs/$(ProcId).log\n".format(od=outdir))
    of.write("Output     = {od}/logs/$(ProcId).out\n".format(od=outdir))
    of.write("Error      = {od}/logs/$(ProcId).error\n".format(od=outdir))
    of.write("getenv      = True\n")
    of.write('environment = "LS_SUBCWD={here}"\n'.format(here=os.environ['PWD']))
    of.write("request_memory = 4000\n")
    of.write("+MaxRuntime = 86400\n")
    of.write('+JobBatchName = "{n}"\n'.format(n=options.jobname))
    if os.environ['USER'] in ['mdunser', 'psilva']:
        of.write('+AccountingGroup = "group_u_CMST3.all"\n')
    elif os.environ['USER'] in ['mciprian']:
        of.write('+AccountingGroup = "group_u_CMS.CAF.ALCA"\n')
    of.write('\n')

    for ivar in range(netaBins):

        cmd = basecmd
        varl = absetaBinning[ivar]
        varh = absetaBinning[ivar+1]
        varcut = " -A onelep varcut 'abs(LepGood1_eta)>{varl} && abs(LepGood1_eta)<{varh}' ".format(varl=str(varl),
                                                                                                    varh=str(varh))
        varFolder = "abseta_{vl}_{vh}".format(vl=str(varl).replace(".","p"),
                                              vh=str(varh).replace(".","p"))
        out_opt = " --pdir {op}/{vf}/ ".format(op=outfolderPlots,
                                               vf=varFolder) 

        cmd = cmd + varcut + out_opt
        

        # if running locally we have to setup the environment
        if options.runlocal:
            os.system("cd {pp}/".format(pp=plotterPath))
            os.system("eval `scramv1 runtime -sh`")
            print "="*30
            print "Option --run-local is active: running following command in local"            
            print "-"*30
            print cmd
            print "="*30
            if options.dryrun: 
                print "Option -d active, not running commands"
            else:
                os.system(cmd)

        # open .sh file to dump the python command
        srcfileName = "{od}src/{vf}.sh".format(od=outdir,vf=varFolder)
        srcfile = open(srcfileName,"w")
        srcfile.write("#! /bin/bash\n\n")
        srcfile.write("cd {pp}/\n".format(pp=plotterPath))
        srcfile.write("eval `scramv1 runtime -sh`\n\n")
        srcfile.write(cmd + "\n\n")
        ## add some protection against core dumps which might write large log files
        ## it used to happen with lxbatch, not sure it is needed with condor
        srcfile.write('''
if ls {pp}/core.* 1> /dev/null 2>&1; then
    echo "Found core.* files in {pp}/."
    echo "Removing them!"
    rm {pp}/core.*
else
    echo "No core.* files found in {pp}/"
fi\n\n
'''.format(pp=plotterPath))
        srcfile.close()  
        # close this .sh file

        ## add this bash script inside the condir submit one
        of.write("arguments = {src}\n".format(src=os.path.abspath(srcfileName)))
        of.write("queue 1\n\n")

        
    ####
    of.close()
    print ""
    print "Condor file saved in " + outf
    if not options.runlocal:   
        if options.dryrun:
            print "Option -d active, not running commands"
            print "condor_submit " + outf
        else:
            print "Executing it:"
            print "condor_submit " + outf
            print "-"*30
            os.system("condor_submit " + outf)
            print "-"*30
            print "Jobs submitted, good luck!"
    print ""
        
