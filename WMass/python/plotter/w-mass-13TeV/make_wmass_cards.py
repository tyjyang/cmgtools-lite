#!/usr/bin/env python
from shutil import copyfile
import re, sys, os, os.path, subprocess, json, ROOT
import numpy as np

def printReturnValueInScript(fileObj, scriptName=None):
    fileObj.write('retval=$?\n') # must be evaluated immediately after the python command in the shell script
    if scriptName != None:
        fileObj.write('echo "Script name = %s"\n' % scriptName)
    fileObj.write('echo "return value = ${retval}"\n')
    fileObj.write('if [[ $retval != 0 ]]; then\n')
    fileObj.write('    exit $retval\n')
    fileObj.write('fi\n')
    fileObj.write('\n')


def getShFile(jobdir, name):
    tmp_srcfile_name = jobdir+'/job_{i}.sh'.format(i=name)
    tmp_srcfile = open(tmp_srcfile_name, 'w')
    tmp_srcfile.write("#! /bin/sh\n")
    tmp_srcfile.write("ulimit -c 0 -S\n")
    tmp_srcfile.write("ulimit -c 0 -H\n")
    tmp_srcfile.write("cd {cmssw};\neval $(scramv1 runtime -sh);\ncd {d};\n".format( d= os.getcwd(), cmssw = os.environ['CMSSW_BASE']))
    return tmp_srcfile_name, tmp_srcfile

def wptBinsScales(i):
    wptbins = [0.0, 1.971, 2.949, 3.838, 4.733, 5.674, 6.684, 7.781, 8.979, 10.303, 11.777, 13.435, 15.332, 17.525, 20.115, 23.245, 27.173, 32.414, 40.151, 53.858, 13000.0]
    if len(wptbins)<2*i:
        print 'you are asking too much from the wpt binning for decorrelation of scales'
    ptlo = wptbins[2*(i-1)]
    pthi = wptbins[2*i]
    return [ptlo, pthi]
    

def getCondorTime(qstr):
    retval = ''
    if   qstr == '1nh':
        retval = 3600
    elif qstr == '8nh':
        retval = 28800
    elif qstr == '1nd':
        retval = 24*3600
    elif qstr == '2nd':
        retval = 48*3600
    elif qstr == '1nw':
        retval = 7*24*3600
    else:
        retval = int(qstr)*60*60
    return int(retval)

def printSysts(systs=[],process="Wmunu"):
    if len(systs):
        print '-'*30
        print "Have these systematics on {p}:".format(p=process)
        print systs
        print '-'*30


NVPTBINS=10 # usually it would be 10, use less for tests, e.g. 2
NPDFSYSTS=100 # Hessian variations (from 1 to 100), use less for tests, e.g. 2 
nominals=[] # array containing the nominal for single process for which we have dedicated corrections (not included in bkg_and_data)
#
# dictionaries associating an array containing systs variations to each relevant process
pdfsysts={} # PDFs variations
qcdsysts={} # QCD scale variations (binned in W boson pt)
inclqcdsysts={} # inclusive QCD scale variations (e.g. for Z)
etaeffsysts={} # uncorrelated efficiency systematics vs eta
fsrsysts={} # (only) the PHOTOS/PYTHIA fsr reweighitng
kamucasysts={} # kalman filter muon momentum scale systematics

#coefficients = ['ac'] + ['a'+str(i) for i in range(8)]
coefficients = [''] # simplifies life for now, but do not use empty array, or it will inhibit some loops
data_eras = ["B", "C", "D", "E", "F", "F_postVFP", "G", "H"]

def getMcaIncl(mcafile,incl_mca='incl_sig'):
    incl_file=''
    mcaf = open(mcafile,'r')
    for l in mcaf.readlines():
        if re.match("\s*#.*", l): continue
        tokens = [t.strip() for t in l.split(':')]
        if len(tokens)<2: continue
        if tokens[0]==incl_mca and "+" in tokens[1]:
            options_str = [t.strip() for t in (l.split(';')[1]).split(',')]
            for o in options_str:
                if "IncludeMca" in o: 
                    incl_file = o.split('=')[1]
            break
    return incl_file

def writeNominalSingleMCA(mcafile,odir,incl_mca='incl_zmumu'):
    incl_file=getMcaIncl(mcafile,incl_mca)
    process = incl_mca.split('_')[1]
    postfix = "_nominal"
    mcafile_nominal = open("{odir}/mca_{p}{pfx}.txt".format(odir=odir,p=process,pfx=postfix), "w")
    mcafile_nominal.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="1." \n')
    nominals.append(process)
    print "written nominal mca relative to ",incl_mca

def writePdfSystsToMCA(mcafile,odir,vec_weight="pdfWeightNNPDF",syst="pdf",incl_mca='incl_wmunu',append=False):
    open("%s/systEnv-dummy.txt" % odir, 'a').close()
    incl_file=getMcaIncl(mcafile,incl_mca)
    if len(incl_file)==0: 
        print "Warning! '%s' include directive not found. Not adding pdf systematics samples to MCA file" % incl_mca
        return
    if append:
        filename = "%s/mca_systs.txt" % odir
        if not os.path.exists(filename): os.system('cp {mca_orig} {mca_syst}'.format(mca_orig=mcafile,mca_syst=filename))
    for i in range(1,NPDFSYSTS+1):
        process = incl_mca.split('_')[1]
        postfix = "_{syst}{idx}".format(syst=syst,idx=i)
        mcafile_syst = open(filename, 'a') if append else open("%s/mca_%s%s.txt" % (odir,process,postfix), "w")
        mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+vec_weight+str(i)+'", PostFix="'+postfix+'" \n')
        if process not in pdfsysts:
            pdfsysts[process] = []
        pdfsysts[process].append(postfix)
    print "written ",syst," systematics relative to ",incl_mca


## this function here is pretty important to do all the theory systematics
## ---
def writeQCDScaleSystsToMCA(mcafile,odir,syst="qcd",incl_mca='incl_wmunu',scales=[],append=False,ptBinned=True,overrideDecorrelation=False):
    open("%s/systEnv-dummy.txt" % odir, 'a').close()
    incl_file=getMcaIncl(mcafile,incl_mca)

    if len(incl_file)==0: 
        print "Warning! '%s' include directive not found. Not adding QCD scale systematics!"
        return

    if append:
        filename = "%s/mca_systs.txt" % odir
        if not os.path.exists(filename): os.system('cp {mca_orig} {mca_syst}'.format(mca_orig=mcafile,mca_syst=filename))    

    scaleVarToNtupleVar = {
        "muRUp"    : "scaleWeightMuR2MuF1",
        "muRDn"    : "scaleWeightMuR05MuF1",
        "muFUp"    : "scaleWeightMuR1MuF2",
        "muFDn"    : "scaleWeightMuR1MuF05",
        "muRmuFUp" : "scaleWeightMuR2MuF2",
        "muRmuFDn" : "scaleWeightMuR05MuF05"
    }

    process = incl_mca.split('_')[1]
    ## make the mcas for the scales and mW and stuff
    for scale in scales:

        ## mW now doesn't just have Up and Down. make many
        if scale == "mW":

            #############################
            # FIXME: STILL TO BE UPDATED
            #############################

            ## loop on all the masses we want. in index of 5 MeV variations each
            masses  = ['mass_{m}'.format(m = j).replace('-','m') for j in range(-20,0)]
            masses += ['mass_0']
            masses += ['mass_p{m}'.format(m = j).replace('-','m') for j in range(1,21)]

            for mass in masses:
                postfix = "_{syst}{mval}".format(syst=scale,mval=mass)
                mcafile_syst = open(filename, 'a') if append else open("%s/mca_%s%s.txt" % (odir,process,postfix), "w")

                ## central mass is 80419 MeV, the closest we have to that is 80420. will scale +- 50 MeV, i.e. 80470 for Up and 80370 for Dn
                fstring = str(mass)
                mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+fstring+'", PostFix="'+postfix+'" \n')
                if process not in qcdsysts:
                    qcdsysts[process] = []
                qcdsysts[process].append(postfix)

        else:
            ## all the others as usual with an Up and Down variation
            for idir in ['Up','Dn']:
                postfix = "_{syst}{idir}".format(syst=scale,idir=idir)
                mcafile_syst = open(filename, 'a') if append else open("%s/mca_%s%s.txt" % (odir,process,postfix), "w")
                if 'muR' in scale or 'muF' in scale:
    
                    scaleVar = scaleVarToNtupleVar["{sc}{idir}".format(sc=scale,idir=idir)]
                    ## for the signal keep as is for the QCD scales
                    if ptBinned:
                        for ipt in range(1,1+NVPTBINS): ## start from 1 to 10
                            for coeff in coefficients if options.decorrelateSignalScales and not overrideDecorrelation else ['']:
                                for pm in ['plus', 'minus']:
                                    ## have to redo the postfix for these
                                    postfix = "_{coeff}{syst}{ipt}{ch}{idir}".format(syst=scale,idir=idir,ipt=ipt,coeff=coeff,ch=pm)
                                    mcafile_syst = open(filename, 'a') if append else open("%s/mca_%s%s.txt" % (odir,process,postfix), "w")
                                    ptcut = wptBinsScales(ipt)
                                    wgtstr = 'TMath::Power({sc}\,(Vpt_preFSR>={ptlo}&&Vpt_preFSR<{pthi}))'.format(sc=scaleVar,ptlo=ptcut[0],pthi=ptcut[1])
                                    mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+wgtstr+'", PostFix="'+postfix+'" \n')
                                    if process not in qcdsysts:
                                        qcdsysts[process] = []
                                    qcdsysts[process].append(postfix)
                    else:
                        ## for the Z only do the unnumbered ones                        
                        mcafile_syst = open(filename, 'a') if append else open("%s/mca_%s%s.txt" % (odir,process,postfix), "w")
                        wgtstr = '{sc}'.format(sc=scaleVar)
                        mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+wgtstr+'", PostFix="'+postfix+'" \n')
                        if process not in inclqcdsysts:
                            inclqcdsysts[process] = []
                        inclqcdsysts[process].append(postfix)
    
                else: ## alphaS is left here
                    ## don't do mW here again, although it should not be found here
                    if "mW" in scale:
                        continue
                    ## ---
                    wgtstr = "pdfWeightNNPDF101" if idir == "Up" else "pdfWeightNNPDF102"
                    mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+wgtstr+'", PostFix="'+postfix+'" \n')
                    if process not in inclqcdsysts:
                        inclqcdsysts[process] = []
                    inclqcdsysts[process].append(postfix)
                    # was like that, but why would alpha go in both arrays? keep commented
                    #if process not in qcdsysts:
                    #    qcdsysts[process] = []
                    #qcdsysts.append(postfix)
    
    print "written ",syst," systematics relative to ",incl_mca

def writeEfficiencyStatErrorSystsToMCA(mcafile,odir,channel,syst="EffStat",incl_mca='incl_wmunu',append=False):
    open("%s/systEnv-dummy.txt" % odir, 'a').close()
    incl_file=getMcaIncl(mcafile,incl_mca)
    if len(incl_file)==0: 
        print "Warning! '%s' include directive not found. Not adding pdf systematics samples to MCA file" % incl_mca
        return
    if append:
        filename = "%s/mca_systs.txt" % odir
        if not os.path.exists(filename): os.system('cp {mca_orig} {mca_syst}'.format(mca_orig=mcafile,mca_syst=filename))
    etalo = -2.5 if channel=='el' else -2.4001
    deta = 0.1; nbins = int(2*abs(etalo)/deta)
    process = incl_mca.split('_')[1]
    for i in range(0,nbins):
        etamin=etalo + i*deta; etamax=etalo + (i+1)*deta;
        # 3 parameters used in the Erf fit vs pt per eta bin
        for ipar in xrange(3):
            postfix = "_ErfPar{ipar}{syst}{idx}".format(syst=syst,ipar=ipar,idx=i+1)
            mcafile_syst = open(filename, 'a') if append else open("%s/mca_%s%s.txt" % (odir,process,postfix), "w")
            weightFcn = 'effSystEtaBins({ipar}\,Muon_pdgId[0]\,Muon_eta[0]\,Muon_pt[0]\,{emin:.1f}\,{emax:.1f})'.format(ipar=ipar,emin=etamin,emax=etamax)
            mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+weightFcn+'", PostFix="'+postfix+'" \n')
            etaeffsysts.append(postfix)
    print "written ",syst," systematics relative to ",incl_mca

# FIXEM: this is probably obsolete, keep for reference for now
def writeFSRSystsToMCA(mcafile,odir,syst="fsr",incl_mca='incl_wmunu',append=False):
    open("%s/systEnv-dummy.txt" % odir, 'a').close()
    incl_file=getMcaIncl(mcafile,incl_mca)
    process = incl_mca.split('_')[1]
    postfix = "_{syst}".format(syst=syst)
    mcafile_syst = open("%s/mca_%s%s.txt" % (odir,process,postfix), "w")
    # this is going to be fully changed, probably
    weightFcn = 'fsrPhotosWeightSimple(GenLepDressed_pdgId[0]\,GenLepDressed_pt[0]\,GenLepBare_pt[0])'
    mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+weightFcn+'", PostFix="'+postfix+'" \n')
    if process not in fsrsysts:
        fsrsysts[process] = []
    fsrsysts[process].append(postfix)
    print "written ",syst," systematics relative to ",incl_mca

# not used for now, let's comment it
# def writePdfSystsToSystFile(filename,sample="W.*",syst="CMS_W_pdf"):
#     SYSTFILEALL=('.').join(filename.split('.')[:-1])+"-all.txt"
#     copyfile(filename,SYSTFILEALL)
#     systfile=open(SYSTFILEALL,"a")
#     for i in range(NPDFSYSTS/2):
#         systfile.write(syst+str(i+1)+"  : "+sample+" : .* : pdf"+str(i+1)+" : templates\n")
#     print "written pdf syst configuration to ",SYSTFILEALL
#     return SYSTFILEALL

# obsolete
def addKaMuCaSysts():
    for idir in ['Up','Dn']:
        for i in range(133):
            kamucasysts.append('kalPtErr{i}{idir}'.format(i=i,idir=idir))
        kamucasysts.append('kalPtClosureErr{idir}'.format(idir=idir))
    print "Added muon momentum scale systematics"

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] mc.txt cuts.txt var bins systs.txt outdir ")
    parser.add_option("-q", "--queue",    dest="queue",     type="string", default=None, help="Run jobs on lxbatch instead of locally");
    parser.add_option("-n", "--job-name",    dest="jobName",     type="string", default=None, help="Select name for condor jobs (by default they are identified by the cluster ID)");
    parser.add_option("-l", "--lumi",    dest="integratedLuminosity",     type="float", default=35.9, help="Integrated luminosity");
    parser.add_option("--dry-run", dest="dryRun",    action="store_true", default=False, help="Do not run the job, only print the command");
    parser.add_option("--decorrelateSignalScales", action="store_true", default=False, help="Treat qcd scales uncorrelated for left/right/long");
    parser.add_option("-s", "--signal-cards",  dest="signalCards",  action="store_true", default=False, help="Make the signal part of the datacards");
    parser.add_option("-b", "--bkgdata-cards", dest="bkgdataCards", action="store_true", default=False, help="Make the background and data part of the datacards");
    parser.add_option("-W", "--weight", dest="weightExpr", default="-W 1", help="Event weight expression (default 1)");
    parser.add_option("-P", "--path", dest="path", type="string",default=None, help="Path to directory with input trees and pickle files");
    parser.add_option("-C", "--channel", dest="channel", type="string", default='mu', help="Channel. either 'el' or 'mu'");
    parser.add_option("--not-unroll2D", dest="notUnroll2D", action="store_true", default=False, help="Do not unroll the TH2Ds in TH1Ds needed for combine (to make 2D plots)");
    parser.add_option("--pdf-syst", dest="addPdfSyst", action="store_true", default=False, help="Add PDF systematics to the signal (need incl_xxx directive in the MCA file)");
    parser.add_option("--qcd-syst", dest="addQCDSyst", action="store_true", default=False, help="Add QCD scale systematics to the signal (need incl_xxx directive in the MCA file)");
    parser.add_option("--qed-syst", dest="addQEDSyst", action="store_true", default=False, help="Add QED scale systematics to the signal (need incl_xxx directive in the MCA file)");
    parser.add_option("--kamuca-syst", dest="addKaMuCaSyst", action="store_true", default=False, help="Add Kalman Filter muon momentum scale correction systematics");
    parser.add_option('-g', "--group-jobs", dest="groupJobs", type=int, default=20, help="group signal jobs so that one job runs multiple makeHistogramsWMass commands");
    parser.add_option('--vpt-weight', dest='procsToPtReweight', action="append", default=[], help="processes to be reweighted according the measured/predicted DY pt. Default is none (possible W,Z).");
    parser.add_option('--wlike', dest='wlike', action="store_true", default=False, help="Make cards for the wlike analysis. Default is wmass");    
    parser.add_option('--add-option', dest="addOptions", type="string", default=None, help="add these options to the option string when running the histograms");
    parser.add_option('--auto-resub', dest='automaticResubmission', action="store_true", default=False, help="Use condor features for automatic job resubmission in case of failures");
    parser.add_option("--n-resub", dest="nResubmissions", type=int, default=2, help="Number of automatic resubmission to be attempted in case of failure when running a job (needs --auto-resub)");
    (options, args) = parser.parse_args()
    
    if len(sys.argv) < 6:
        parser.print_usage()
        quit()
        
    print ""
    print "="*30
    print ""
    FASTTEST=''
    #FASTTEST='--max-entries 1000 '
    T=options.path
    print "used trees from: ",T
    J=2
    MCA = args[0]
    CUTFILE = args[1]
    fitvar = args[2]
    binning = args[3]
    SYSTFILE = args[4]
    
    luminosity = options.integratedLuminosity
    
    if not os.path.exists("cards/"):
        os.makedirs("cards/")
    outdir="cards/"+args[5]
    
    if not os.path.exists(outdir): 
        os.mkdir(outdir)
    if options.queue:
        foldersToCreate = ["/jobs", "/mca", "/logs", "/outs", "/errs"]
        for f in foldersToCreate:
            if not os.path.exists(outdir+f): 
                os.mkdir(outdir+f)
    
    # copy some cfg for bookkeeping
    os.system("cp %s %s" % (CUTFILE, outdir))
    os.system("cp %s %s" % (MCA, outdir))
    
    ## save template binning (eta on X, pt on y axis)
    ptEta_binfile = open(outdir+'/binningPtEta.txt','w')
    ptEta_binfile.write("#Template binning: eta-pt on x-y axis\n")
    ptEta_binfile.write("reco: "+binning+"\n")
    ptEta_binfile.write('\n')
    ptEta_binfile.close()
    
    if options.addPdfSyst:
        # write the additional systematic samples in the MCA file
        writePdfSystsToMCA(MCA,outdir+"/mca",incl_mca='incl_wmunu') # on W + jets 
        writePdfSystsToMCA(MCA,outdir+"/mca",incl_mca='incl_wtaunu') # on WTau
        writePdfSystsToMCA(MCA,outdir+"/mca",incl_mca='incl_zmumu') # on DY + jets
        writePdfSystsToMCA(MCA,outdir+"/mca",incl_mca='incl_ztautau') # overkilling?
        # write the complete systematics file (this was needed when trying to run all systs in one job)
        # SYSTFILEALL = writePdfSystsToSystFile(SYSTFILE)
    if options.addQCDSyst:
        scales = ['muR','muF',"muRmuF", "alphaS"]
        if options.wlike:
            wmasses = []
            zmasses = [] # ["mW"]
            ptBinnedScalesForW = False
        else:
            wmasses = [] #["mW"]
            zmasses = []
            ptBinnedScalesForW = True
        # overrideDecorrelation is used to avoid decorrelating systs by angular coefficients
        writeQCDScaleSystsToMCA(MCA,outdir+"/mca",scales=scales+wmasses,incl_mca='incl_wmunu',
                                ptBinned=ptBinnedScalesForW, overrideDecorrelation=False)
        writeQCDScaleSystsToMCA(MCA,outdir+"/mca",scales=scales        ,incl_mca='incl_wtaunu',
                                ptBinned=ptBinnedScalesForW, overrideDecorrelation=True)
        writeQCDScaleSystsToMCA(MCA,outdir+"/mca",scales=scales+zmasses,incl_mca='incl_zmumu',
                                ptBinned=options.wlike, overrideDecorrelation=True)
        writeQCDScaleSystsToMCA(MCA,outdir+"/mca",scales=scales        ,incl_mca='incl_ztautau',
                                ptBinned=options.wlike, overrideDecorrelation=True)
    if options.addQEDSyst:
        if options.wlike:
            writeFSRSystsToMCA(MCA,outdir+"/mca",incl_mca='incl_zmumu') # DY + jets
        else:
            writeFSRSystsToMCA(MCA,outdir+"/mca") # on W + jets
    
    # apparently needed to make mca for nominal samples to run faster, to be checked
    procs_nomi_mca = ['wtaunu','ztautau'] + (['wmunu'] if options.wlike else ['zmumu'])
    for proc in procs_nomi_mca:
        if proc not in nominals: 
            writeNominalSingleMCA(MCA,outdir+"/mca",incl_mca='incl_{p}'.format(p=proc))
    
    if options.addKaMuCaSyst:
        addKaMuCaSysts()

    # not needed if we fit eta/pt. Will be needed if we fit another variable correlated with eta/pt
    # writeEfficiencyStatErrorSystsToMCA(MCA,outdir+"/mca",options.channel)
    
    ARGS=" ".join([MCA,CUTFILE,"'"+fitvar+"' "+"'"+binning+"'",SYSTFILE])
    BASECONFIG=os.path.dirname(MCA)
    ## use rel paths if options.queue:
    ## use rel paths     ARGS = ARGS.replace(BASECONFIG,os.getcwd()+"/"+BASECONFIG)
    OPTIONS=" -P "+T+" --s2v -j "+str(J)+" -l "+str(luminosity)+" -f --obj Events "+FASTTEST
    if not options.notUnroll2D:
        OPTIONS+=" --2d-binning-function unroll2Dto1D "
    if options.addOptions:
        OPTIONS += options.addOptions
    #print "-"*30
    #print OPTIONS
    #print "-"*30

    if options.wlike:
        POSCUT=" -A alwaystrue positive 'event%2 != 0' "
        NEGCUT=" -A alwaystrue negative 'event%2 == 0' "    
        SIGPROC='Zmumu'
        SIGSUFFIX='zmumu'
    else:
        POSCUT=" -A alwaystrue positive 'Muon_charge[0]>0' "
        NEGCUT=" -A alwaystrue negative 'Muon_charge[0]<0' "        
        SIGPROC='Wmunu'
        SIGSUFFIX='wmunu'
    antiSIGPROC = 'Zmumu' if SIGPROC=='Wmunu' else 'Wmunu'

    #print 'these are angular coefficients ', coefficients
    print ""
    print ""

    fullJobList = set()
    if options.signalCards:
        print "MAKING SIGNAL PART!"
        sigsyst = ['']
        if SIGSUFFIX in pdfsysts:
            sigsyst += pdfsysts[SIGSUFFIX]
        if SIGSUFFIX in qcdsysts:
            sigsyst += qcdsysts[SIGSUFFIX]
        if SIGSUFFIX in inclqcdsysts:
            sigsyst += inclqcdsysts[SIGSUFFIX]
        if SIGSUFFIX in etaeffsysts:
            sigsyst += etaeffsysts[SIGSUFFIX]
        if SIGSUFFIX in fsrsysts:
            sigsyst += fsrsysts[SIGSUFFIX]
        sigsyst += kamucasysts
        printSysts(systs=sigsyst,process=SIGPROC)

        ## loop on both charges
        for charge in ['plus', 'minus']:
            antich = 'plus' if charge == 'minus' else 'minus'
            ycut = POSCUT if charge=='plus' else NEGCUT

            for coeff in coefficients:
                print "Making histograms for signal process {coeff} with charge {ch} ".format(coeff=coeff,ch=charge)
                ## get the list of all other coefficients
                anticoeff = [proc for proc in coefficients if not proc==coeff]
                ## loop on all the variations
                for ivar,var in enumerate(sigsyst):
                    ## exclude the anti charge and the anti coefficients
                    if antich in var: 
                        #print "skipping because ",antich," is in ",var
                        continue
                    if len(anticoeff) and any(i in var for i in anticoeff) and options.decorrelateSignalScales: 
                        continue ## this might work. but who knows...
                    SYST_OPTION = ''
                    ## 0 is nominal, all the theory variations are non-zero
                    if ivar==0: 
                        IARGS = ARGS
                    else: 
                        ## keep kamuca part here for now, but will need to be updated
                        if 'kalPt' in var:
                            IARGS = re.sub(r'LepGood(\d+)_pt',r'LepGood\g<1>_kalPt{systvar}'.format(systvar=var.replace('kalPt','')),ARGS)
                            SYST_OPTION = ' --replace-var LepGood1_pt LepGood1_{systvar} --replace-var LepGood2_pt LepGood2_{systvar} '.format(systvar=var)
                        else:
                            IARGS = ARGS.replace(MCA,"{outdir}/mca/mca_{p}{syst}.txt".format(outdir=outdir,p=SIGSUFFIX,syst=var))
                            IARGS = IARGS.replace(SYSTFILE,"{outdir}/mca/systEnv-dummy.txt".format(outdir=outdir))
                        # let's not print too many things
                        # print "Running the systematic: ",var
                    ## ---

                    ## now go and make the submit commands
                    xpsel=' --xp "{sig}_{antich}.*,{antisig}.*,Top,DiBosons,Ztautau.*,Wtaunu.*,data.*"'.format(sig=SIGPROC,antisig=antiSIGPROC,antich=antich,ch=charge)                        
                    if len(anticoeff):
                        excl_anticoeff = ','.join(SIGPROC+"_"+charge+'_'+ac+'.*' for ac in anticoeff)
                        xpsel += " --xp {ahel}".format(ahel=excl_anticoeff)
                    xpsel += " --asimov "

                    ## make specific names and such
                    syst = '' if ivar==0 else var
                    ## for kamuca corrections, there is not any _sigsuffix, so add an "_" to distinguish better the files
                    if 'kalPt' in var: syst = '_{sfx}_{var}'.format(sfx=SIGSUFFIX,var=var)
                    dcname = "{sig}_{charge}{coeff}{syst}".format(sig=SIGPROC, charge=charge, coeff="" if coeff=='' else ('_'+coeff),syst=syst)

                    ## for Wlike, must transform process name for signal (Zmumu) into Zmumu_plus or Zmumu_minus.
                    ## one could directly change the name of the histogram inside the output root file, after it has been created, but maybe better to do it now.
                    ## This can be achieved either by using the IncludeMca statement and adding the '_plus' or '_minus' postfix to the process name, or by using the --pg option to change the signal process name
                    ## for now we will try the second option, where the new process name is basically dcname (because one has to include the general case with systematics)
                    if options.wlike:
                        #xpsel += " -p Zmumu_{ch} --pg 'Zmumu_{ch} := Zmumu' ".format(ch=charge)
                        syst_postfix = "" if ivar==0 else syst
                        xpsel += " -p {dc} --pg '{dc} := Zmumu{s}' ".format(dc=dcname,s=syst_postfix)


                    ## reweight the boson-pT here
                    zptWeight = 'dyptWeight(Vpt_preFSR,{isZ})'.format(isZ=0 if SIGPROC=="Wmunu" else 1)

                    ## construct the full weight to be applied
                    fullWeight = options.weightExpr+'*'+zptWeight if SIGPROC in options.procsToPtReweight else options.weightExpr
                    BIN_OPTS = OPTIONS + " -W '" + fullWeight+ "'" + " -o "+dcname+" --od "+outdir + xpsel + ycut + SYST_OPTION

                    ## do not do by pdf weight pdfmatch = re.search('pdf(\d+)',var)
                    ## do not do by pdf weight if pdfmatch:
                    ## do not do by pdf weight     BIN_OPTS += " -W 'pdfRatioHel(abs(prefsrw_y),prefsrw_pt,prefsrw_costcs,evt,{pol},{ipdf})' ".format(pol=helmap[helicity],ipdf=pdfmatch.group(1))
                    ## ---

                    ## now we accumulate all the commands to be run for all processes
                    if options.queue:
                        mkShCardsCmd = "python makeHistogramsWMass.py {args} \n".format(dir = os.getcwd(), args = IARGS+" "+BIN_OPTS)
                        fullJobList.add(mkShCardsCmd)
                    ## ---
    
    ## this is for all the background cards, excluding the Z(W) for Wmass(Wlike), Wtaunu, and Ztautau
    if options.bkgdataCards:
        print "MAKING BKG and DATA PART:\n"
        # will split data from other backgrounds, and also split data in multiple eras for faster jobs
        ## again loop on both charges
        for charge in ['plus','minus']:
            print "Making histograms for other background processes with charge ", charge
            chargecut = POSCUT if charge=='plus' else NEGCUT
            ## remove W or Z signal processes and others that aren't needed            
            xpsel=' --xp "{sig}.*,{antisig}.*,Wtaunu.*,Ztautau.*,data.*" '.format(sig=SIGPROC,antisig=antiSIGPROC) 
            ## now make the names of the cards etc
            dcname = "otherBkgHisto_{charge}".format(charge=charge)
            BIN_OPTS=OPTIONS + " -W '" + options.weightExpr + "'" + " -o "+dcname+" --od "+outdir + xpsel + chargecut
            ## accumulate the commands to be run
            if options.queue:
                mkShCardsCmd = "python makeHistogramsWMass.py {args} \n".format(dir = os.getcwd(), args = ARGS+" "+BIN_OPTS)
                fullJobList.add(mkShCardsCmd)
            ##
            ## ---
            ## now data for different eras
            ## see /eos/cms/store/cmst3/group/wmass/w-mass-13TeV/postNANO/dec2020/DATA_NanoV8/
            for era in data_eras:
                print "Making histograms for data {e} with charge {ch}".format(e=era,ch=charge)
                excluded_eras = ["data_{e}".format(e=e) for e in data_eras if e != era]
                psel=" -p 'data,data_fakes' --pg 'data := data_{era}' --pg 'data_fakes := data_{era}_fakes' --xp '{x}'".format(era=era, x=",".join(excluded_eras))
                ## now make the names of the cards etc
                dcname = "dataHisto_{era}_{charge}".format(era=era,charge=charge)
                BIN_OPTS=OPTIONS + " -W '" + options.weightExpr + "'" + " -o "+dcname+" --od "+outdir + psel + chargecut
                ## accumulate the commands to be run
                if options.queue:
                    mkShCardsCmd = "python makeHistogramsWMass.py {args} \n".format(args = ARGS+" "+BIN_OPTS)
                    fullJobList.add(mkShCardsCmd)
                ## ---


    
    ## this is for background cards which include the Z(W) samples for Wmass (Wlike)
        antiSIGSUFFIX = 'zmumu' if SIGPROC=='Wmunu' else 'wmunu'
        antisigsyst = [''] # nominal
        if antiSIGSUFFIX in pdfsysts: 
            antisigsyst += pdfsysts[antiSIGSUFFIX]        
        if antiSIGSUFFIX in qcdsysts: 
            antisigsyst += qcdsysts[antiSIGSUFFIX]
        if antiSIGSUFFIX in inclqcdsysts:
            antisigsyst += inclqcdsysts[antiSIGSUFFIX]
        printSysts(systs=antisigsyst,process=antiSIGPROC)
        ## loop on both charges
        for charge in ['plus','minus']:
            print "Making histograms for",antiSIGPROC," process with charge ", charge
            chcut = POSCUT if charge=='plus' else NEGCUT
            antich = 'plus' if charge == 'minus' else 'minus'
            ## loop on the Z(W) related systematics for Wmass (Wlike)
            for ivar,var in enumerate(antisigsyst):
    
                ## nominal first, then the variations
                if ivar==0: 
                    # the following would be faster with DY-only, but it misses processes for possible systematics, like the lines for the Z_lepeff systs that we had in the past
                    IARGS = ARGS.replace(MCA,"{outdir}/mca/mca_{p}_nominal.txt".format(outdir=outdir,p=antiSIGSUFFIX))
                else: 
                    IARGS = ARGS.replace(MCA,"{outdir}/mca/mca_{p}{syst}.txt".format(outdir=outdir,p=antiSIGSUFFIX,syst=var))
                    IARGS = IARGS.replace(SYSTFILE,"{outdir}/mca/systEnv-dummy.txt".format(outdir=outdir))
                    # print "Running Zmumu with systematic: ",var

                ## ---
                # exclude everything that does not start with antiSIGPROC 
                # for Wmass the antisignal is the Zmumu, and the charge is not in process name in the mca file
                # for Wlike, the antisignal is Wmunu and has to pick the samples with proper charge
                psel=' -p "^{antisig}_{ch}.*" --asimov '.format(antisig=antiSIGPROC,ch="" if SIGPROC=='Wmunu' else charge)
                syst = '' if ivar==0 else var

                ## make names for the files and outputs
                dcname = "{antisig}_{charge}{syst}".format(antisig=antiSIGPROC, charge=charge,syst=syst)

                ## construct the final weight. reweight also the DY to whatever new pT spectrum we want
                zptWeight = 'dyptWeight(Vpt_preFSR,{isZ})'.format(isZ=1 if SIGPROC=='Wmunu' else 0)
                fullWeight = options.weightExpr+'*'+zptWeight if antiSIGPROC in options.procsToPtReweight else options.weightExpr
                BIN_OPTS=OPTIONS + " -W '" + fullWeight + "'" + " -o "+dcname+" --od "+outdir + psel + chcut
                ## ---

                ## accumulate the commands to be run
                if options.queue:
                    mkShCardsCmd = "python makeHistogramsWMass.py {args} \n".format(dir = os.getcwd(), args = IARGS+" "+BIN_OPTS)
                    fullJobList.add(mkShCardsCmd)
                ## ---
    
    ## repetition for Wtau, but better to keep Z/Tau cards separated (faster jobs)
        wtausyst = [''] # nominal
        if "wtaunu" in pdfsysts:
            wtausyst += pdfsysts["wtaunu"]    
        if "wtaunu" in qcdsysts:
            wtausyst += qcdsysts["wtaunu"]
        if "wtaunu" in inclqcdsysts:
            wtausyst += inclqcdsysts["wtaunu"]        
        printSysts(systs=wtausyst,process="Wtaunu")
        ## loop on systmatic variations
        for charge in ['plus','minus']:
            print "Making histograms for Wtaunu process with charge ", charge
            ## loop on both charges
            chcut = POSCUT if charge=='plus' else NEGCUT
            antich = 'plus' if charge == 'minus' else 'minus'
            for ivar,var in enumerate(wtausyst):
        
                ## 0 is nominal again
                if ivar==0: 
                    IARGS = ARGS.replace(MCA,"{outdir}/mca/mca_wtaunu_nominal.txt".format(outdir=outdir))
                else: 
                    IARGS = ARGS.replace(MCA,"{outdir}/mca/mca_wtaunu{syst}.txt".format(outdir=outdir,syst=var))
                    IARGS = IARGS.replace(SYSTFILE,"{outdir}/mca/systEnv-dummy.txt".format(outdir=outdir))
                    # print "Running the Wtaunu with systematic: ",var

                # exclude everything that does not start with Wtaunu with the proper charge
                psel=' -p "^Wtaunu_{ch}.*" --asimov '.format(antisig=antiSIGPROC,ch=charge)
                syst = '' if ivar==0 else var

                ## make names for files etc again
                dcname = "Wtaunu_{charge}{syst}".format(charge=charge,syst=syst)
                ## construct full reweighting weight. boson-pT reweighting here again
                zptWeight = 'dyptWeight(Vpt_preFSR,0)'
                fullWeight = options.weightExpr+'*'+zptWeight if 'Wtaunu' in options.procsToPtReweight else options.weightExpr
                BIN_OPTS=OPTIONS + " -W '" + fullWeight + "'" + " -o "+dcname+" --od "+outdir + psel + chcut
                ## ---

                ## accumulate the commands to be run
                if options.queue:
                    mkShCardsCmd = "python makeHistogramsWMass.py {args} \n".format(dir = os.getcwd(), args = IARGS+" "+BIN_OPTS)
                    fullJobList.add(mkShCardsCmd)
                ## ---

#########
########

    ## repetition for Ztautau
        ztausyst = [''] # nominal
        if "ztautau" in pdfsysts:
            ztausyst += pdfsysts["ztautau"]
        if "ztautau" in qcdsysts:
            ztausyst += qcdsysts["ztautau"]
        if "ztautau" in inclqcdsysts:
            ztausyst += inclqcdsysts["ztautau"]
        printSysts(systs=ztausyst,process="Ztautau")
        ## loop on systmatic variations
        for charge in ['plus','minus']:
            print "Making histograms for Ztautau process with charge ", charge
            ## loop on both charges
            chcut = POSCUT if charge=='plus' else NEGCUT
            antich = 'plus' if charge == 'minus' else 'minus'
            for ivar,var in enumerate(ztausyst):
        
                ## 0 is nominal again
                if ivar==0: 
                    IARGS = ARGS.replace(MCA,"{outdir}/mca/mca_ztautau_nominal.txt".format(outdir=outdir))
                else: 
                    IARGS = ARGS.replace(MCA,"{outdir}/mca/mca_ztautau{syst}.txt".format(outdir=outdir,syst=var))
                    IARGS = IARGS.replace(SYSTFILE,"{outdir}/mca/systEnv-dummy.txt".format(outdir=outdir))
                    # print "Running the Ztautau with systematic: ",var

                # exclude everything that does not start with Ztautau with the proper charge
                psel=' -p "^Ztautau_.*" --asimov '
                syst = '' if ivar==0 else var

                ## make names for files etc again
                dcname = "Ztautau_{charge}{syst}".format(charge=charge,syst=syst)
                ## construct full reweighting weight. boson-pT reweighting here again
                zptWeight = 'dyptWeight(Vpt_preFSR,1)'
                fullWeight = options.weightExpr+'*'+zptWeight if 'Ztautau' in options.procsToPtReweight else options.weightExpr
                BIN_OPTS=OPTIONS + " -W '" + fullWeight + "'" + " -o "+dcname+" --od "+outdir + psel + chcut
                ## ---

                ## accumulate the commands to be run
                if options.queue:
                    mkShCardsCmd = "python makeHistogramsWMass.py {args} \n".format(dir = os.getcwd(), args = IARGS+" "+BIN_OPTS)
                    fullJobList.add(mkShCardsCmd)
                ## ---


#######
#######
        
    ## now loop on all the accumulated jobs that we have from above
    if len(fullJobList):
        reslist = list(fullJobList)
    
        ## maximum 6000 output cards and files per directory. to avoid too many files errors
        reslistnew = []
        npart = 0
        
        for ic, pc in enumerate(reslist):
            nc = pc.replace('--od {od}'.format(od=outdir), '--od {od}/part{n}'.format(od=outdir,n=npart))
            reslistnew.append(nc)
            if not (ic+1)%6000:
                npart += 1
        reslist = reslistnew
        ## ---
    
        ## split the list into non-bkg and bkg_and_data
        bkglist = []
        datalist = []
        siglist = []
        for n,i in enumerate(reslist):
            #print "%d)\n %s\n\n" % (n,i)
            if '-o otherBkgHisto_' in i:
                bkglist.append(i)
            elif '-o dataHisto_' in i:
                datalist.append(i)
            else:
                siglist.append(i)
    
        print "\n#sig = %s\n#bkg = %d\n#data = %d" %(len(siglist),len(bkglist),len(datalist))
        print "Note: #sig is the number of real signal jobs when option -s was specified, otherwise it is simply the number of jobs that are not data or other backgrounds (currently they include W/Z->tau and anti-signal) \n"

        ## get the number of total submit jobs from the bunching
        nj = len(siglist)
        print 'full number of python commands to submit', nj+len(bkglist)+len(datalist)
        print '   ... grouping #sig into bunches of', options.groupJobs
        
        if not nj%options.groupJobs:
            njobs = int(nj/options.groupJobs)
        else: 
            njobs = int(nj/options.groupJobs) + 1
        ## ---
        
        jobdir = outdir+'/jobs/'
        logdir = outdir+'/logs/'
        outdirCondor = outdir+'/outs/'
        errdir = outdir+'/errs/'

        subcommands = []
    
        sourcefiles = []
        ## a bit awkward, but this keeps the bkg and data jobs separate. do the backgrounds first for speedup
        for ib in bkglist:
            pm = 'plus' if '-A alwaystrue positive' in ib else 'minus'
            tmp_srcfile_name, tmp_srcfile = getShFile(jobdir, 'bkg_'+pm)
            tmp_srcfile.write(ib)
            if options.automaticResubmission: 
                printReturnValueInScript(tmp_srcfile, tmp_srcfile_name)
            tmp_srcfile.close()
            sourcefiles.append(tmp_srcfile_name)

        for ib in datalist:
            pm = 'plus' if '-A alwaystrue positive' in ib else 'minus'
            era = ib.split("data := data_")[1].split("'")[0]
            tmp_srcfile_name, tmp_srcfile = getShFile(jobdir, 'data_'+era+'_'+pm)
            tmp_srcfile.write(ib)
            if options.automaticResubmission: 
                printReturnValueInScript(tmp_srcfile, tmp_srcfile_name)
            tmp_srcfile.close()
            sourcefiles.append(tmp_srcfile_name)
    
        ## now do the others.
        tag = "sig_" if options.signalCards else "otherSig_"
        for ij in range(njobs):
            tmp_srcfile_name, tmp_srcfile = getShFile(jobdir, tag+str(ij))
            tmp_n = options.groupJobs
            while len(siglist) and tmp_n:
                tmp_pycmd = siglist[0]
                tmp_srcfile.write(tmp_pycmd+'\n')
                siglist.remove(tmp_pycmd)
                tmp_n -= 1
                if options.automaticResubmission: 
                    printReturnValueInScript(tmp_srcfile, tmp_srcfile_name)
            tmp_srcfile.close()
            sourcefiles.append(tmp_srcfile_name)
        ## ---
    
        ## now on to the usual condor magic
        dummy_exec = open(jobdir+'/dummy_exec.sh','w')
        dummy_exec.write('#!/bin/bash\n')
        dummy_exec.write('script=$1\n')
        dummy_exec.write('bash $script\n')
        dummy_exec.write('retval=$?\n')
        # make exit code for success explicit, it can be used by condor to resubmit jobs
        if options.automaticResubmission:
            dummy_exec.write('attempt=$2\n')
            dummy_exec.write('echo "attempt = ${attempt}"\n')
        dummy_exec.write('echo "script = ${script}"\n')
        dummy_exec.write('echo "final return value = ${retval}"\n')
        dummy_exec.write('exit $retval\n')
        dummy_exec.close()
    
        condor_file_name = jobdir+'/condor_submit_'+('background' if options.bkgdataCards else 'signal')+'.condor'
        condor_file = open(condor_file_name,'w')
        condor_file.write('Universe = vanilla\n')
        condor_file.write('Executable = {de}\n'.format(de=os.path.abspath(dummy_exec.name)))
        condor_file.write('use_x509userproxy = true\n')
        # keep one .log in case of resubmission, text is appended there
        # but create a different .err and .out when resubmitting failures automatically
        condor_file.write('''
Log         = {ld}/{sb}_$(ClusterId)_$(ProcId).log
Output      = {od}/{sb}_$(ClusterId)_$(ProcId){text}.out
Error       = {ed}/{sb}_$(ClusterId)_$(ProcId){text}.error
getenv      = True
environment = "LS_SUBCWD={here}"
request_memory = 2000
+MaxRuntime = {rt}\n'''.format(ld=os.path.abspath(logdir), od=os.path.abspath(outdirCondor), 
                               ed=os.path.abspath(errdir), 
                               sb="sig" if options.signalCards else "bkg",
                               text="_trial$$([NumJobStarts])" if options.automaticResubmission else "",
                               rt=getCondorTime(options.queue), here=os.environ['PWD'] ) )
        ## ---
    
        if options.jobName != None:
            condor_file.write('+JobBatchName = "{n}"\n'.format(n=options.jobName))
        ## special permissions for people in CMG
        if os.environ['USER'] in ['mdunser', 'psilva', 'bendavid', 'kelong']:
            condor_file.write('+AccountingGroup = "group_u_CMST3.all"\n')
        if os.environ['USER'] in ['mciprian']:
            condor_file.write('+AccountingGroup = "group_u_CMS.CAF.ALCA"\n')
        condor_file.write('\n\n')

        # can retry up tp N times, if the job fails
        # when it fails, keep it on hold for 1 minute before releasing it again 
        # (in case the error was due to temporary problems in either the filesystem or eos itself)
        if options.automaticResubmission: 
            condor_file.write('''
max_retries = {n}
success_exit_code = 0
retry_until = (ExitCode == 0) && (NumJobStarts < {n1})
+OnExitHold   = (ExitCode != 0) && (NumJobStarts < {n1})
periodic_release =  (NumJobStarts < {n1}) && ((time() - EnteredCurrentStatus) > 60)\n\n'''.format(n=options.nResubmissions,n1=options.nResubmissions+1))

        for sf in sourcefiles:
            condor_file.write('arguments = {sf} {trial}\nqueue 1 \n\n'.format(sf=os.path.abspath(sf),trial="$$([NumJobStarts])" if options.automaticResubmission else ""))
        condor_file.close()
        ## ---
    
        ## now just finish up, either submitting or just printing the commands
        print 'I have {n} jobs to submit!'.format(n=len(subcommands))
        if options.dryRun:
            print 'running dry, printing the commands...'
            print 'condor_submit ',condor_file_name
            #for cmd in subcommands:
            #    print cmd
        else:
            print 'submitting for real...'
            os.system('condor_submit '+condor_file_name)
        ## ---

    print 'done'
