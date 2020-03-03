#!/usr/bin/env python
from shutil import copyfile
import re, sys, os, os.path, subprocess, json, ROOT
import numpy as np
# import some parameters from wmass_parameters.py, they are also used by other scripts
from wmass_parameters import *


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


NPDFSYSTS=60 # Hessian variations of NNPDF 3.0
nominals=[] # array containing the nominal for single process for which we have dedicated corrections (not included in bkg_and_data)
pdfsysts=[] # array containing the PDFs signal variations
qcdsysts=[] # array containing the QCD scale signal variations
inclqcdsysts=[] # array containing the inclusive QCD scale variations for Z
etaeffsysts=[] # array containing the uncorrelated efficiency systematics vs eta
fsrsysts=[] # array containing (only) the PHOTOS/PYTHIA fsr reweighitng
kamucasysts=[] # array containing the kalman filter muon momentum scale systematics
#coefficients = ['ac'] + ['a'+str(i) for i in range(8)]

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

def writeNominalSingleMCA(mcafile,odir,incl_mca='incl_dy'):
    incl_file=getMcaIncl(mcafile,incl_mca)
    proc = incl_mca.split('_')[1]
    postfix = "_{proc}_nominal".format(proc=proc)
    mcafile_nominal = open("{odir}/mca{pfx}.txt".format(odir=odir,pfx=postfix), "w")
    mcafile_nominal.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="1." \n')
    nominals.append(proc)
    print "written nominal mca relative to ",incl_mca

def writePdfSystsToMCA(mcafile,odir,vec_weight="hessWgt",syst="pdf",incl_mca='incl_sig',append=False):
    open("%s/systEnv-dummy.txt" % odir, 'a').close()
    incl_file=getMcaIncl(mcafile,incl_mca)
    if len(incl_file)==0: 
        print "Warning! '%s' include directive not found. Not adding pdf systematics samples to MCA file" % incl_mca
        return
    if append:
        filename = "%s/mca_systs.txt" % odir
        if not os.path.exists(filename): os.system('cp {mca_orig} {mca_syst}'.format(mca_orig=mcafile,mca_syst=filename))
    for i in range(1,NPDFSYSTS+1):
        postfix = "_{proc}_{syst}{idx}".format(proc=incl_mca.split('_')[1],syst=syst,idx=i)
        mcafile_syst = open(filename, 'a') if append else open("%s/mca%s.txt" % (odir,postfix), "w")
        mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+vec_weight+str(i)+'", PostFix="'+postfix+'" \n')
        pdfsysts.append(postfix)
    print "written ",syst," systematics relative to ",incl_mca


## this function here is pretty important to do all the theory systematics
## ---
def writeQCDScaleSystsToMCA(mcafile,odir,syst="qcd",incl_mca='incl_sig',scales=[],append=False,ptBinned=True,overrideDecorrelation=False):
    open("%s/systEnv-dummy.txt" % odir, 'a').close()
    incl_file=getMcaIncl(mcafile,incl_mca)

    if len(incl_file)==0: 
        print "Warning! '%s' include directive not found. Not adding QCD scale systematics!"
        return

    if append:
        filename = "%s/mca_systs.txt" % odir
        if not os.path.exists(filename): os.system('cp {mca_orig} {mca_syst}'.format(mca_orig=mcafile,mca_syst=filename))    

    ## make the mcas for the scales and mW and stuff
    for scale in scales:

        ## mW now doesn't just have Up and Down. make many
        if scale == "mW":
            ## loop on all the masses we want. in index of 5 MeV variations each
            masses  = ['mass_{m}'.format(m = j).replace('-','m') for j in range(-20,0)]
            masses += ['mass_0']
            masses += ['mass_p{m}'.format(m = j).replace('-','m') for j in range(1,21)]

            for mass in masses:
                postfix = "_{proc}_{syst}{mval}".format(proc=incl_mca.split('_')[1],syst=scale,mval=mass)
                mcafile_syst = open(filename, 'a') if append else open("%s/mca%s.txt" % (odir,postfix), "w")

                ## central mass is 80419 MeV, the closest we have to that is 80420. will scale +- 50 MeV, i.e. 80470 for Up and 80370 for Dn
                fstring = str(mass)
                mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+fstring+'", PostFix="'+postfix+'" \n')
                qcdsysts.append(postfix)

        else:
            ## all the others as usual with an Up and Down variation
            for idir in ['Up','Dn']:
                postfix = "_{proc}_{syst}{idir}".format(proc=incl_mca.split('_')[1],syst=scale,idir=idir)
                mcafile_syst = open(filename, 'a') if append else open("%s/mca%s.txt" % (odir,postfix), "w")
                if 'muR' in scale or 'muF' in scale:
    
                    ## for the signal keep as is for the QCD scales
                    if ptBinned:
                        for ipt in range(1,11): ## start from 1 to 10
                            for coeff in coefficients if options.decorrelateSignalScales and not overrideDecorrelation else ['']:
                                for pm in ['plus', 'minus']:
                                    ## have to redo the postfix for these
                                    postfix = "_{proc}_{coeff}{syst}{ipt}{ch}{idir}".format(proc=incl_mca.split('_')[1],syst=scale,idir=idir,ipt=ipt,coeff=coeff,ch=pm)
                                    mcafile_syst = open(filename, 'a') if append else open("%s/mca%s.txt" % (odir,postfix), "w")
                                    ptcut = wptBinsScales(ipt)
                                    wgtstr = 'TMath::Power(qcd_{sc}{idir}\,({wv}_pt>={ptlo}&&{wv}_pt<{pthi}))'.format(sc=scale,idir=idir,wv=options.wvar,ptlo=ptcut[0],pthi=ptcut[1])
                                    mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+wgtstr+'", PostFix="'+postfix+'" \n')
                                    qcdsysts.append(postfix)
                    else:
                        ## for the Z only do the unnumbered ones
                        postfix = "_{proc}_{syst}{idir}".format(proc=incl_mca.split('_')[1],syst=scale,idir=idir)
                        mcafile_syst = open(filename, 'a') if append else open("%s/mca%s.txt" % (odir,postfix), "w")
                        wgtstr = 'qcd_{sc}{idir}'.format(sc=scale,idir=idir)
                        mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+wgtstr+'", PostFix="'+postfix+'" \n')
                        inclqcdsysts.append(postfix)
    
                else: ## alphaS is left here. keep as is i guess
                    ## don't do mW here again
                    if "mW" in scale:
                        continue
                    ## ---
                    mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="qcd_'+scale+idir+'", PostFix="'+postfix+'" \n')
                    qcdsysts.append(postfix)
                    inclqcdsysts.append(postfix)
    
    print "written ",syst," systematics relative to ",incl_mca

def writeEfficiencyStatErrorSystsToMCA(mcafile,odir,channel,syst="EffStat",incl_mca='incl_sig',append=False):
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
    for i in range(0,nbins):
        etamin=etalo + i*deta; etamax=etalo + (i+1)*deta;
        # 3 parameters used in the Erf fit vs pt per eta bin
        for ipar in xrange(3):
            postfix = "_{proc}_ErfPar{ipar}{syst}{idx}".format(proc=incl_mca.split('_')[1],syst=syst,ipar=ipar,idx=i+1)
            mcafile_syst = open(filename, 'a') if append else open("%s/mca%s.txt" % (odir,postfix), "w")
            weightFcn = 'effSystEtaBins({ipar}\,LepGood1_pdgId\,LepGood1_eta\,LepGood1_pt\,{emin:.1f}\,{emax:.1f})'.format(ipar=ipar,emin=etamin,emax=etamax)
            mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+weightFcn+'", PostFix="'+postfix+'" \n')
            etaeffsysts.append(postfix)
    print "written ",syst," systematics relative to ",incl_mca

def writeFSRSystsToMCA(mcafile,odir,syst="fsr",incl_mca='incl_sig',append=False):
    open("%s/systEnv-dummy.txt" % odir, 'a').close()
    incl_file=getMcaIncl(mcafile,incl_mca)
    postfix = "_{proc}_{syst}".format(proc=incl_mca.split('_')[1],syst=syst)
    mcafile_syst = open("%s/mca%s.txt" % (odir,postfix), "w")
    weightFcn = 'fsrPhotosWeightSimple(GenLepDressed_pdgId[0]\,GenLepDressed_pt[0]\,GenLepBare_pt[0])'
    mcafile_syst.write(incl_mca+postfix+'   : + ; IncludeMca='+incl_file+', AddWeight="'+weightFcn+'", PostFix="'+postfix+'" \n')
    fsrsysts.append(postfix)
    print "written ",syst," systematics relative to ",incl_mca

def writePdfSystsToSystFile(filename,sample="W.*",syst="CMS_W_pdf"):
    SYSTFILEALL=('.').join(filename.split('.')[:-1])+"-all.txt"
    copyfile(filename,SYSTFILEALL)
    systfile=open(SYSTFILEALL,"a")
    for i in range(NPDFSYSTS/2):
        systfile.write(syst+str(i+1)+"  : "+sample+" : .* : pdf"+str(i+1)+" : templates\n")
    print "written pdf syst configuration to ",SYSTFILEALL
    return SYSTFILEALL

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
    parser.add_option("-l", "--lumi",    dest="integratedLuminosity",     type="float", default=35.9, help="Integrated luminosity");
    parser.add_option("--dry-run", dest="dryRun",    action="store_true", default=False, help="Do not run the job, only print the command");
    parser.add_option("--decorrelateSignalScales", action="store_true", default=False, help="Treat qcd scales uncorrelated for left/right/long");
    parser.add_option("-s", "--signal-cards",  dest="signalCards",  action="store_true", default=False, help="Make the signal part of the datacards");
    parser.add_option("-b", "--bkgdata-cards", dest="bkgdataCards", action="store_true", default=False, help="Make the background and data part of the datacards");
    parser.add_option("-W", "--weight", dest="weightExpr", default="-W 1", help="Event weight expression (default 1)");
    parser.add_option("-P", "--path", dest="path", type="string",default=None, help="Path to directory with input trees and pickle files");
    parser.add_option("-C", "--channel", dest="channel", type="string", default='el', help="Channel. either 'el' or 'mu'");
    parser.add_option("--not-unroll2D", dest="notUnroll2D", action="store_true", default=False, help="Do not unroll the TH2Ds in TH1Ds needed for combine (to make 2D plots)");
    parser.add_option("--pdf-syst", dest="addPdfSyst", action="store_true", default=False, help="Add PDF systematics to the signal (need incl_sig directive in the MCA file)");
    parser.add_option("--qcd-syst", dest="addQCDSyst", action="store_true", default=False, help="Add QCD scale systematics to the signal (need incl_sig directive in the MCA file)");
    parser.add_option("--qed-syst", dest="addQEDSyst", action="store_true", default=False, help="Add QED scale systematics to the signal (need incl_sig directive in the MCA file)");
    parser.add_option("--kamuca-syst", dest="addKaMuCaSyst", action="store_true", default=False, help="Add Kalman Filter muon momentum scale correction systematics");
    parser.add_option('-g', "--group-jobs", dest="groupJobs", type=int, default=20, help="group signal jobs so that one job runs multiple makeShapeCards commands");
    parser.add_option('-w', "--wvar", type="string", default='prefsrw', help="switch between genw and prefsrw. those are the only options (default %default)");
    parser.add_option('--vpt-weight', dest='procsToPtReweight', action="append", default=[], help="processes to be reweighted according the measured/predicted DY pt. Default is none (possible W,TauDecaysW,Z).");
    parser.add_option('--wlike', dest='wlike', action="store_true", default=False, help="Make cards for the wlike analysis. Default is wmass");
    parser.add_option('--inclusive', dest='inclusive', action="store_true", default=False, help="do not reweight angular coefficients. make it all inclusive.");
    (options, args) = parser.parse_args()
    
    if not options.wvar in ['genw', 'prefsrw']:
        print 'the W variable has to be either "genw" or "prefsrw". exiting...'
        quit()
    
    if len(sys.argv) < 6:
        parser.print_usage()
        quit()
        
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

    coefficients = ['inc'] if options.inclusive else ['ac'] + ['a'+str(i) for i in range(8)]
    
    luminosity = options.integratedLuminosity
    
    if not os.path.exists("cards/"):
        os.makedirs("cards/")
    outdir="cards/"+args[5]
    
    if not os.path.exists(outdir): os.mkdir(outdir)
    if options.queue and not os.path.exists(outdir+"/jobs"): 
        os.mkdir(outdir+"/jobs")
        os.mkdir(outdir+"/mca")
    
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
        writePdfSystsToMCA(MCA,outdir+"/mca") # on W + jets 
        writePdfSystsToMCA(MCA,outdir+"/mca",incl_mca='incl_dy') # on DY + jets
        writePdfSystsToMCA(MCA,outdir+"/mca",incl_mca='incl_wtau') # on WTau
        # write the complete systematics file (this was needed when trying to run all systs in one job)
        # SYSTFILEALL = writePdfSystsToSystFile(SYSTFILE)
    if options.addQCDSyst:
        scales = ['muR','muF',"muRmuF", "alphaS"]
        if options.wlike:
            wmasses = []
            zmasses = ["mW"]
            ptBinnedScalesForW = False
        else:
            wmasses = ["mW"]
            zmasses = []
            ptBinnedScalesForW = True
        writeQCDScaleSystsToMCA(MCA,outdir+"/mca",scales=scales+wmasses,incl_mca='incl_sig' ,ptBinned=ptBinnedScalesForW)
        writeQCDScaleSystsToMCA(MCA,outdir+"/mca",scales=scales+zmasses,incl_mca='incl_dy'  ,ptBinned=False)
        writeQCDScaleSystsToMCA(MCA,outdir+"/mca",scales=scales        ,incl_mca='incl_wtau',ptBinned=False,overrideDecorrelation=True)
    if options.addQEDSyst:
        if options.wlike:
            writeFSRSystsToMCA(MCA,outdir+"/mca",incl_mca='incl_dy') # DY + jets
        else:
            writeFSRSystsToMCA(MCA,outdir+"/mca") # on W + jets
    
    if len(pdfsysts+qcdsysts)>1:
        for proc in ['dy','wtau']:
            if proc not in nominals: writeNominalSingleMCA(MCA,outdir+"/mca",incl_mca='incl_{p}'.format(p=proc))
    
    if options.addKaMuCaSyst:
        addKaMuCaSysts()

    # not needed if we fit eta/pt. Will be needed if we fit another variable correlated with eta/pt
    # writeEfficiencyStatErrorSystsToMCA(MCA,outdir+"/mca",options.channel)
    
    ARGS=" ".join([MCA,CUTFILE,"'"+fitvar+"' "+"'"+binning+"'",SYSTFILE])
    BASECONFIG=os.path.dirname(MCA)
    ## use rel paths if options.queue:
    ## use rel paths     ARGS = ARGS.replace(BASECONFIG,os.getcwd()+"/"+BASECONFIG)
    OPTIONS=" -P "+T+" --s2v -j "+str(J)+" -l "+str(luminosity)+" -f --obj tree "+FASTTEST
    if not os.path.exists(outdir): os.makedirs(outdir)
    OPTIONS+=" -F Friends '{P}/friends/tree_Friend_{cname}.root' "
    if not options.notUnroll2D:
        OPTIONS+=" --2d-binning-function unroll2Dto1D "
    

    if options.wlike:
        POSCUT=" -A alwaystrue positive 'evt%2 != 0' "
        NEGCUT=" -A alwaystrue negative 'evt%2 == 0' "    
        SIGPROC='Z'
        SIGSUFFIX='dy'
    else:
        POSCUT=" -A alwaystrue positive 'LepGood1_charge>0' "
        NEGCUT=" -A alwaystrue negative 'LepGood1_charge<0' "        
        SIGPROC='W'
        SIGSUFFIX='sig'
    antiSIGPROC = 'Z' if SIGPROC=='W' else 'W'

    print 'these are angular coefficients ', coefficients

    fullJobList = set()
    if options.signalCards:
        print "MAKING SIGNAL PART!"
        wsyst = ['']+[x for x in pdfsysts+qcdsysts+inclqcdsysts+etaeffsysts+fsrsysts if SIGSUFFIX in x]
        wsyst += kamucasysts
        ## loop on all the variations
        for ivar,var in enumerate(wsyst):
            SYST_OPTION = ''
            ## loop on all the angular coefficients
            for coeff in coefficients:
                ## get the list of all other coefficients
                anticoeff = [proc for proc in coefficients if not proc==coeff]
                if options.inclusive:
                    anticoeff = ['ac'] + ['a'+str(i) for i in range(8)] 
                if any(i in var for i in anticoeff) and options.decorrelateSignalScales: continue ## this might work. but who knows...
                ## loop on both charges
                for charge in ['plus', 'minus']:
                    antich = 'plus' if charge == 'minus' else 'minus'
                    ## 0 is nominal, all the theory variations are non-zero
                    if ivar==0: 
                        IARGS = ARGS
                    else: 
                        if 'kalPt' in var:
                            IARGS = re.sub(r'LepGood(\d+)_pt',r'LepGood\g<1>_kalPt{systvar}'.format(systvar=var.replace('kalPt','')),ARGS)
                            SYST_OPTION = ' --replace-var LepGood1_pt LepGood1_{systvar} --replace-var LepGood2_pt LepGood2_{systvar} '.format(systvar=var)
                        else:
                            IARGS = ARGS.replace(MCA,"{outdir}/mca/mca{syst}.txt".format(outdir=outdir,syst=var))
                            IARGS = IARGS.replace(SYSTFILE,"{outdir}/mca/systEnv-dummy.txt".format(outdir=outdir))
                        print "Running the systematic: ",var
                    ## ---
                    job_group =  [] ## group

                    ## now go and make the submit commands
                    print "Making card for signal process {coeff} with charge {ch} ".format(coeff=coeff,ch=charge)
                    ycut = POSCUT if charge=='plus' else NEGCUT
                    ## exclude the anti charge and the anti coefficients
                    if antich in var: 
                        print "skipping because ",antich," is in ",var
                        continue
                    excl_anticoeff = ','.join(SIGPROC+charge+'_'+ac+'.*' for ac in anticoeff)
                    xpsel=' --xp "{sig}{antich}.*,{ahel},Flips,{antisig}.*,Top,DiBosons,TauDecaysW.*,data.*" --asimov '.format(sig=SIGPROC,antisig=antiSIGPROC,antich=antich,ch=charge,ahel=excl_anticoeff)
                    ## ---

                    ## make necessary directories
                    if not os.path.exists(outdir): 
                        os.mkdir(outdir)
                    if options.queue and not os.path.exists(outdir+"/jobs"): 
                        os.mkdir(outdir+"/jobs")
                    ## ---

                    ## make specific names and such
                    syst = '' if ivar==0 else var
                    ## for kamuca corrections, there is not any _sigsuffix, so add an "_" to distinguish better the files
                    if 'kalPt' in var: syst = '_{sfx}_{var}'.format(sfx=SIGSUFFIX,var=var)
                    dcname = "{sig}{charge}_{coeff}_{channel}{syst}".format(sig=SIGPROC, charge=charge, coeff=coeff, channel=options.channel,syst=syst)
                    
                    ## reweight the boson-pT here
                    zptWeight = 'dyptWeight({wvar}_pt,{isZ})'.format(wvar=options.wvar,isZ=0 if SIGPROC=='W' else 1)

                    ## construct the full weight to be applied
                    fullWeight = options.weightExpr+'*'+zptWeight if SIGPROC in options.procsToPtReweight else options.weightExpr
                    BIN_OPTS = OPTIONS + " -W '" + fullWeight+ "'" + " -o "+dcname+" --od "+outdir + xpsel + ycut + SYST_OPTION

                    ## do not do by pdf weight pdfmatch = re.search('pdf(\d+)',var)
                    ## do not do by pdf weight if pdfmatch:
                    ## do not do by pdf weight     BIN_OPTS += " -W 'pdfRatioHel(abs(prefsrw_y),prefsrw_pt,prefsrw_costcs,evt,{pol},{ipdf})' ".format(pol=helmap[helicity],ipdf=pdfmatch.group(1))
                    ## ---

                    ## now we accumulate all the commands to be run for all processes
                    if options.queue:
                        mkShCardsCmd = "python makeShapeCards.py {args} \n".format(dir = os.getcwd(), args = IARGS+" "+BIN_OPTS)
                        fullJobList.add(mkShCardsCmd)

                    ## i would suggest not ever going into this else statement ... 
                    else:
                        cmd = "python makeShapeCards.py "+IARGS+" "+BIN_OPTS
                        if options.dryRun: print cmd
                        else:
                            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
                            out, err = p.communicate() 
                            result = out.split('\n')
                            for lin in result:
                                if not lin.startswith('#'):
                                    print(lin)
                    ## ---
    
    ## this is for all the background cards, excluding the Z(W) for Wmass(Wlike)
    if options.bkgdataCards:
        print "MAKING BKG and DATA PART:\n"
        ## again loop on both charges
        for charge in ['plus','minus']:
            ## remove W or Z signal processes and others that aren't needed
            xpsel=' --xp "{sig}.*" '.format(sig=SIGPROC)

            if len(pdfsysts+qcdsysts)>1: # 1 is the nominal 
                xpsel+=' --xp "{antisig}.*,TauDecaysW.*" '.format(antisig=antiSIGPROC)   # adding .* to tau will be necessary if we ever decide to use lepeff and XXscale on that as well
            ## ---

            ## now make the names of the cards etc
            chargecut = POSCUT if charge=='plus' else NEGCUT
            dcname = "bkg_and_data_{channel}_{charge}".format(channel=options.channel, charge=charge)
            BIN_OPTS=OPTIONS + " -W '" + options.weightExpr + "'" + " -o "+dcname+" --od "+outdir + xpsel + chargecut

            ## accumulate the commands to be run
            if options.queue:
                mkShCardsCmd = "python makeShapeCards.py {args} \n".format(dir = os.getcwd(), args = ARGS+" "+BIN_OPTS)
                fullJobList.add(mkShCardsCmd)

            ## i would suggest not ever going into this else statement ... 
            else:
                cmd = "python makeShapeCards.py "+ARGS+" "+BIN_OPTS
                if options.dryRun: print cmd
                else:
                    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
                    out, err = p.communicate() 
                    result = out.split('\n')
                    for lin in result:
                        if not lin.startswith('#'):
                            print(lin)
            ## ---
    
    ## this is for background cards which include the Z(W) samples for Wmass (Wlike)
    if options.bkgdataCards and len(pdfsysts+inclqcdsysts)>1:
        antiSIGSUFFIX = 'dy' if SIGPROC=='W' else 'sig'
        antisig_syst = ['_'+antiSIGSUFFIX+'_nominal']+[x for x in pdfsysts+inclqcdsysts if antiSIGSUFFIX in x] # for the Z in W mass the QCDscales are unbinned, for Wlike the W QCD scales are also unbinned => only need inclqcdsysts

        ## loop on the Z(W) related systematics for Wmass (Wlike)
        for ivar,var in enumerate(antisig_syst):
            ## loop on both charges
            for charge in ['plus','minus']:
                antich = 'plus' if charge == 'minus' else 'minus'
                
                ## nominal first, then the variations
                if ivar==0: 
                    # the following would be faster with DY-only, but it misses the lines for the Z_lepeff systs
                    # IARGS = ARGS.replace(MCA,"{outdir}/mca/mca_dy_nominal.txt".format(outdir=outdir))
                    IARGS = ARGS
                else: 
                    IARGS = ARGS.replace(MCA,"{outdir}/mca/mca{syst}.txt".format(outdir=outdir,syst=var))
                    IARGS = IARGS.replace(SYSTFILE,"{outdir}/mca/systEnv-dummy.txt".format(outdir=outdir))
                    print "Running the DY with systematic: ",var

                ## ---
                print "Making card for ",antiSIGSUFFIX, " process with charge ", charge
                chcut = POSCUT if charge=='plus' else NEGCUT
                xpsel=' --xp "[^{antisig}]*" --asimov '.format(antisig=antiSIGPROC)
                syst = '' if ivar==0 else var

                ## make names for the files and outputs
                dcname = "{antisig}_{channel}_{charge}{syst}".format(antisig=antiSIGPROC,channel=options.channel, charge=charge,syst=syst)

                ## construct the final weight. reweight also the DY to whatever new pT spectrum we want
                zptWeight = 'dyptWeight({wvar}_pt,{isZ})'.format(wvar=options.wvar,isZ=1 if SIGPROC=='W' else 0)
                fullWeight = options.weightExpr+'*'+zptWeight if antiSIGPROC in options.procsToPtReweight else options.weightExpr
                BIN_OPTS=OPTIONS + " -W '" + fullWeight + "'" + " -o "+dcname+" --od "+outdir + xpsel + chcut
                ## ---

                ## accumulate the commands to be run
                if options.queue:
                    mkShCardsCmd = "python makeShapeCards.py {args} \n".format(dir = os.getcwd(), args = IARGS+" "+BIN_OPTS)
                    fullJobList.add(mkShCardsCmd)
                ## i would suggest not ever going into this else statement ... 
                else:
                    cmd = "python makeShapeCards.py "+IARGS+" "+BIN_OPTS
                    if options.dryRun: print cmd
                    else:
                        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
                        out, err = p.communicate() 
                        result = out.split('\n')
                        for lin in result:
                            if not lin.startswith('#'):
                                print(lin)
                ## ---
    
    ## repetition for WTau, but better to keep Z/Tau cards separated (faster jobs)
    if options.bkgdataCards and len(pdfsysts+qcdsysts)>1:
        wtausyst = ['_wtau_nominal']+[x for x in pdfsysts+qcdsysts+inclqcdsysts if 'wtau' in x] # for the W mass the QCD scales are binned, for Wlike they are unbinned, so need qcdsysts+inclqcdsysts to be general
        ## loop on systmatic variations
        for ivar,var in enumerate(wtausyst):
            ## loop on both charges
            for charge in ['plus','minus']:
                antich = 'plus' if charge == 'minus' else 'minus'

                ## 0 is nominal again
                if ivar==0: 
                    IARGS = ARGS.replace(MCA,"{outdir}/mca/mca_wtau_nominal.txt".format(outdir=outdir))
                else: 
                    IARGS = ARGS.replace(MCA,"{outdir}/mca/mca{syst}.txt".format(outdir=outdir,syst=var))
                    IARGS = IARGS.replace(SYSTFILE,"{outdir}/mca/systEnv-dummy.txt".format(outdir=outdir))
                    print "Running the WTau with systematic: ",var

                print "Making card for WTau process with charge ", charge
                chcut = POSCUT if charge=='plus' else NEGCUT
                xpsel=' --xp "[^TauDecaysW]*" --asimov '
                syst = '' if ivar==0 else var

                ## make names for files etc again
                dcname = "TauDecaysW_{channel}_{charge}{syst}".format(channel=options.channel, charge=charge,syst=syst)
                ## construct full reweighting weight. boson-pT reweighting here again
                zptWeight = 'dyptWeight({wvar}_pt,0)'.format(wvar=options.wvar)
                fullWeight = options.weightExpr+'*'+zptWeight if 'TauDecaysW' in options.procsToPtReweight else options.weightExpr
                BIN_OPTS=OPTIONS + " -W '" + fullWeight + "'" + " -o "+dcname+" --od "+outdir + xpsel + chcut
                ## ---

                ## accumulate the commands to be run
                if options.queue:
                    mkShCardsCmd = "python makeShapeCards.py {args} \n".format(dir = os.getcwd(), args = IARGS+" "+BIN_OPTS)
                    fullJobList.add(mkShCardsCmd)
                ## i would suggest not ever going into this else statement ...
                else:
                    cmd = "python makeShapeCards.py "+IARGS+" "+BIN_OPTS
                    if options.dryRun: print cmd
                    else:
                        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
                        out, err = p.communicate() 
                        result = out.split('\n')
                        for lin in result:
                            if not lin.startswith('#'):
                                print(lin)
                ## ---
    
    ## not sure why this function is here, but who cares
    def getShFile(jobdir, name):
        tmp_srcfile_name = jobdir+'/job_{t}_{i}.sh'.format(i=name,t='sig' if options.signalCards else 'bkg')
        tmp_srcfile = open(tmp_srcfile_name, 'w')
        tmp_srcfile.write("#! /bin/sh\n")
        tmp_srcfile.write("ulimit -c 0 -S\n")
        tmp_srcfile.write("ulimit -c 0 -H\n")
        tmp_srcfile.write("cd {cmssw};\neval $(scramv1 runtime -sh);\ncd {d};\n".format( d= os.getcwd(), cmssw = os.environ['CMSSW_BASE']))
        return tmp_srcfile_name, tmp_srcfile
    
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
        bkglist = [i for i in reslist if     'bkg_and_data' in i]
        reslist = [i for i in reslist if not 'bkg_and_data' in i]
    
        ## get the number of total submit jobs from the bunching
        nj = len(reslist)
        print 'full number of python commands to submit', nj+len(bkglist)
        print '   ... grouping them into bunches of', options.groupJobs
        
        if not nj%options.groupJobs:
            njobs = int(nj/options.groupJobs)
        else: 
            njobs = int(nj/options.groupJobs) + 1
        ## ---
        
        jobdir = outdir+'/jobs/'
        os.system('mkdir -p '+jobdir)
        subcommands = []
    
        sourcefiles = []
        ## a bit awkward, but this keeps the bkg and data jobs separate. do the backgrounds first for speedup
        for ib in bkglist:
            pm = 'plus' if 'positive' in ib else 'minus'
            tmp_srcfile_name, tmp_srcfile = getShFile(jobdir, 'bkg_'+pm)
            tmp_srcfile.write(ib)
            tmp_srcfile.close()
            sourcefiles.append(tmp_srcfile_name)
    
        ## now do the others.
        for ij in range(njobs):
            tmp_srcfile_name, tmp_srcfile = getShFile(jobdir, ij)
            tmp_n = options.groupJobs
            while len(reslist) and tmp_n:
                tmp_pycmd = reslist[0]
                tmp_srcfile.write(tmp_pycmd)
                reslist.remove(tmp_pycmd)
                tmp_n -= 1
            tmp_srcfile.close()
            sourcefiles.append(tmp_srcfile_name)
        ## ---
    
        ## now on to the usual condor magic
        dummy_exec = open(jobdir+'/dummy_exec.sh','w')
        dummy_exec.write('#!/bin/bash\n')
        dummy_exec.write('bash $*\n')
        dummy_exec.close()
    
        condor_file_name = jobdir+'/condor_submit_'+('background' if options.bkgdataCards else 'signal')+'.condor'
        condor_file = open(condor_file_name,'w')
        condor_file.write('''Universe = vanilla
    Executable = {de}
    use_x509userproxy = true
    Log        = {jd}/$(ProcId).log
    Output     = {jd}/$(ProcId).out
    Error      = {jd}/$(ProcId).error
    getenv      = True
    environment = "LS_SUBCWD={here}"
    request_memory = 4000
    +MaxRuntime = {rt}\n
    '''.format(de=os.path.abspath(dummy_exec.name), jd=os.path.abspath(jobdir), rt=getCondorTime(options.queue), here=os.environ['PWD'] ) )
        ## ---

        ## special permissions for people in CMG
        if os.environ['USER'] in ['mdunser', 'psilva', 'bendavid', 'kelong']:
            condor_file.write('+AccountingGroup = "group_u_CMST3.all"\n\n\n')
        for sf in sourcefiles:
            condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(sf)))
        condor_file.close()
        ## ---
    
        ## now just finish up, either submitting or just printing the commands
        print 'i have {n} jobs to submit!'.format(n=len(subcommands))
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
