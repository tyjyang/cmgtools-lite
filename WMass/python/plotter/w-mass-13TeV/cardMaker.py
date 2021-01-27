#!/bin/env python

import ROOT, sys,os,re,json, copy, math, root_numpy, array
from rollingFunctions import roll1Dto2D, dressed2D, unroll2Dto1D
import makeSystematicMultipliers as ms
import utilities
utilities = utilities.util()
def isExcludedNuisance(excludeNuisances=[], name="", keepNuisances=[]):
    if len(excludeNuisances) and any(re.match(x,name) for x in excludeNuisances):
        if len(keepNuisances) and any(re.match(x,name) for x in keepNuisances):
            return False
        else:
            print ">>>>> Excluding nuisance: ", name
            return True
    else:
        return False

class CardMaker:
    def __init__(self,options,charge,flavor):
        self._options = options
        self.charge = charge
        self.flavor = flavor
        self.data_eras = ["B", "C", "D", "E", "F", "F_postVFP", "G", "H"]
        self.excludeNuisances = []
        self.keepNuisances = []
        if len(options.excludeNuisances):
            self.excludeNuisances = options.excludeNuisances.split(",")
        if len(options.keepNuisances):
            self.keepNuisances = options.keepNuisances.split(",")
                
        ### import the right dictionary process - systematic list depending on wmass/wlike
        #pathToImport = os.environ['CMSSW_BASE']+'/src/CMGTools/WMass/python/plotter/w-mass-13TeV/'
        pathToImport = os.path.dirname(sys.argv[0])+'/'
        if self._options.wmass:
            pathToImport += 'wmass_'+flavor
            self.boson = 'W{fl}nu'.format(fl=self.flavor)
        else:
            pathToImport += 'wlike_'+flavor
            self.boson = 'Z{fl}{fl}'.format(fl=self.flavor)
        self.systematics  = self.getSystList()
        self.centralfiles = self.getCentralProcesses()
        self.shapesfile = os.path.join(self._options.inputdir,self.boson+'_{ch}_shapes.root'.format(ch=charge))
        self.cardfile   = os.path.join(self._options.inputdir,self.boson+'_{ch}_card.txt'   .format(ch=charge))
        self.systFile = pathToImport+'/systsFit.txt'
        # self.centralHistograms = self.getCentralHistograms() # to implement


    def getCentralProcesses(self):        
        #acoeffs = ['ac']+['a{ic}'.format(ic=i) for i in range(8)]
        acoeffs = ['']
        # signal
        sig_proc          = ['{boson}_{charge}{ic}'.format(charge=self.charge,ic=("_"+coeff) if len(coeff) else "",boson=self.boson) for coeff in acoeffs]
        # antisignal
        other_boson_procs = ['{otherboson}_{charge}'.format(charge=self.charge,otherboson='Zmumu' if self._options.wmass else 'Wmunu')]
        # other stuff (for data use then global one (not split by era, it is created afterwards)
        other_files = ["Wtaunu","Ztautau","otherBkgHisto", "dataHisto"]
        otherprocs = ["{op}_{charge}".format(op=otherproc,charge=self.charge) for otherproc in other_files] 
        return sig_proc + other_boson_procs + otherprocs

    def getSystList(self):
        ### systematics that are applied to both W and Z, for mu and tau decays
        NPDFSYSTS=2
        baseSysts  = ['pdf{i}'.format(i=ipdf) for ipdf in range(1,1+NPDFSYSTS)]
        baseSysts += ['alphaS{idir}'.format(idir=idir) for idir in ['Up','Dn']]
        # FIXME: need a better way to decide if and how these are binned, depending on what was done in make_wmass_cards.py
        qcdSystsUnbin  = ['{qcdpar}{idir}'.format(qcdpar=par,idir=idir) for par in ['muR','muF','muRmuF'] for idir in ['Up','Dn']]     
        NVPTBINS = 2
        qcdSystsPtChargebin  = ['{qcdpar}{ipt}{ch}{idir}'.format(qcdpar=par,ipt=ipt,ch=self.charge,idir=idir) for par in ['muR','muF','muRmuF'] for ipt in range(1,1+NVPTBINS) for idir in ['Up','Dn']]     
        #baseSysts += ['kalPtErr{i}{idir}'.format(i=istat,idir=idir) for istat in range(133) for idir in ['Up','Dn']]
        #baseSysts += ['kalPtClosureErr{idir}'.format(idir=idir) for idir in ['Up','Dn']]
        ### W/Z mass points
        MASSVARIATIONS = [10* i for i in range(1,3)]
        massPoints = ['massShift{v}MeV{d}".format(v=mvar,d=idir)' for mvar in MASSVARIATIONS for idir in ['Up','Down']]
        ### other systematics
        # otherSysts = ['fsr']
        otherSysts = []
        ### build the systematic dictionary for the root files
        systsCards = {self.boson: baseSysts + qcdSystsPtChargebin + massPoints + otherSysts,
                      'Zmumu' if self._options.wmass else 'Wmunu': baseSysts + qcdSystsUnbin,
                      'Wtaunu': baseSysts + qcdSystsPtChargebin + massPoints,
                      'Ztautau': baseSysts + qcdSystsUnbin}
        return systsCards

    def applySystematics(self, dictWithCentrals, systEnv): #process, syst_regexp, syst_file, outfile_name):
        
        f_systs = ROOT.TFile(options.inputdir+'/systematicMultipliers.root', 'READ')
        systs_lok = f_systs.GetListOfKeys()
    
        allVars = []
        dictProcSyst = {}

        ## loop on all the systematics in the multipliers file
        for syst in systs_lok:
            systName = syst.GetName()
            if '_2d' in systName: ## forget about the 2d histograms
                continue
            ## loop on all the systematics in the systFit.txt file
            for dummy,(syst_regexp,proc_regexp,group,scale) in systEnv.items():
                if re.match(syst_regexp, systName):
                    tmp_syst = f_systs.Get(systName)
                    if 'TH2' in tmp_syst.ClassName():
                        continue
                    ## loop on all the processes in the dictWithCentrals
                    for (proc,cen),proc_hist in dictWithCentrals.items():
                        if ('plus' in proc and 'minus' in systName) or ('minus' in proc and 'plus' in systName):
                            continue
                        if re.match(proc_regexp,proc):
                            #print "hname = %s     systname = %s" % (proc_hist.GetName(),systName)
                            tmp_var  = copy.deepcopy(proc_hist.Clone(proc_hist.GetName()+'_'+systName))
                            tmp_var.Multiply(tmp_syst)
                            allVars.append(tmp_var)
                            dictProcSyst[(proc, systName)] = None
    
        outfile = ROOT.TFile(options.inputdir+'/{f}_{ch}_shapes.root.multiplierSystematics'.format(f=self.boson,ch=self.charge), 'RECREATE')
        
        for v in allVars:
            v.Write()
        outfile.Close()

        return dictProcSyst


    def mirrorShape(self,nominal,alternate,newname,alternateShapeOnly=False,use2xNomiIfAltIsZero=False):
        alternate.SetName("%sUp" % newname)
        if alternateShapeOnly:
            alternate.Scale(nominal.Integral()/alternate.Integral())
        mirror = nominal.Clone("%sDown" % newname)
        for b in xrange(1,nominal.GetNbinsX()+1):
            y0 = nominal.GetBinContent(b)
            yA = alternate.GetBinContent(b)
            # geometric mean
            # yM = y0
            # if yA != 0:
            #     yM = y0*y0/yA
            # elif yA == 0:
            #     if use2xNomiIfAltIsZero: 
            #         yM = 2. * y0
            #     else: 
            #         yM = 0
            # arithmetic mean
            yM = max(0,2*y0-yA)
            mirror.SetBinContent(b, yM)
        if alternateShapeOnly:
            # keep same normalization
            mirror.Scale(nominal.Integral()/mirror.Integral())
        return (alternate,mirror)

    def mergeChunks(self, dryRun=False):
        ## prepare the relevant files. only the datacards and the correct charge
        allfiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(self._options.inputdir,followlinks=True) for f in fn if f.endswith('.input.root')]
        
        datafiles = [f for f in allfiles if any('dataHisto_{e}_{ch}.input.root'.format(e=era,ch=self.charge) in f for era in self.data_eras)]
        datafolder = os.path.dirname(os.path.abspath(datafiles[0]))+"/"
        # start merging data (contains data and fakes)
        outf = "dataHisto_{ch}.input.root".format(ch=self.charge)
        haddcmd = "hadd -f -k {df}{of} {i}".format(df=datafolder,of=outf,i=" ".join(datafiles))
        print "Merging files with data and fakes"
        print haddcmd
        if not dryRun: os.system(haddcmd)

        ## get only the central processes
        files = [f for f in allfiles if any(os.path.basename(f)==centralfile+'.input.root' for centralfile in self.centralfiles)]
        files = sorted(files) #, key = lambda x: int(x.rstrip('.input.root').split('_')[1]))
    
        #acoeffs = [self.charge+'_ac']+['{charge}_a{ic}'.format(charge=self.charge,ic=i) for i in range(8)]

        tmpfiles = []
        for proc in self.centralfiles:
            print 'processing process: ',proc
            boson = proc.split("_")[0]
            bosonSysts = self.systematics[boson] if boson in self.systematics else []
            ## the input ROOT file is not built from the inputdir since it may come from different partX subdir
            central_input = [f for f in files if os.path.basename(f)==proc+'.input.root'][0]
            if len(bosonSysts):
                syst_inputs   = sorted([f for f in allfiles if any(os.path.basename(f)=='{proc}_{syst}.input.root'.format(proc=proc,syst=syst) for syst in bosonSysts)])
            else:
                syst_inputs = []

            nominals = {}
            for irf,rf in enumerate([central_input]+syst_inputs):
                print '\twith nominal/systematic file: ',rf
                tf = ROOT.TFile.Open(rf)
                tmpfile = os.path.join(self._options.inputdir,'tmp_{proc}_{sys}.root'.format(proc=proc,sys=irf))
                of=ROOT.TFile(tmpfile,'recreate')
                tmpfiles.append(tmpfile)
                # remove the duplicates also
                plots = {}
                for e in tf.GetListOfKeys() :
                    name=e.GetName()
                    obj=e.ReadObj()
                    # x_data must be removed, it was already copied as x_data_obs, which is what we need
                    if name == "x_data" and 'data' in proc: 
                        continue
                    # exclude the asimov data_obs which is in any single card, apart the file with the real data
                    if name.endswith('data_obs') and 'data' not in proc: 
                        continue
                    if name.endswith("Dn"):
                        name = name.replace('Dn','Down')                    
                    ### bkg_and_data is in central_input, with no systs. All the systematics (Up/Down) are just copied into the target file
                    if irf==0:
                        if name not in plots:
                            plots[name] = obj.Clone(name+"_clone")
                            nominals[name] = obj.Clone(name+"0")
                            nominals[name].SetDirectory(None)
                            #print 'replacing old %s with %s' % (name,name)
                            plots[name].Write(name)
                    else:
                        if any(sysname in name for sysname in ['pdf','fsr']): # these changes by default shape and normalization. Each variation should be symmetrized wrt nominal
                            pfx = '_'.join(name.split("_")[:-1]) # was -2, but names changed
                            if 'pdf' in name:
                                (alternate,mirror) = self.mirrorShape(nominals[pfx],obj,name,alternateShapeOnly=self._options.pdfShapeOnly)
                            else:
                                (alternate,mirror) = self.mirrorShape(nominals[pfx],obj,name,alternateShapeOnly=True)
                            for alt in [alternate,mirror]:
                                if alt.GetName() not in plots:
                                    plots[alt.GetName()] = alt.Clone(alt.GetName()+"_clone")
                                    plots[alt.GetName()].Write(alt.GetName())
                        elif re.match('.*massShift.*',name):
                            plots[name] = obj.Clone(name+"_clone")
                            plots[name].Write(name)
                        else:
                            if name not in plots:
                                plots[name] = obj.Clone(name+"_clone")
                                plots[name].Write(name)
                of.Close()

        if self._options.mergeRoot:
            haddcmd = 'hadd -f {of}.baseSystematics {indir}/tmp_*.root'.format(of=self.shapesfile, indir=self._options.inputdir )
            os.system(haddcmd)
            os.system('rm {indir}/tmp_*.root'.format(indir=self._options.inputdir))
            
        print "DONE adding input chunks."        


    def activateSystMatrix(self,systFile):
        activeSysts = {}
        for line in open(systFile, 'r'):
            if re.match("\s*#.*", line): continue
            line = re.sub("#.*","",line).strip()
            if len(line) == 0: continue
            field = [f.strip() for f in line.split(':')]
            if len(field) < 4:
                raise RuntimeError, "Malformed line %s in file %s"%(line.strip(),systFile)
            elif len(field) == 4:
                (name, sytPatt, procPatt, group, scaling) = field[:4] + [1]
            else:
                (name, sytPatt, procPatt, group, scaling) = field[:5]
            activeSysts[name] = (re.compile(sytPatt+"$"),re.compile(procPatt+"$"),group,scaling)
        return activeSysts

    def writeDatacard(self,systFile=''):
        ### this should be the file with .baseSystematics + .AllOtherSysts from the generic syst adder.
        ### for now just use the .baseSystematics file
        shapesFile = self.shapesfile+'.baseSystematics'
        systematicsFile = self.systFile if len(systFile)==0 else systFile
        print "reading shapes from ",shapesFile," and activating systematics found in: ",systematicsFile
        tf = ROOT.TFile.Open(shapesFile)
        ### to save memory, the value of this is filled with the histogram only for central ones
        ### the central histo is needed to create the systs "on-the-fly"
        procAndHistos = {}
        for e in tf.GetListOfKeys() :
            name=e.GetName()
            obj=e.ReadObj()
            syst = name.split('_')[-1]
            if not any(idir in syst for idir in ['Up','Down']):
                if re.match('((m|p)\d+|0)',syst): # this was probably for some obsolete systs
                    (proc,syst) = ('_'.join(name.split('_')[1:-2]),'_'.join(name.split('_')[-2:]))
                    histo = None
                else:
                    #print "hname = %s" % name
                    (proc,syst) = ('_'.join(name.split('_')[1:]),'central')
                    histo = obj.Clone() # this changes the name adding _clone! resetting name below
                    histo.SetDirectory(None)
                    histo.SetName(name)
            else:
                (proc,syst) = ('_'.join(name.split('_')[1:-1]),name.split('_')[-1])
                histo = None
            procAndHistos[(proc,syst)] = histo
        tf.Close()

        processes   = [proc for (proc,syst) in procAndHistos if (syst=='central' and proc!='data' and proc!='data_obs')]
        signalprocs = sorted([p for p in filter(lambda x: x.startswith('{boson}_{charge}'.format(boson=self.boson,charge=self.charge)),processes)])
        signalProcAndInd = [(p,-1*i) for i,p in enumerate(signalprocs)]
        otherprocs = [p for p in filter(lambda x: x not in signalprocs, processes)]
        otherProcAndInd = [(p,i+1) for i,p in enumerate(otherprocs)]
        allprocs = signalProcAndInd + otherProcAndInd

        ### this is needed to keep the sorted matrix of the syst/process
        procindices = {}
        for i,p in enumerate(allprocs):
            procindices[p[0]] = i

        ### build the systematics matrix matching the process name with the systematic name in the activeFitSysts.txt file
        appliedSystMatrix = {}
        systGroups = {}
        systEnv = self.activateSystMatrix(systematicsFile)

        ## here apply all the systematics from the systFit in combination 
        ## with the systematicMultipliers.root file
        procAndCentralHistos = dict(filter(lambda x: x[0][1] == 'central', procAndHistos.items()))
        newProcAndSysts = self.applySystematics(procAndCentralHistos, systEnv)

        procAndHistos.update(newProcAndSysts)
        ## done applying systematics from the multipliers file

        for (proc,syst),histo in procAndHistos.iteritems():
            ### skip the central histograms
            if syst=='central': continue
            ### the mW steps need to be implemented when we know how it is implemented ion combinetf
            if 'massShift' in syst: 
                pass
            ### for each syst, there is Up/Down, just write the line once
            if any([syst.endswith('Up'),syst.endswith('Down')]):
                systName = syst.replace('Up','').replace('Down','')
                if systName not in appliedSystMatrix:
                    appliedSystMatrix[systName] = ['-' for p in allprocs]
                for systKey,(sytPatt, procPatt, group, scaling) in systEnv.iteritems():
                    if sytPatt.match(systName) and procPatt.match(proc):
                        process_index = procindices[proc]
                        appliedSystMatrix[systName][process_index] = str(scaling)
                        if group not in systGroups: 
                            systGroups[group] = [systName]
                        else:                       
                            if systName not in systGroups[group]: 
                                systGroups[group].append(systName)

        if options.mergeRoot or options.remakeMultipliers:
            os.system('hadd -f {of} {of}.baseSystematics {of}.multiplierSystematics'.format(of=self.shapesfile))

        ### now write the datacard
        channel = '{boson}_{charge}'.format(boson='Wmunu' if self._options.wmass else 'Zmumu',charge=self.charge)
        datacard = open(self.cardfile,'w')
        ### -- header --
        datacard.write("imax 1\n")
        datacard.write("jmax *\n")
        datacard.write("kmax *\n")
        datacard.write('##----------------------------------\n') 
        datacard.write('shapes *  *  {shapes} x_$PROCESS x_$PROCESS_$SYSTEMATIC\n'.format(shapes=os.path.abspath(self.shapesfile)))
        datacard.write('##----------------------------------\n')
        datacard.write('bin    {channel}'.format(channel=channel)+'\n')
        datacard.write('observation    {obsyield}\n'.format(obsyield=procAndHistos[('data_obs','central')].Integral()))
        datacard.write('##----------------------------------\n')
        klen = 10
        kpatt = " %%%ds "  % klen
        datacard.write('bin                     %s' % ' '.join([kpatt % channel  for (proc,ind) in allprocs])+'\n')
        datacard.write('process                 %s' % ' '.join([kpatt % proc     for (proc,ind) in allprocs])+'\n')
        datacard.write('process                 %s' % ' '.join([kpatt % str(ind) for (proc,ind) in allprocs])+'\n')
        datacard.write('rate                    %s' % ' '.join([kpatt % '-1'     for (proc,ind) in allprocs])+'\n')
        datacard.write('##----------------------------------\n')

        ### -- single systematics --
        excludedSysts = []
        for systName,values in sorted(appliedSystMatrix.iteritems()):
            if isExcludedNuisance(self.excludeNuisances, systName, self.keepNuisances): 
                excludedSysts.append(systName)
                continue
            ### do not write systematic lines if it is not applied to any process
            if any([v!='-' for v in values]):
                datacard.write( '%-15s   shape %s\n' % (systName," ".join([kpatt % v for v in values])))

        ### -- groups --
        datacard.write('\n\n')
        for group,systs in systGroups.iteritems():
            sortedsysts = [x for x in systs if x not in excludedSysts]
            sortedsysts = sorted(sortedsysts)
            datacard.write( '%-15s   group = %s\n' % (group," ".join(sortedsysts)) )
        datacard.close()

        print "Done. Datacard in ",self.cardfile


def prepareChargeFit(options, charges=["plus"]):

    suffix = 'card' if options.freezePOIs else 'card_withXsecMask'
    if options.fitSingleCharge: 
    # this makes sense only when a single datacard for a given charge is different from those that would be used
    # for the combination (e.g. because to facilitate the combination some lines that require both charges are 
    # added to the datacard for a specific charge)
        suffix += "_singleCharge{ch}".format(ch=charges[0])

    datacards=[]; 
    channels=[]
    binname = "Wmunu" if options.wmass else "Zmumu"
    for charge in charges:
        datacards.append(os.path.abspath(options.inputdir)+"/{b}_{ch}_card.txt".format(b=binname,ch=charge))
        channels.append('{b}_{ch}'.format(b=binname,ch=charge))
        maskedChannels = ['InAcc']
        #if something_happens:
        #    maskedChannels.append('OutAcc')
        if not options.freezePOIs:
            for mc in maskedChannels:
                datacards.append(os.path.abspath(options.inputdir)+"/{b}_{ch}_xsec_{maskchan}_card.txt".format(b=binname,ch=charge,maskchan=mc))
                channels.append('{b}_{ch}_xsec_{maskchan}'.format(b=binname,ch=charge,maskchan=mc))

    print "="*20
    print "Looking for these cards"
    print "-"*20
    for d in datacards:
        print d
    print "="*20

    ### prepare the combineCards and txt2hdf5 commands
    if True: #sum([os.path.exists(card) for card in datacards])==len(datacards):
        if options.fitSingleCharge:
            print "I am going to run fit for single charge {ch}".format(ch=charges[0])
        else:
            print "Cards for W+ and W- done. Combining them now..."

        combinedCard = os.path.abspath(options.inputdir)+"/"+binname+'_'+suffix+'.txt'
        ccCmd = 'combineCards.py --noDirPrefix '+' '.join(['{ch}={dcfile}'.format(ch=channels[i],dcfile=card) for i,card in enumerate(datacards)])+' > '+combinedCard

        txt2hdf5Cmd = 'text2hdf5.py {cf} '.format(cf=combinedCard)

        ## here making the TF meta file
        if not options.freezePOIs:
            # actually it seems this would work also for frozen POIs, to be checked
            maskchan = [' --maskedChan {b}_{charge}_xsec_{maskchan}'.format(b=binname,charge=ch,maskchan=mc) for ch in charges for mc in maskedChannels]
            txt2hdf5Cmd += ' {maskch} --X-allow-no-background'.format(maskch=' '.join(maskchan))

        if not options.doSystematics:
            txt2hdf5Cmd = txt2hdf5Cmd + " -S 0 "
        if len(options.postfix):
            txt2hdf5Cmd = txt2hdf5Cmd + " --postfix " + options.postfix
        if options.clipSystVariations > 0.0:
            txt2hdf5Cmd = txt2hdf5Cmd + " --clipSystVariations " + str(options.clipSystVariations)
        if options.clipSystVariationsSignal > 0.0:
            txt2hdf5Cmd = txt2hdf5Cmd + " --clipSystVariationsSignal " + str(options.clipSystVariationsSignal)

        ## run the commands: need cmsenv in the combinetf release
        print 
        print ccCmd
        print 
        if not options.justFit:
            os.system(ccCmd)
        print "Combined card in ",combinedCard
        print
        print txt2hdf5Cmd
        print 
        if not options.skip_text2hdf5: 
            print "Running text2hdf5.py, it might take time ..."
            os.system(txt2hdf5Cmd)

        metafilename = combinedCard.replace('.txt','.hdf5')
        if len(options.postfix):
            metafilename = metafilename.replace('.hdf5','_%s.hdf5' % options.postfix)

        bbboptions = " --binByBinStat "
        if not options.noCorrelateXsecStat: bbboptions += "--correlateXsecStat "
        combineCmd = 'combinetf.py -t -1 {bbb} {metafile} --saveHists --computeHistErrors --doh5Output '.format(metafile=metafilename, bbb="" if options.noBBB else bbboptions)
        if options.freezePOIs:
            combineCmd += " --POIMode none"
            if options.doImpactsOnMW:
                combineCmd += " --doImpacts  "
        else:
            if options.doSystematics:
                combineCmd += " --doImpacts  "

        fitdir_data = "{od}/fit/data/".format(od=os.path.abspath(options.inputdir))
        fitdir_Asimov = "{od}/fit/hessian/".format(od=os.path.abspath(options.inputdir))
        if not os.path.exists(fitdir_data):
            print "Creating folder", fitdir_data
            os.system("mkdir -p " + fitdir_data)
        if not os.path.exists(fitdir_Asimov):
            print "Creating folder", fitdir_Asimov
            os.system("mkdir -p " + fitdir_Asimov)
        print ""
        fitPostfix = "" if not len(options.postfix) else ("_"+options.postfix)

        print "Use the following command to run combine (add --seed <seed> to specify the seed, if needed). See other options in combinetf.py"
        print
        combineCmd_data = combineCmd.replace("-t -1 ","-t 0 ")
        combineCmd_data = combineCmd_data + " --postfix Data{pf}_bbb{b} --outputDir {od} ".format(pf=fitPostfix, od=fitdir_data, b="0" if options.noBBB else "1_cxs0" if options.noCorrelateXsecStat else "1_cxs1")
        combineCmd_Asimov = combineCmd + " --postfix Asimov{pf}_bbb{b} --outputDir {od} ".format(pf=fitPostfix, od=fitdir_Asimov,  b="0" if options.noBBB else "1_cxs0" if options.noCorrelateXsecStat else "1_cxs1")
        print combineCmd_data
        if not options.skip_combinetf and not options.skipFitData:
            os.system(combineCmd_data)
        print ""
        print combineCmd_Asimov            
        if not options.skip_combinetf and not options.skipFitAsimov:
            os.system(combineCmd_Asimov)


def combineCharges(options):
    prepareChargeFit(options, charges=['plus','minus'])


def makeAllSystematicMultipliers(indir):
    generic_systs = [ {'name': 'fakesPt1', 'pt': [[20,70]], 'eta': [[-2.4,-2.1]], 'size': [0.95] },
                      {'name': 'fakesPt2', 'pt': [[20,70]], 'eta': [[-2.1,-1.5]], 'size': [0.98] },
                      {'name': 'fakesPt3', 'pt': [[20,70]], 'eta': [[-1.5, 1.5]], 'size': [0.85] },
                      {'name': 'fakesPt4', 'pt': [[20,70]], 'eta': [[ 1.5, 2.4]], 'size': [0.77] },
                      {'name': 'lnN30'   , 'pt': [[20,70]], 'eta': [[-2.5, 2.5]], 'size': [1.30] },
                      {'name': 'lumi'    , 'pt': [[20,70]], 'eta': [[-2.5, 2.5]], 'size': [1.025] } ]
    
    allHistos = []
    
    binningFile = indir+'/binningPtEta.txt'

    ## first make the ones for fakes etc. defined from the generic_systs dictionary
    ## ============================================================================
    
    for s in generic_systs:
        allHistos += [ i for i in ms.makeGenericMultipliers(binningFile, s['name'], s['pt'], s['eta'], s['size'])]

    ## then make the ones for the effStat parameters
    ## =============================================
    ## this needs files systEff_trgmu_{plus|minus}_mu.root in the directory given

    effStatHistos = ms.makeEffStatMultipliers(binningFile, 
                    '/afs/cern.ch/work/m/mdunser/public/cmssw/w-mass-13TeV/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/../postprocessing/data/leptonSF/new2016_madeSummer2018/')

    allHistos += effStatHistos
    
    ## write it out into file called systematicMultipliers.root
    outfile = ROOT.TFile(indir+'/systematicMultipliers.root', 'RECREATE')
    for h in allHistos:
        h.Write()
    outfile.Close()

## marc applySystematics('data_fakes', '.*fakesPt2|.*lnN30.*|.*EffStat1.*', 'testEffSyst.root', 'appliedSysts.root')


if __name__ == "__main__":    

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [self._options]')
    parser.add_option('-m','--merge-root', dest='mergeRoot', default=False, action='store_true', help='Merge the root files with the inputs also')
    parser.add_option('-i','--input', dest='inputdir', default='', type='string', help='input directory with all the cards inside')
    parser.add_option('-f','--flavor', dest='flavor', default='mu', type='string', help='lepton flavor (mu,el)')
    parser.add_option('-C','--charge', dest='charge', default='plus,minus', type='string', help='process given charge. default is both')
    parser.add_option(     '--pdf-shape-only'   , dest='pdfShapeOnly' , default=False, action='store_true', help='Normalize the mirroring of the pdfs to central rate.')
    parser.add_option(     '--remake-syst-multipliers'   , dest='remakeMultipliers' , default=False, action='store_true', help='Remake all the systematic multipliers.')
    parser.add_option(     '--comb'   , dest='combineCharges' , default=False, action='store_true', help='Combine W+ and W-, if single cards are done')
    parser.add_option(     '--postfix',    dest='postfix', type="string", default="", help="Postfix for .hdf5 file created with text2hdf5.py when combining charges");
    parser.add_option('--wmass', dest='wmass', action="store_true", default=False, help="Make cards for the wmass analysis. Default is wlike");
    # options for card maker and fit
    parser.add_option("-S",  "--doSystematics", type=int, default=1, help="enable systematics when running text2hdf5.py (-S 0 to disable them)")
    parser.add_option(       "--exclude-nuisances", dest="excludeNuisances", default="", type="string", help="Pass comma-separated list of regular expressions to exclude some systematics")
    parser.add_option(       "--keep-nuisances", dest="keepNuisances", default="", type="string", help="Pass comma-separated list of regular expressions to keep some systematics, overriding --exclude-nuisances. Can be used to keep only 1 syst while excluding all the others")
    parser.add_option("", "--clipSystVariations", type=float, default=-1.,  help="Clipping of syst variations, passed to text2hdf5.py")
    parser.add_option("", "--clipSystVariationsSignal", type=float, default=-1.,  help="Clipping of signal syst variations, passed to text2hdf5.py")
    parser.add_option('--fp','--freezePOIs'  , dest='freezePOIs'   , default=False, action='store_true', help='run tensorflow with --freezePOIs (for the pdf only fit)')
    parser.add_option(       '--no-bbb'  , dest='noBBB', default=False, action='store_true', help='Do not use bin-by-bin uncertainties')
    parser.add_option(       '--no-correlate-xsec-stat'  , dest='noCorrelateXsecStat', default=False, action='store_true', help='Do not use option --correlateXsecStat when using bin-by-bin uncertainties ')
    parser.add_option(       '--no-text2hdf5'  , dest='skip_text2hdf5', default=False, action='store_true', help='skip running text2hdf5.py at the end, only prints command (useful if hdf5 file already exists, or for tests)')
    parser.add_option(       '--no-combinetf'  , dest='skip_combinetf', default=False, action='store_true', help='skip running combinetf.py at the end, just print command (useful for tests)')
    parser.add_option(       '--just-fit', dest='justFit' , default=False, action='store_true', help='Go directly to fit part (can also skip text2hdf5 with --no-text2hdf5)')
    parser.add_option(       '--skip-fit-data', dest='skipFitData' , default=False, action='store_true', help='If True, fit only Asimov')
    parser.add_option(       '--skip-fit-asimov', dest='skipFitAsimov' , default=False, action='store_true', help='If True, fit only data')
    parser.add_option(      '--fit-single-charge', dest='fitSingleCharge', default=False, action='store_true', help='Prepare datacard for single-charge fit (groups with other charge are skipped')

    # --impacts-mW might not be needed or might not work, to be checked
    parser.add_option(      '--impacts-mW', dest='doImpactsOnMW', default=False, action='store_true', help='Set up cards to make impacts of nuisances on mW')
    (options, args) = parser.parse_args()
    
    charges = options.charge.split(',')

    if options.combineCharges:
        if options.fitSingleCharge:
            print "Error: options --fit-single-charge and --comb are incompatible. Abort"
            quit()
        if len(charges) != 2:
            print "Error: --comb requires two charges, use -C 'plus,minus' and try again"
            quit()

    if options.remakeMultipliers:
        print 'making all the systematic multipliers'
        makeAllSystematicMultipliers(options.inputdir)
        print 'done making all the systematic multipliers'
    

    for charge in charges:
        if not options.justFit:
            cm = CardMaker(options,charge,options.flavor)
            if options.mergeRoot:
                cm.mergeChunks()
            cm.writeDatacard()
        if options.fitSingleCharge:
            prepareChargeFit(options, charges=[charge])
            print "-"*30
            print "Done fitting charge {ch}".format(ch=charge)
            print "="*30

    if options.combineCharges and len(charges)==2:
        combineCharges(options)                
        print "-"*30
        print "Done fitting charge combination"
        print "="*30
