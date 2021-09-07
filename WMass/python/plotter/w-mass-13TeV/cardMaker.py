#!/bin/env python

import os, re, copy, math, array
#import root_numpy
import makeSystematicMultipliers as ms
import utilities
import argparse

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

utilities = utilities.util()

def safeSystem(cmd, dryRun=False, quitOnFail=True):
    print(cmd)
    if not dryRun:
        res = os.system(cmd)
        if res:
            print('-'*30)
            print("safeSystem(): error occurred when executing the following command. Aborting")
            print(cmd)
            print('-'*30)
            if quitOnFail:
                quit()
        return res
    else:
        return 0
        
def sortSystsForDatacard(params):

    params = sorted(params, key= lambda x: int(x.replace('pdf','')) if 'pdf' in x else 101 if 'alphaS' in x else 0)
    params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if ('muRmuF' in x and x != "muRmuF")  else 0)
    params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if (''.join([j for j in x if not j.isdigit()]) == 'muR' and x != "muR") else 0)
    params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if (''.join([j for j in x if not j.isdigit()]) == 'muF' and x != "muF") else 0)
    params = sorted(params, key= lambda x: utilities.getNFromString(x) if 'EffStat' in x else 0)
    return params
    
def isExcludedNuisance(excludeNuisances=[], name="", keepNuisances=[]):
    if len(excludeNuisances) and any(re.match(x,name) for x in excludeNuisances):
        if len(keepNuisances) and any(re.match(x,name) for x in keepNuisances):
            return False
        else:
            print(">>>>> Excluding nuisance: ", name)
            return True
    else:
        return False

class CardMaker:
    def __init__(self,options,charge):
        self._options = options
        self.charge = charge
        self.flavor = self._options.flavor
        #self.data_eras = ["B", "C", "D", "E", "F", "F_postVFP", "G", "H"] # not needed I think, FIXME
        self.excludeNuisances = []
        self.keepNuisances = []
        if len(options.excludeNuisances):
            self.excludeNuisances = options.excludeNuisances.split(",")
        if len(options.keepNuisances):
            self.keepNuisances = options.keepNuisances.split(",")
                
        pathToImport = os.path.abspath("./w-mass-13TeV/") + "/" # should work
        if self._options.isWlike:
            pathToImport += 'wlike_{fl}'.format(fl=self.flavor)
            self.boson = 'Z{fl}{fl}'.format(fl=self.flavor)
            self.otherboson = 'W{fl}nu'.format(fl=self.flavor)
        else:
            pathToImport += 'wmass_{fl}'.format(fl=self.flavor)
            self.boson = 'W{fl}nu'.format(fl=self.flavor)
            self.otherboson = 'Z{fl}{fl}'.format(fl=self.flavor)
        #self.systematics  = self.getSystList()
        #self.centralfiles = self.getCentralProcesses()
        shapeFilePostfix = ""
        if self._options.shapeFilePostfix:
            shapeFilePostfix = "_" + self._options.shapeFilePostfix
        self.shapesfile = os.path.join(self._options.inputdir,self.boson+'_{ch}_shapes{sfp}.root'.format(ch=charge, sfp=shapeFilePostfix))
        self.cardfile   = os.path.join(self._options.inputdir+self._options.cardFolder,self.boson+'_{ch}_card.txt'   .format(ch=charge))
        self.systFile = pathToImport+'/systsFit.txt'
        print("Copying syst configuration file")
        copyCmd = "cp %s %s" % (self.systFile, self._options.inputdir+self._options.cardFolder)
        #print(copyCmd)
        safeSystem(copyCmd, dryRun=self._options.dryRun)

        # self.centralHistograms = self.getCentralHistograms() # to implement


    # def getCentralProcesses(self):        
    #     #acoeffs = ['ac']+['a{ic}'.format(ic=i) for i in range(8)]
    #     acoeffs = ['']
    #     # signal
    #     sig_proc          = ['{boson}_{charge}{ic}'.format(charge=self.charge,ic=("_"+coeff) if len(coeff) else "",boson=self.boson) for coeff in acoeffs]
    #     # antisignal
    #     other_boson_procs = ['{otherboson}_{charge}'.format(charge=self.charge,otherboson=self.otherboson)]
    #     # other stuff (for data use then global one (not split by era, it is created afterwards)
    #     other_files = ["Wtaunu","Ztautau","otherBkgHisto", "dataHisto"]
    #     otherprocs = ["{op}_{charge}".format(op=otherproc,charge=self.charge) for otherproc in other_files] 
    #     return sig_proc + other_boson_procs + otherprocs

    # def getSystList(self):
    #     ### systematics that are applied to both W and Z, for mu and tau decays
    #     NPDFSYSTS=2
    #     baseSysts  = ['pdf{i}'.format(i=ipdf) for ipdf in range(1,1+NPDFSYSTS)]
    #     baseSysts += ['alphaS{idir}'.format(idir=idir) for idir in ['Up','Dn']]
    #     # FIXME: need a better way to decide if and how these are binned, depending on what was done in make_wmass_cards.py
    #     qcdSystsUnbin  = ['{qcdpar}{idir}'.format(qcdpar=par,idir=idir) for par in ['muR','muF','muRmuF'] for idir in ['Up','Dn']]     
    #     NVPTBINS = 2
    #     qcdSystsPtChargebin  = ['{qcdpar}{ipt}{ch}{idir}'.format(qcdpar=par,ipt=ipt,ch=self.charge,idir=idir) for par in ['muR','muF','muRmuF'] for ipt in range(1,1+NVPTBINS) for idir in ['Up','Dn']]     
    #     #baseSysts += ['kalPtErr{i}{idir}'.format(i=istat,idir=idir) for istat in range(133) for idir in ['Up','Dn']]
    #     #baseSysts += ['kalPtClosureErr{idir}'.format(idir=idir) for idir in ['Up','Dn']]
    #     ### W/Z mass points
    #     MASSVARIATIONS = [10* i for i in range(1,3)]
    #     massPoints = ['massShift{v}MeV{d}".format(v=mvar,d=idir)' for mvar in MASSVARIATIONS for idir in ['Up','Down']]
    #     ### other systematics
    #     # otherSysts = ['fsr']
    #     otherSysts = []
    #     ### build the systematic dictionary for the root files
    #     systsCards = {self.boson: baseSysts + qcdSystsPtChargebin + massPoints + otherSysts,
    #                   self.otherboson: baseSysts + qcdSystsUnbin,
    #                   'Wtaunu': baseSysts + qcdSystsPtChargebin + massPoints,
    #                   'Ztautau': baseSysts + qcdSystsUnbin}
    #     return systsCards


    # this needs to be revised
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
                            #print("hname = %s     systname = %s" % (proc_hist.GetName(),systName))
                            tmp_var  = copy.deepcopy(proc_hist.Clone(proc_hist.GetName()+'_'+systName))
                            tmp_var.Multiply(tmp_syst)
                            allVars.append(tmp_var)
                            dictProcSyst[(proc, systName)] = None
    
        outfile = ROOT.TFile(options.inputdir+'/{f}_{ch}_shapes.root.multiplierSystematics'.format(f=self.boson,ch=self.charge), 'RECREATE')
        
        for v in allVars:
            v.Write()
        outfile.Close()

        return dictProcSyst

    
    def activateSystMatrix(self, systFile):
        activeSysts = {}
        for line in open(systFile, 'r'):
            if re.match("\s*#.*", line): continue
            line = re.sub("#.*","",line).strip()
            if len(line) == 0: continue
            field = [f.strip() for f in line.split(':')]
            if len(field) < 5:
                raise RuntimeError("Malformed line %s in file %s"%(line.strip(),systFile))
            elif len(field) == 5:
                (name, sytPatt, procPatt, group, sysType, scaling) = field[:5] + [1]
            else:
                (name, sytPatt, procPatt, group, sysType, scaling) = field[:6]
            activeSysts[name] = (re.compile(sytPatt+"$"), re.compile(procPatt+"$"), group, sysType, scaling)
        return activeSysts

    def writeDatacard(self, systFile=''):
        ### this should be the file with .baseSystematics + .AllOtherSysts from the generic syst adder.
        ### for now just use the .baseSystematics file
        #shapesFile = self.shapesfile+'.baseSystematics'
        shapesFile = self.shapesfile
        systematicsFile = self.systFile if len(systFile)==0 else systFile
        print("reading shapes from ",shapesFile," and activating systematics found in: ",systematicsFile)
        tf = ROOT.TFile.Open(shapesFile)
        ### to save memory, the value of this is filled with the histogram only for central ones
        ### the central histo is needed to create the systs "on-the-fly"
        procAndHistos = {}
        for e in tf.GetListOfKeys() :
            name=e.GetName()
            obj=e.ReadObj()
            syst = name.split('_')[-1]
            if not any(idir in syst for idir in ['Up','Down']):
                if re.match('((m|p)\d+|0)',syst): # this was probably for some obsolete systs like pt-scale systs decorrelated between 2 pt bins
                    (proc,syst) = ('_'.join(name.split('_')[1:-2]),'_'.join(name.split('_')[-2:]))
                    histo = None
                else:
                    #print("hname = %s" % name)
                    (proc,syst) = ('_'.join(name.split('_')[1:]),'central')
                    histo = obj.Clone() # this changes the name adding _clone! resetting name below
                    histo.SetDirectory(0)
                    histo.SetName(name)
            else:
                (proc,syst) = ('_'.join(name.split('_')[1:-1]), name.split('_')[-1])
                histo = None
            procAndHistos[(proc,syst)] = histo
        tf.Close()

        processes   = [proc for (proc,syst) in procAndHistos if (syst == 'central' and proc != 'data' and proc != "data_obs" and not proc.startswith("pseudodata") and proc != self._options.dataname)]
        if self._options.allProcessesBackground:
            allprocs = [(p,i+1) for i,p in enumerate(processes)]
        else:
            # for W, might treat both Wmunu_plus and Wmunu_minus as signal
            # as it is now, the card combination will fail, because Wmunu_plus would be signal in one card and background in the other
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
        systTypes = {} # used to decide which kind of systematic we have (shapes, lnN, shapesNoConstraint, ...)
        systGroups = {}
        systEnv = self.activateSystMatrix(systematicsFile)

        ## here apply all the systematics from the systFit in combination 
        ## with the systematicMultipliers.root file
        procAndCentralHistos = dict(filter(lambda x: x[0][1] == 'central', procAndHistos.items()))

        ## skip for now, might not need
        #newProcAndSysts = self.applySystematics(procAndCentralHistos, systEnv)
        #procAndHistos.update(newProcAndSysts)
        ## done applying systematics from the multipliers file
        
        # do all non-lnN nuisances, for which histograms exist
        # this will do nothing on lnN nuisances, which have no histogram
        # if some systs have to be decorrelated by charge, the histograms shoudl already exist with appropriate name
        # and the syst should already inherit the name with charge inside
        for (proc,syst),histo in iter(procAndHistos.items()):
            ### skip the central histograms
            if syst=='central': continue
            if any(syst.endswith(x) for x in ["Up", "Down"]):
                systName = syst.replace('Up','').replace('Down','')
                if systName not in appliedSystMatrix:
                    appliedSystMatrix[systName] = ['-' for p in allprocs]
                for systKey,(systPatt, procPatt, group, sysType, scaling) in iter(systEnv.items()):
                    if systPatt.match(systName) and procPatt.match(proc):
                        process_index = procindices[proc]
                        appliedSystMatrix[systName][process_index] = str(scaling)
                        systTypes[systName] = sysType
                        if group != "SKIPGROUP":
                            if group not in systGroups: 
                                systGroups[group] = [systName]
                            else:                       
                                if systName not in systGroups[group]: 
                                    systGroups[group].append(systName)

        # distinguish lnN nuisances from non-lnN (for which histograms exist)
        for systKey,(systPatt, procPatt, group, sysType, scaling) in iter(systEnv.items()):
            if sysType != "lnN":
                continue
            systName = systKey
            if isExcludedNuisance(self.excludeNuisances, systName, self.keepNuisances):
                continue
            if systName not in appliedSystMatrix:
                appliedSystMatrix[systName] = ['-' for p in allprocs]
            for proc,index in iter(procindices.items()):
                if procPatt.match(proc):
                    appliedSystMatrix[systName][index] = str(scaling)
            systTypes[systName] = sysType
            if group not in systGroups: 
                systGroups[group] = [systName]
            else:                       
                if systName not in systGroups[group]: 
                    systGroups[group].append(systName)
                                
        #if options.mergeRoot or options.remakeMultipliers:
        if self._options.remakeMultipliers:
            safeSystem('hadd -f {of} {of}.baseSystematics {of}.multiplierSystematics'.format(of=self.shapesfile), dryRun=self._options.dryRun)

        ### now write the datacard
        channel = '{boson}_{charge}'.format(boson=self.boson,charge=self.charge)
        datacard = open(self.cardfile,'w')
        ### -- header --
        datacard.write("imax 1\n")
        datacard.write("jmax *\n")
        datacard.write("kmax *\n")
        datacard.write('##----------------------------------\n') 
        datacard.write('shapes *  *  {shapes} x_$PROCESS x_$PROCESS_$SYSTEMATIC\n'.format(shapes=os.path.abspath(self.shapesfile)))
        datacard.write('##----------------------------------\n')
        datacard.write('bin    {channel}'.format(channel=channel)+'\n')
        datacard.write('observation    {obsyield}\n'.format(obsyield=procAndHistos[(self._options.dataname,'central')].Integral()))
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
        for systName in sortSystsForDatacard(list(appliedSystMatrix.keys())):
            values = appliedSystMatrix[systName]
            if isExcludedNuisance(self.excludeNuisances, systName, self.keepNuisances): 
                excludedSysts.append(systName)
                continue
            ### do not write systematic lines if it is not applied to any process
            if any([v != '-' for v in values]):
                datacard.write( '%-15s %-10s  %s\n' % (systName, systTypes[systName], " ".join([kpatt % v for v in values])))

        ### -- groups --
        datacard.write('\n\n')
        for group,systs in iter(systGroups.items()):
            sortedsysts = [x for x in systs if x not in excludedSysts]
            if len(sortedsysts) == 0:
                continue
            sortedsysts = sortSystsForDatacard(sortedsysts)
            datacard.write( '%-15s   group = %s\n' % (group," ".join(sortedsysts)) )

        # now noiGroup on mass nuisance parameter
        datacard.write( '%-15s   noiGroup = %s\n' % ("massnoi", self._options.massNuisance))

        datacard.close()
        print("Done. Datacard in ",self.cardfile)


def prepareChargeFit(options, charges=["plus"]):

    cardSubfolderFullName = options.inputdir + options.cardFolder
    postfix = options.postfix
    cardkeyname = 'card' if options.freezePOIs else 'card_withXsecMask'
    if options.fitSingleCharge: 
    # this cardkeyname is needed only when a single datacard for a given charge is different from those that would be used
    # for the combination (e.g. because to facilitate the combination some lines that require both charges are 
    # added to the datacard for a specific charge)
        cardkeyname += "_singleCharge{ch}".format(ch=charges[0])
        if postfix == "":
            postfix = "only{ch}".format(ch=charges[0])
        else:
            postfix = postfix + "_only{ch}".format(ch=charges[0])
        
    datacards=[]; 
    channels=[]
    binname = "Z{fl}{fl}".format(fl=options.flavor) if options.isWlike else "W{fl}nu".format(fl=options.flavor)

    for charge in charges:
        datacards.append(os.path.abspath(cardSubfolderFullName)+"/{b}_{ch}_card.txt".format(b=binname,ch=charge))
        channels.append('{b}_{ch}'.format(b=binname,ch=charge))
        maskedChannels = ['InAcc']
        #if something_happens:
        #    maskedChannels.append('OutAcc')
        if not options.freezePOIs:
            for mc in maskedChannels:
                datacards.append(os.path.abspath(cardSubfolderFullName)+"/{b}_{ch}_xsec_{maskchan}_card.txt".format(b=binname,ch=charge,maskchan=mc))
                channels.append('{b}_{ch}_xsec_{maskchan}'.format(b=binname,ch=charge,maskchan=mc))

    print('='*30)
    print("Looking for these cards")
    print('-'*30)
    for d in datacards:
        print(d)
    print('='*30)

    ### prepare the combineCards and txt2hdf5 commands
    if sum([os.path.exists(card) for card in datacards]) == len(datacards):
        if options.fitSingleCharge:
            print("I am going to run fit for single charge {ch}".format(ch=charges[0]))
        else:
            print("Cards for W+ and W- done. Combining them now...")

        combinedCard = "{d}/{b}_{s}.txt".format(d=os.path.abspath(cardSubfolderFullName), b=binname, s=cardkeyname)
        ccCmd = "combineCards.py --noDirPrefix {cards} > {combinedCard} ".format(cards=' '.join(['{ch}={dcfile}'.format(ch=channels[i],dcfile=card) for i,card in enumerate(datacards)]), combinedCard=combinedCard)
        ## run the commands: need cmsenv in the combinetf release
        print 
        print 
        if options.justFit:
            print(ccCmd)
        else:
            safeSystem(ccCmd, dryRun=options.dryRun)
        print("Combined card in ",combinedCard)
        print

        txt2hdf5Cmd = 'text2hdf5.py {cf} --dataset {dn}'.format(cf=combinedCard, dn=options.dataname)

        ## here making the TF meta file
        if not options.freezePOIs:
            # actually it seems this would work also for frozen POIs, to be checked
            maskchan = [' --maskedChan {b}_{charge}_xsec_{maskchan}'.format(b=binname,charge=ch,maskchan=mc) for ch in charges for mc in maskedChannels]
            txt2hdf5Cmd += ' {maskch} --X-allow-no-background'.format(maskch=' '.join(maskchan))

        if options.allProcessesBackground:
            txt2hdf5Cmd += " --X-allow-no-signal"
            
        if not options.doSystematics:
            txt2hdf5Cmd = txt2hdf5Cmd + " -S 0 "
        if len(postfix):
            txt2hdf5Cmd = txt2hdf5Cmd + " --postfix " + postfix
        if options.clipSystVariations > 0.0:
            txt2hdf5Cmd = txt2hdf5Cmd + " --clipSystVariations " + str(options.clipSystVariations)
        if options.clipSystVariationsSignal > 0.0:
            txt2hdf5Cmd = txt2hdf5Cmd + " --clipSystVariationsSignal " + str(options.clipSystVariationsSignal)

        print 
        if options.skip_text2hdf5: 
            print(txt2hdf5Cmd)
        else:
            print("Running text2hdf5.py, it might take time ...")
            safeSystem(txt2hdf5Cmd, dryRun=options.dryRun)
            
        metafilename = combinedCard.replace('.txt','.hdf5')
        if len(postfix):
            metafilename = metafilename.replace('.hdf5','_%s.hdf5' % postfix)

        bbboptions = " --binByBinStat "
        if not options.noCorrelateXsecStat: bbboptions += "--correlateXsecStat "
        combineCmd = 'combinetf.py -t -1 {bbb} {metafile} --saveHists --computeHistErrors --doh5Output '.format(metafile=metafilename, bbb="" if options.noBBB else bbboptions)
        if options.combinetfOption:
            combineCmd += " %s" % options.combinetfOption
        if options.freezePOIs:
            combineCmd += " --POIMode none"
            if options.doImpactsOnMW:
                combineCmd += " --doImpacts  "
        else:
            if options.doSystematics:
                combineCmd += " --doImpacts  "

        fitdir_data = "{od}/fit/data/".format(od=os.path.abspath(cardSubfolderFullName))
        fitdir_Asimov = "{od}/fit/hessian/".format(od=os.path.abspath(cardSubfolderFullName))
        if not os.path.exists(fitdir_data):
            print("Creating folder", fitdir_data)
            safeSystem("mkdir -p " + fitdir_data, dryRun=options.dryRun)
        if not os.path.exists(fitdir_Asimov):
            print("Creating folder", fitdir_Asimov)
            safeSystem("mkdir -p " + fitdir_Asimov, dryRun=options.dryRun)
        print("")
        fitPostfix = "" if not len(postfix) else ("_"+postfix)

        print("Use the following command to run combine (add --seed <seed> to specify the seed, if needed). See other options in combinetf.py")
        print
        combineCmd_data = combineCmd.replace("-t -1 ","-t 0 ")
        combineCmd_data = combineCmd_data + " --postfix Data{pf}_bbb{b} --outputDir {od} ".format(pf=fitPostfix, od=fitdir_data, b="0" if options.noBBB else "1_cxs0" if options.noCorrelateXsecStat else "1_cxs1")
        combineCmd_Asimov = combineCmd + " --postfix Asimov{pf}_bbb{b} --outputDir {od} ".format(pf=fitPostfix, od=fitdir_Asimov,  b="0" if options.noBBB else "1_cxs0" if options.noCorrelateXsecStat else "1_cxs1")
        if not options.skip_combinetf and not options.skipFitData:
            safeSystem(combineCmd_data, dryRun=options.dryRun)
        else:
            print(combineCmd_data)
        print
        if not options.skip_combinetf and not options.skipFitAsimov:
            safeSystem(combineCmd_Asimov, dryRun=options.dryRun)
        else:
            print(combineCmd_Asimov)

            
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

    # effStatHistos = ms.makeEffStatMultipliers(binningFile, 
    #                 '/afs/cern.ch/work/m/mdunser/public/cmssw/w-mass-13TeV/CMSSW_9_4_12/src/CMGTools/WMass/python/plotter/../postprocessing/data/leptonSF/new2016_madeSummer2018/')

    # allHistos += effStatHistos
    
    ## write it out into file called systematicMultipliers.root
    outfile = ROOT.TFile(indir+'/systematicMultipliers.root', 'RECREATE')
    for h in allHistos:
        h.Write()
    outfile.Close()

## marc applySystematics('data_fakes', '.*fakesPt2|.*lnN30.*|.*EffStat1.*', 'testEffSyst.root', 'appliedSysts.root')


if __name__ == "__main__":    

    parser = argparse.ArgumentParser()
    #parser.add_argument('-m','--merge-root', dest='mergeRoot', default=False, action='store_true', help='Merge the root files with the inputs also')
    parser.add_argument('-i','--input', dest='inputdir', default='', type=str, help='input directory with all the cards inside')
    parser.add_argument('--cf','--card-folder', dest='cardFolder', default='nominal', type=str, help='Subfolder created inside inputdir to store all cards and fit results in a given configuration, for better bookkeeping when doing tests')
    parser.add_argument('-f','--flavor', dest='flavor', default='mu', type=str, choices=["mu","el"], help='lepton flavor (mu,el)')
    parser.add_argument('-c','--charge', dest='charge', default='plus,minus', type=str, help='process given charge. default is both')
    parser.add_argument(     '--remake-syst-multipliers'   , dest='remakeMultipliers' , default=False, action='store_true', help='Remake all the systematic multipliers.')
    parser.add_argument(     '--comb'   , dest='combineCharges' , default=False, action='store_true', help='Combine W+ and W-, if single cards are done')
    parser.add_argument(     '--postfix',    dest='postfix', type=str, default="", help="Postfix for .hdf5 file created with text2hdf5.py when combining charges");
    parser.add_argument(     '--shape-file-postfix',    dest='shapeFilePostfix', type=str, default="", help="Postfix to identify root file with histograms (e.g. Wmunu_plus_shapes_{postfix}.root). Useful to select different root files or for tests");
    parser.add_argument('--wlike', dest='isWlike', action="store_true", default=False, help="Make cards for the W-like analysis. Default is Wmass");
    # options for card maker and fit
    parser.add_argument("-S",  "--doSystematics", type=int, default=1, help="enable systematics when running text2hdf5.py (-S 0 to disable them)")
    parser.add_argument("-x",  "--exclude-nuisances", dest="excludeNuisances", default="", type=str, help="Pass comma-separated list of regular expressions to exclude some systematics")
    parser.add_argument("-k",  "--keep-nuisances", dest="keepNuisances", default="", type=str, help="Pass comma-separated list of regular expressions to keep some systematics, overriding --exclude-nuisances. Can be used to keep only 1 syst while excluding all the others")
    parser.add_argument(       "--clipSystVariations", type=float, default=-1.,  help="Clipping of syst variations, passed to text2hdf5.py")
    parser.add_argument(       "--clipSystVariationsSignal", type=float, default=-1.,  help="Clipping of signal syst variations, passed to text2hdf5.py")
    parser.add_argument('--fp','--freezePOIs'  , dest='freezePOIs'   , default=False, action='store_true', help='run tensorflow with --freezePOIs')
    parser.add_argument(       '--no-bbb'  , dest='noBBB', default=False, action='store_true', help='Do not use bin-by-bin uncertainties')
    parser.add_argument(       '--no-correlate-xsec-stat'  , dest='noCorrelateXsecStat', default=False, action='store_true', help='Do not use option --correlateXsecStat when using bin-by-bin uncertainties ')
    parser.add_argument('-d',  '--dry-run'  , dest='dryRun', default=False, action='store_true', help='Do not execute command to make cards or fit')
    parser.add_argument(       '--no-text2hdf5'  , dest='skip_text2hdf5', default=False, action='store_true', help='skip running text2hdf5.py at the end, only prints command (useful if hdf5 file already exists, or for tests)')
    parser.add_argument(       '--no-combinetf'  , dest='skip_combinetf', default=False, action='store_true', help='skip running combinetf.py at the end, just print command (useful for tests)')
    parser.add_argument(       '--just-fit', dest='justFit' , default=False, action='store_true', help='Go directly to fit part (can also skip text2hdf5 with --no-text2hdf5)')
    parser.add_argument(       '--skip-fit-data', dest='skipFitData' , default=False, action='store_true', help='If True, fit only Asimov')
    parser.add_argument(       '--skip-fit-asimov', dest='skipFitAsimov' , default=False, action='store_true', help='If True, fit only data')
    parser.add_argument(      '--fit-single-charge', dest='fitSingleCharge', default=False, action='store_true', help='Prepare datacard for single-charge fit. For each charge, a postfix is appended to option --postfix, so no need to add the charge explicitly')
    # --impacts-mW might not be needed or might not work, to be checked
    parser.add_argument(      '--impacts-mW', dest='doImpactsOnMW', default=False, action='store_true', help='Set up cards to make impacts of nuisances on mW')
    parser.add_argument(     '--mass-nuis',    dest='massNuisance', type=str, default="massShift50MeV", help="Use only this mass nuisance in the datacard and to define noiGroup. It overrides exclusion from option --exclude-nuisances");     parser.add_argument(      '--all-proc-background', dest='allProcessesBackground', default=False, action='store_true', help='Set all processes as backgrounds (e.g. when running the fit with fixed poi)')
    parser.add_argument("-D", "--dataset",  dest="dataname", default="data_obs",  type=str,  help="Name of the observed dataset (pass name without x_ in the beginning). Useful to fit another pseudodata histogram")
    parser.add_argument("--combinetf-option",  dest="combinetfOption", default="",  type=str,  help="Pass other options to combinetf")
    args = parser.parse_args()

    if not args.dryRun:
        try:
            cmssw = os.environ['CMSSW_BASE']
        except:
            cmssw = ""
        if cmssw == "":
            print("\n")
            print("Error: to use combinetf you need to activate cmsenv from a release.")
            print("You should work from a cmssw-cc7 singularity environment to get the release.")
            print("Aborting ...")
            print("\n")
            quit()

    if not args.inputdir.endswith("/"):
        args.inputdir += "/"
            
    if args.cardFolder:
        testName = args.cardFolder.rstrip("/")
        if testName in ["plus", "minus"]:
            # these names are already used, they contain the configurations used to make histograms
            print("Warning: folder specified with --card-folder cannot be named as plus|minus/ (you used %s)" % args.cardFolder)
            quit()
        if not args.cardFolder.endswith("/"):
            args.cardFolder += "/"
        cardFolderFullName = args.inputdir + args.cardFolder
        if not os.path.exists(cardFolderFullName):
            print("Creating folder", cardFolderFullName)
            safeSystem("mkdir -p " + cardFolderFullName, dryRun=False) # always create this folder, even for tests
            fcmd = open(cardFolderFullName+"cardMaker_command.txt", "w")
            fcmd.write("%s\n\n" % " ".join(sys.argv))
            fcmd.close()

    options = args

    charges = options.charge.split(',')

    if options.combineCharges:
        if options.fitSingleCharge:
            print("Error: options --fit-single-charge and --comb are incompatible. Abort")
            quit()
        if len(charges) != 2:
            print("Error: --comb requires two charges, use -C 'plus,minus' and try again")
            quit()
            
    if args.massNuisance:
        if not args.excludeNuisances:
            args.excludeNuisances = "massShift.*"
        else:
            args.excludeNuisances = "{orig},massShift.*".format(orig=args.excludeNuisances)
        if not args.keepNuisances:
            args.keepNuisances = args.massNuisance
        elif not re.match(args.massNuisance, args.keepNuisances):
            args.keepNuisances = "{orig},{m}".format(orig=args.keepNuisances, m=args.massNuisance)
        
    # skip for now, to be tested
    if options.remakeMultipliers:
        print('making all the systematic multipliers')
        makeAllSystematicMultipliers(options.inputdir)
        print('done making all the systematic multipliers')
    
    for charge in charges:
        if not options.justFit:
            cm = CardMaker(options, charge)
            #if options.mergeRoot:
            #    cm.mergeChunks()
            cm.writeDatacard()
        if options.fitSingleCharge:
            prepareChargeFit(options, charges=[charge])
            print('-'*30)
            print("Done fitting charge {ch}".format(ch=charge))
            print('-'*30)

    if options.combineCharges and len(charges)==2:
        combineCharges(options)                
        print('-'*30)
        print("Done fitting charge combination")
        print('-'*30)
