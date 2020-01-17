#!/bin/env python

import ROOT, sys,os,re,json, copy, math, root_numpy, array
from rollingFunctions import roll1Dto2D, dressed2D, unroll2Dto1D
import makeSystematicMultipliers as ms
import utilities
utilities = utilities.util()

class CardMaker:
    def __init__(self,options,charge,flavor):
        self._options = options
        self.charge = charge
        self.flavor = flavor

        ### import the right dictionary process - systematic list depending on wmass/wlike
        #pathToImport = os.environ['CMSSW_BASE']+'/src/CMGTools/WMass/python/plotter/w-mass-13TeV/'
        pathToImport = os.path.dirname(sys.argv[0])+'/'
        if self._options.wmass:
            pathToImport += 'wmass_'+flavor
            self.boson = 'W'
        else:
            pathToImport += 'wlike_'+flavor
            self.boson = 'Z'
        self.systematics  = self.getSystList()
        self.centralfiles = self.getCentralProcesses()
        self.shapesfile = os.path.join(self._options.inputdir,self.flavor+'_{ch}_shapes.root'.format(ch=charge))
        self.cardfile   = os.path.join(self._options.inputdir,self.flavor+'_{ch}_card.txt'   .format(ch=charge))
        self.systFile = pathToImport+'/systsFit.txt'

    def getCentralProcesses(self):
        acoeffs = ['ac']+['a{ic}'.format(ic=i) for i in range(8)]
        sig_proc          = ['{boson}{charge}_{ic}_{flav}'.format(charge=self.charge,ic=coeff,flav=self.flavor,boson=self.boson) for coeff in acoeffs]
        other_boson_procs = ['{otherboson}_{flav}_{charge}'.format(flav=self.flavor,charge=self.charge,otherboson='Z' if self._options.wmass else 'W')]
        otherprocs = ['bkg_and_data_{flav}_{charge}'.format(flav=self.flavor,charge=self.charge)] 
        return sig_proc + other_boson_procs + otherprocs

    def getSystList(self):
        ### systematics that are applied to both W and Z
        baseSysts  = ['pdf{i}'.format(i=ipdf) for ipdf in range(1,61)]
        baseSysts += ['{qcdpar}{idir}'.format(qcdpar=par,idir=idir) for par in ['alphaS','muR','muF','muRmuF'] for idir in ['Up','Dn']]     
        ### W/Z mass points
        massPoints = ['mWmass_0']+['mWmass_{sgn}{i}'.format(sgn=sgn,i=imass) for sgn in['m','p'] for imass in range(1,21)]
        ### other systematics
        otherSysts = ['fsr']
        ### build the systematic dictionary for the root files
        systsCards = {self.boson: baseSysts + massPoints + otherSysts,
                      'Z' if self._options.wmass else 'W': baseSysts}
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
                            tmp_var  = copy.deepcopy(proc_hist.Clone(proc_hist.GetName()+'_'+systName))
                            tmp_var.Multiply(tmp_syst)
                            allVars.append(tmp_var)
                            dictProcSyst[(proc, systName)] = None
    
        outfile = ROOT.TFile(options.inputdir+'/{f}_{ch}_shapes.root.multiplierSystematics'.format(f=self.flavor,ch=self.charge), 'RECREATE')
        
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

    def mergeChunks(self):
        ## prepare the relevant files. only the datacards and the correct charge
        allfiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(self._options.inputdir,followlinks=True) for f in fn if f.endswith('.input.root')]
        
        ## get only the central processes
        files = [f for f in allfiles if any(os.path.basename(f)==centralfile+'.input.root' for centralfile in self.centralfiles)]
        files = sorted(files) #, key = lambda x: int(x.rstrip('.input.root').split('_')[1]))
    
        acoeffs = [self.charge+'_ac']+['{charge}_a{ic}'.format(charge=self.charge,ic=i) for i in range(8)]

        tmpfiles = []
        for proc in self.centralfiles:
            boson = proc[0]
            bosonSysts = self.systematics[boson] if boson in self.systematics else []
            ### this should be probably removed from make_wmass_cards.py, but now it is there and should use it
            suffix = {'W': 'sig', 'Z': 'dy'}
            ## the input ROOT file is not built from the inputdir since it may come from different partX subdir
            central_input = [f for f in files if os.path.basename(f)==proc+'.input.root'][0]
            syst_inputs   = sorted([f for f in allfiles if any(os.path.basename(f)=='{proc}_{sfx}_{syst}.input.root'.format(proc=proc,sfx=suffix[boson],syst=syst) for syst in bosonSysts)])
            print 'processing process: ',proc
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
                    # exclude the asimov data_obs which is in any single card, apart the file with the real data
                    if name.endswith('data_obs') and 'data' not in proc: continue

                    name = name.replace('Dn','Down')                    
                    ### remove the _dy_ / _sig_ from the process name which is present only in the mWmass "systematic" variations
                    cleanName = name.replace('_sig_','_').replace('_dy_','_')

                    ### bkg_and_data is in central_input, with no systs. All the systematics (Up/Down) are just copied into the target file
                    if irf==0:
                        if name not in plots:
                            plots[name] = obj.Clone(cleanName)
                            nominals[name] = obj.Clone(cleanName+"0")
                            nominals[name].SetDirectory(None)
                            #print 'replacing old %s with %s' % (name,name)
                            plots[name].Write()
                    else:
                        if any(sysname in name for sysname in ['pdf','fsr']): # these changes by default shape and normalization. Each variation should be symmetrized wrt nominal
                            pfx = '_'.join(name.split("_")[:-2])
                            if 'pdf' in name:
                                patt = re.compile('(pdf)(\d+)')
                                tokens = patt.findall(name)
                                sysname = tokens[0][0]; isys = int(tokens[0][1])
                                cleanName = "{pfx}_{sysname}{isys}".format(pfx=pfx,sysname=sysname,isys=isys)
                                (alternate,mirror) = self.mirrorShape(nominals[pfx],obj,cleanName,self._options.pdfShapeOnly)
                            else:
                                tokens = name.split("_"); pfx = '_'.join(tokens[:-2]); syst = tokens[-1]
                                cleanName = "{pfx}_{syst}".format(pfx=pfx,syst=syst)
                                (alternate,mirror) = self.mirrorShape(nominals[pfx],obj,cleanName,alternateShapeOnly=True)
                            for alt in [alternate,mirror]:
                                if alt.GetName() not in plots:
                                    plots[alt.GetName()] = alt.Clone()
                                    plots[alt.GetName()].Write()
                        elif re.match('.*mWmass.*',name): # this is a special "syst". No need to have Up/Down
                            plots[name] = obj.Clone(cleanName)
                            plots[name].Write()
                        else:
                            tokens = name.split("_"); pfx = '_'.join(tokens[:-2]); syst = tokens[-1].replace('Dn','Down')
                            cleanName = "{pfx}_{syst}".format(pfx=pfx,syst=syst)
                            if name not in plots:
                                plots[name] = obj.Clone(cleanName)
                                plots[name].Write()
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
                if re.match('((m|p)\d+|0)',syst):
                    (proc,syst) = ('_'.join(name.split('_')[1:-2]),'_'.join(name.split('_')[-2:]))
                    histo = None
                else:
                    (proc,syst) = ('_'.join(name.split('_')[1:]),'central')
                    histo = obj.Clone()
                    histo.SetDirectory(None)
            else:
                (proc,syst) = ('_'.join(name.split('_')[1:-1]),name.split('_')[-1])
                histo = None
            procAndHistos[(proc,syst)] = histo
        tf.Close()

        processes   = [proc for (proc,syst) in procAndHistos if (syst=='central' and proc!='data' and proc!='data_obs')]
        signalprocs = sorted([p for p in filter(lambda x: x.startswith('{boson}{charge}'.format(boson=self.boson,charge=self.charge)),processes)])
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
            if 'mWmass' in syst: 
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
                        if group not in systGroups: systGroups[group] = [systName]
                        else:                       
                            if systName not in systGroups[group]: systGroups[group].append(systName)

        if options.mergeRoot or options.remakeMultipliers:
            os.system('hadd -f {of} {of}.baseSystematics {of}.multiplierSystematics'.format(of=self.shapesfile))

        ### now write the datacard
        channel = '{boson}{charge}'.format(boson='W' if self._options.wmass else 'Z',charge=self.charge)
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
        for systName,values in appliedSystMatrix.iteritems():
            ### do not write systematic lines if it is not applied to any process
            if any([v!='-' for v in values]):
                datacard.write( '%-15s   shape %s\n' % (systName," ".join([kpatt % v for v in values])))

        ### -- groups --
        datacard.write('\n\n')
        for group,systs in systGroups.iteritems():
            datacard.write( '%-15s   group = %s\n' % (group," ".join(sorted(systs))) )
        datacard.close()

        print "Done. Datacard in ",self.cardfile

def combineCharges(options):
    ### prepare the combineCards and txt2hdf5 commands
    charges = ['plus','minus']
    singleChargeCards = ['{idir}/{flav}_{charge}_card.txt'.format(idir=options.inputdir,flav=options.flavor,charge=ch) for ch in charges]
    if any([not os.path.exists(card) for card in singleChargeCards]):
        raise RuntimeError, "At least one between %s does not exist. Cannot combine." % " ".join(singleChargeCards)
    print "Cards for + and - ready. Combining them now..."
    combinedCard = '{idir}/{flav}_card.txt'.format(idir=options.inputdir,flav=options.flavor)
    cards = {}
    for ch in charges:
        cards[ch] = '{idir}/{flav}_{charge}_card.txt'.format(idir=os.path.abspath(options.inputdir),flav=options.flavor,charge=ch)
    ccCmd = 'combineCards.py --noDirPrefix '+' '.join(['{ch}={dcfile}'.format(ch=ch,dcfile=card) for ch,card in cards.iteritems()])+' > '+combinedCard
    txt2hdf5Cmd = 'text2hdf5.py {cf} --clipSystVariations {varmax} {pfx}'.format(cf=combinedCard,varmax=1.3,pfx='--postfix '+ options.postfix if len(options.postfix) else '')

    ## run the commands: need cmsenv in the combinetf release
    print ccCmd
    os.system(ccCmd)
    print "Combined card in ",combinedCard
    print txt2hdf5Cmd
    os.system(txt2hdf5Cmd)
    hdf5name=combinedCard.replace('.txt','.hdf5')
    print "hdf5 file in ",hdf5name

    if len(options.postfix):
        hdf5name = hdf5name.replace('.hdf5','_%s.hdf5' % options.postfix)
        combineCmd = 'combinetf.py -t -1 --binByBinStat {hdf5file}'.format(hdf5file=hdf5name)


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
    parser = OptionParser(usage='%prog [self._options] cards/card*.txt')
    parser.add_option('-m','--merge-root', dest='mergeRoot', default=False, action='store_true', help='Merge the root files with the inputs also')
    parser.add_option('-i','--input', dest='inputdir', default='', type='string', help='input directory with all the cards inside')
    parser.add_option('-f','--flavor', dest='flavor', default='mu', type='string', help='lepton flavor (mu,el)')
    parser.add_option('-C','--charge', dest='charge', default='plus,minus', type='string', help='process given charge. default is both')
    parser.add_option(     '--pdf-shape-only'   , dest='pdfShapeOnly' , default=False, action='store_true', help='Normalize the mirroring of the pdfs to central rate.')
    parser.add_option(     '--remake-syst-multipliers'   , dest='remakeMultipliers' , default=False, action='store_true', help='Remake all the systematic multipliers.')
    parser.add_option(     '--comb'   , dest='combineCharges' , default=False, action='store_true', help='Combine W+ and W-, if single cards are done')
    parser.add_option(     '--postfix',    dest='postfix', type="string", default="", help="Postfix for .hdf5 file created with text2hdf5.py when combining charges");
    parser.add_option('--wmass', dest='wmass', action="store_true", default=False, help="Make cards for the wmass analysis. Default is wlike");

    (options, args) = parser.parse_args()
    
    charges = options.charge.split(',')

    if options.remakeMultipliers:
        print 'making all the systematic multipliers'
        makeAllSystematicMultipliers(options.inputdir)
        print 'done making all the systematic multipliers'
    

    for charge in charges:
        cm = CardMaker(options,charge,options.flavor)
        if options.mergeRoot:
            cm.mergeChunks()
        cm.writeDatacard()
    
    if options.combineCharges and len(charges)==2:
        combineCharges(options)
