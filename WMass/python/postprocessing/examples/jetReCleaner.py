import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from CMGTools.WMass.postprocessing.framework.datamodel import Collection
from CMGTools.WMass.postprocessing.framework.eventloop import Module
import ROOT, os

class JetReCleaner(Module):
    def __init__(self,label,jetCollection="Jet", saveJetPrefireSF=False):
        self.label = "" if (label in ["",None]) else ("_"+label)
        self.jetCollection = jetCollection
        # for jet prefiring
        self.saveJetPrefireSF = saveJetPrefireSF        
        self.SFfile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/data/scaleFactors/muons/Map_Jet_L1FinOReff_bxm1_looseJet_SingleMuon_Run2016B-H.root', 'read')
        self.teff = self.SFfile.Get('prefireEfficiencyMap')
        self.eff2d = self.teff.CreateHistogram()

        #self.vars = ("pt","eta","phi","mass","btagCSV")
        self.vars = ("pt","eta","phi","mass","id","rawPt","puId")
        if "jetReCleanerHelper.cc_cc.so" not in ROOT.gSystem.GetLibraries():
            print "Load C++ Worker"
            ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/postprocessing/helpers/jetReCleanerHelper.cc+" % os.environ['CMSSW_BASE'])
        self._worker = ROOT.JetReCleanerHelper(0.4)
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree) # initReaders must be called in beginFile
        self.out = wrappedOutputTree
        self.out.branch("n"+self.jetCollection+self.label, "I")
        for V in self.vars:
            self.out.branch(self.jetCollection+self.label+"_"+V, "F", lenVar="n"+self.jetCollection+self.label)
        if self.saveJetPrefireSF:
            self.out.branch("jetPrefireSF_"+self.jetCollection+self.label, "F")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def initReaders(self,tree): # this function gets the pointers to Value and ArrayReaders and sets them in the C++ worker class
        #for B in "nLepGood", "nJet": setattr(self, B, tree.valueReader(B))
        #for B in "nLepGood", str("n"+self.jetCollection): setattr(self, B, tree.valueReader(B))
        # need generic name for jet collection attached to self, while the reader must be initialized with the proper collection
        setattr(self, "nLepGood", tree.valueReader("nLepGood"))
        setattr(self, "nJetDummy", tree.valueReader("n"+self.jetCollection))
        # should not be needed, isData is already an attribute of the event
        #if self.saveJetPrefireSF and not hasattr(self,"isData"):
        #    setattr(self, "isData", tree.valueReader("isData"))

        for B in "eta", "phi" : setattr(self,"LepGood_"+B, tree.arrayReader("LepGood_"+B))
        for B in "eta", "phi" : setattr(self,"JetDummy_"+B, tree.arrayReader(self.jetCollection+"_"+B))
        self._worker.setLeptons(self.nLepGood,self.LepGood_eta,self.LepGood_phi)
        self._worker.setJets(self.nJetDummy,self.JetDummy_eta,self.JetDummy_phi)
        for v in self.vars:
            if not hasattr(self,"JetDummy_"+v): setattr(self,"JetDummy_"+v, tree.arrayReader(self.jetCollection+"_"+v))
        self._ttreereaderversion = tree._ttreereaderversion # self._ttreereaderversion must be set AFTER all calls to tree.valueReader or tree.arrayReader

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ret={}
        jets = Collection(event,self.jetCollection)
        isDataEvent = event.isData
        for V in self.vars:
            branch = getattr(self, "JetDummy_"+V)
            ret[self.jetCollection+self.label+"_"+V] = [getattr(j,V) for j in jets]
        if event._tree._ttreereaderversion > self._ttreereaderversion: # do this check at every event, as other modules might have read further branches
            self.initReaders(event._tree)
        # do NOT access other branches in python between the check/call to initReaders and the call to C++ worker code
        ## Algo
        cleanJets = self._worker.run()
        ## Output
        self.out.fillBranch('n'+self.jetCollection+self.label, len(cleanJets))
        for V in self.vars:
            self.out.fillBranch(self.jetCollection+self.label+"_"+V, [ ret[self.jetCollection+self.label+"_"+V][j] for j in cleanJets ])
        if self.saveJetPrefireSF:
            if isDataEvent:
                sf = 1.0
            else:
                effs = []
                for j in cleanJets:
                    absjeta = abs(ret[self.jetCollection+self.label+"_eta"][j])
                    jpt = ret[self.jetCollection+self.label+"_pt"][j]
                    # map arrives up to 3.5, and jets above that cannot prefire (even above 3.0, but a jet spread energy in a large cone, so the map is not empty between 3.0 and 3.5)
                    # jet with pt > 30 do not prefire significantly as well
                    if absjeta < 3.5 and jpt > 30.0: 
                        etabin = max(1, min(self.eff2d.GetNbinsX(),self.eff2d.GetXaxis().FindFixBin(absjeta)))
                        ptbin  = max(1, min(self.eff2d.GetNbinsY(),self.eff2d.GetYaxis().FindFixBin(jpt)))
                        eff = self.eff2d.GetBinContent(etabin, ptbin)
                        effs.append(eff)

                if not len(effs):
                    sf = 1.
                elif len(effs) == 1:
                    sf = (1.-effs[0])
                elif len(effs) >= 2:
                    effs = sorted(effs)
                    effs = effs [:2]
                    sf = 1.0 - (effs[0] + effs[1] - effs[0]*effs[1]) # SF = 1-p, p=probability that at least one jet prefires (let's stop at the first 2 with larger prefiring probability)                
            self.out.fillBranch("jetPrefireSF_"+self.jetCollection+self.label, sf)
        # end of --> if self.saveJetPrefireSF:

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

jetReCleaner = lambda : JetReCleaner(label="Clean", saveJetPrefireSF=True)
jetAllReCleaner = lambda : JetReCleaner(label="Clean", jetCollection="JetAll", saveJetPrefireSF=True)
