import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from CMGTools.WMass.postprocessing.framework.datamodel import Collection
from CMGTools.WMass.postprocessing.framework.eventloop import Module
import ROOT, os

class JetReCleaner(Module):
    def __init__(self,label,jetCollection="Jet"):
        self.label = "" if (label in ["",None]) else ("_"+label)
        self.jetCollection = jetCollection
        #self.vars = ("pt","eta","phi","mass","btagCSV")
        self.vars = ("pt","eta","phi","mass")
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
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def initReaders(self,tree): # this function gets the pointers to Value and ArrayReaders and sets them in the C++ worker class
        #for B in "nLepGood", "nJet": setattr(self, B, tree.valueReader(B))
        #for B in "nLepGood", str("n"+self.jetCollection): setattr(self, B, tree.valueReader(B))
        # need generic name for jet collection attached to self, while the reader must be initialized with the proper collection
        setattr(self, "nLepGood", tree.valueReader("nLepGood"))
        setattr(self, "nJetDummy", tree.valueReader("n"+self.jetCollection))
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
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

jetReCleaner = lambda : JetReCleaner(label="Clean")
jetAllReCleaner = lambda : JetReCleaner(label="Clean", jetCollection="JetAll")
