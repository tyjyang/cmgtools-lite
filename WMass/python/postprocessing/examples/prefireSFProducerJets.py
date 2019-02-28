import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from CMGTools.WMass.postprocessing.framework.datamodel import Collection
from CMGTools.WMass.postprocessing.framework.eventloop import Module
import ROOT, os

class PrefireSFProducerJets(Module):
    def __init__(self):
        self.SFfile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/data/scaleFactors/muons/Map_Jet_L1FinOReff_bxm1_looseJet_SingleMuon_Run2016B-H.root', 'read')
        self.teff = self.SFfile.Get('prefireEfficiencyMap')
        self.eff2d = self.teff.CreateHistogram()
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        ### self.initReaders(inputTree) # initReaders must be called in beginFile
        self.out = wrappedOutputTree
        self.out.branch("jetsPrefireSF", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        #print 'new event\n===================================='


        if not event.isData:

            ##leps = Collection(event, "LepGood")
            jets = Collection(event, "Jet")
            jfwd = Collection(event, "JetFwd")

            effs = []
            for j in jets:
                etabin = max(1, min(self.eff2d.GetNbinsX(),self.eff2d.GetXaxis().FindFixBin(abs(j.eta))))
                ptbin  = max(1, min(self.eff2d.GetNbinsY(),self.eff2d.GetYaxis().FindFixBin(j.pt)))
                eff = self.eff2d.GetBinContent(etabin, ptbin)
                effs.append(eff)

            for j in jfwd:
                etabin = max(1, min(self.eff2d.GetNbinsX(),self.eff2d.GetXaxis().FindFixBin(abs(j.eta))))
                ptbin  = max(1, min(self.eff2d.GetNbinsY(),self.eff2d.GetYaxis().FindFixBin(j.pt)))
                eff = self.eff2d.GetBinContent(etabin, ptbin)
                effs.append(eff)

            if not len(effs):
                sf = 1.
            elif len(effs) == 1:
                sf = (1.-effs[0])
            elif len(effs) >= 2:
                effs = sorted(effs)
                effs = effs [:2]
                sf = (1.-effs[0])*(1.-effs[1])
        else:
            sf = 1.

        
        ## Output
        self.out.fillBranch('jetsPrefireSF', sf)
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

prefireSFProducerJets = lambda : PrefireSFProducerJets()
