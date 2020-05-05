import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
kreator = ComponentCreator()


# TTbar cross section: NNLO, https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO (172.5)
ZJToMuMu_powhegMiNNLO_pythia8_testProd        = kreator.makeMCComponent("ZJToMuMu_powhegMiNNLO_pythia8_testProd"       , "/ZJToMuMu_TuneCUETP8M1_13TeV-powheg-MiNNLO-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"       , "CMS", ".*root", 2008.4)
ZJToMuMu_powhegMiNNLO_pythia8_photos_testProd = kreator.makeMCComponent("ZJToMuMu_powhegMiNNLO_pythia8_photos_testProd", "/ZJToMuMu_TuneCUETP8M1_13TeV-powheg-MiNNLO-pythia8-photos/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "CMS", ".*root", 2008.4)



WJetsToLNu_94X        = kreator.makeMCComponent("WJetsToLNu_94X"       , "/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"       , "CMS", ".*root", 3* 20508.9)

WJetsToLNu_94X_ext        = kreator.makeMCComponent("WJetsToLNu_94X_ext"       , "/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM"       , "CMS", ".*root", 3* 20508.9)

ZJToMuMu_mWPilot_powhegMiNNLO_pythia8_photos = kreator.makeMCComponent("ZJToMuMu_mWPilot_powhegMiNNLO_pythia8_photos", "/ZJToMuMu_mWPilot_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "CMS", ".*root", 1990.0)


# for private usage with some files
myFiles = [f.strip() for f in open("%s/src/CMGTools/WMass/cfg/ZJToMuMu_mWPilot_listOn28April2020_first100files.txt" % os.environ['CMSSW_BASE'], "r")]
PARTIAL_ZJToMuMu_mWPilot = kreator.makePrivateMCComponent('PARTIAL_ZJToMuMu_mWPilot',  '/ZJToMuMu_mWPilot_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', myFiles, 1990.0 )


### ----------------------------- summary ----------------------------------------

#mcSamples = [ZJToMuMu_powhegMiNNLO_pythia8_testProd, ZJToMuMu_powhegMiNNLO_pythia8_photos_testProd ]
mcSamples = [ZJToMuMu_mWPilot_powhegMiNNLO_pythia8_photos ]

samples = mcSamples

### ---------------------------------------------------------------------

from CMGTools.TTHAnalysis.setup.Efficiencies import *
dataDir = "$CMSSW_BASE/src/CMGTools/TTHAnalysis/data"

#Define splitting
for comp in mcSamples:
    comp.isMC = True
    comp.isData = False
    comp.splitFactor = 250 #  if comp.name in [ "WJets", "DY3JetsM50", "DY4JetsM50","W1Jets","W2Jets","W3Jets","W4Jets","TTJetsHad" ] else 100
    comp.puFileMC=dataDir+"/puProfile_Summer12_53X.root"
    comp.puFileData=dataDir+"/puProfile_Data12.root"
    comp.efficiency = eff2012

if __name__ == "__main__":
    from CMGTools.RootTools.samples.tools import runMain
    runMain(samples,localobjs=locals())
