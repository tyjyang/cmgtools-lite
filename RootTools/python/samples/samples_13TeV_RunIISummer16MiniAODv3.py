import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
kreator = ComponentCreator()


#ZJToMuMu_powhegMiNNLO_pythia8_testProd        = kreator.makeMCComponent("ZJToMuMu_powhegMiNNLO_pythia8_testProd"       , "/ZJToMuMu_TuneCUETP8M1_13TeV-powheg-MiNNLO-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"       , "CMS", ".*root", 2008.4)
#ZJToMuMu_powhegMiNNLO_pythia8_photos_testProd = kreator.makeMCComponent("ZJToMuMu_powhegMiNNLO_pythia8_photos_testProd", "/ZJToMuMu_TuneCUETP8M1_13TeV-powheg-MiNNLO-pythia8-photos/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "CMS", ".*root", 2008.4)



WJetsToLNu_94X        = kreator.makeMCComponent("WJetsToLNu_94X"       , "/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"       , "CMS", ".*root", 3* 20508.9)

WJetsToLNu_94X_ext        = kreator.makeMCComponent("WJetsToLNu_94X_ext"       , "/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM"       , "CMS", ".*root", 3* 20508.9)

WamcAtNLO = [WJetsToLNu_94X, WJetsToLNu_94X_ext]

ZJToMuMu_mWPilot_powhegMiNNLO_pythia8_photos = kreator.makeMCComponent("ZJToMuMu_mWPilot_powhegMiNNLO_pythia8_photos", "/ZJToMuMu_mWPilot_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "CMS", ".*root", 1990.0)


# for private usage with some files
#myFiles = [f.strip() for f in open("%s/src/CMGTools/WMass/cfg/ZJToMuMu_mWPilot_listOn28April2020_first100files.txt" % os.environ['CMSSW_BASE'], "r")]
#PARTIAL_ZJToMuMu_mWPilot = kreator.makePrivateMCComponent('PARTIAL_ZJToMuMu_mWPilot',  '/ZJToMuMu_mWPilot_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', myFiles, 1990.0 )


#######################
# top (cross section taken from samples_13TeV_RunIISummer16MiniAODv2.py when possible)
# just use TT->2l2n and TT->lnqq, choosing powheg when available

#TTJets = kreator.makeMCComponent("TTJets", "/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 831.76)

#TTJets_SingleLeptonFromTbar     = kreator.makeMCComponent("TTJets_SingleLeptonFromTbar"    , "/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 831.76*(3*0.108)*(1-3*0.108) )
#TTJets_SingleLeptonFromTbar_ext = kreator.makeMCComponent("TTJets_SingleLeptonFromTbar_ext", "/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "CMS", ".*root", 831.76*(3*0.108)*(1-3*0.108) )

#TTJets_SingleLeptonFromT     = kreator.makeMCComponent("TTJets_SingleLeptonFromTbar"    , "/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 831.76*(3*0.108)*(1-3*0.108) )

#TTJets_SingleLeptonFromT_ext     = kreator.makeMCComponent("TTJets_SingleLeptonFromTbar"    , "/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "CMS", ".*root", 831.76*(3*0.108)*(1-3*0.108) )

# no fully hadronic, not needed
TTTo2L2Nu_powheg     = kreator.makeMCComponent("TTTo2L2Nu_powheg"    , "/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 831.76*(3*0.108)*(3*0.108) )

TTToSemilepton_powheg     = kreator.makeMCComponent("TTToSemilepton_powheg"    , "/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 2.0*831.76*(3*0.108)*(1-3*0.108) )

# now single top (no fully hadronic decays when possible)
TToLeptons_sch_amcatnlo = kreator.makeMCComponent("TToLeptons_sch_amcatnlo", "/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "CMS", ".*root", 3.36 )

TBarToLeptons_tch_powheg = kreator.makeMCComponent("TBarToLeptons_tch_powheg", "/ST_t-channel_antitop_5f_leptDecays_TuneCUETP8M1_13TeV-powhegV2-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "CMS", ".*root", 80.95*0.108*3)

TToLeptons_tch_powheg = kreator.makeMCComponent("TToLeptons_tch_powheg", "/ST_t-channel_top_5f_leptDecays_TuneCUETP8M1_13TeV-powhegV2-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "CMS", ".*root", 136.02*0.108*3)

T_tWch_noFullyHad_powheg = kreator.makeMCComponent("T_tWch_noFullyHad_powheg"     , "/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "CMS", ".*root",19.55)

TBar_tWch_noFullyHad_powheg      = kreator.makeMCComponent("TBar_tWch_noFullyHad_powheg"     , "/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "CMS", ".*root",19.55)

Top = [
    TTTo2L2Nu_powheg,
    TTToSemilepton_powheg,
    TToLeptons_sch_amcatnlo,
    TBarToLeptons_tch_powheg,
    TToLeptons_tch_powheg,
    T_tWch_noFullyHad_powheg,
    TBar_tWch_noFullyHad_powheg
]

#### QCD mu enriched
# qcd muenr
QCD_Mu15 = kreator.makeMCComponent("QCD_Mu15", "/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 720.65e6*0.00042)

#######################
# Dibosons (cross section taken from samples_13TeV_RunIISummer16MiniAODv2.py when possible)
WW = kreator.makeMCComponent("WW", "/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 63.21 * 1.82)
WW_ext = kreator.makeMCComponent("WW_ext", "/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "CMS", ".*root", 63.21 * 1.82)
WZ = kreator.makeMCComponent("WZ", "/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 47.13)
WZ_ext = kreator.makeMCComponent("WZ_ext", "/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "CMS", ".*root", 47.13)
ZZ = kreator.makeMCComponent("ZZ", "/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 16.523)
ZZ_ext = kreator.makeMCComponent("ZZ_ext", "/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "CMS", ".*root", 16.523)

DiBosons = [
    WW,
    WW_ext,
    WZ,
    WZ_ext,
    ZZ,
    ZZ_ext
]


# cross section from https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns#Diboson

# WWTo2L2Nu = kreator.makeMCComponent("WWTo2L2Nu", "/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 10.481 )
# WWToLNuQQ = kreator.makeMCComponent("WWToLNuQQ", "/WWToLNuQQ_13TeV-powheg/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 43.53 )
# WWToLNuQQ_ext = kreator.makeMCComponent("WWToLNuQQ_ext", "/WWToLNuQQ_13TeV-powheg/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "CMS", ".*root", 43.53 )

# ZZTo2L2Nu = kreator.makeMCComponent("ZZTo2L2Nu", "/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 0.564)
# ZZTo2L2Nu_ext = kreator.makeMCComponent("ZZTo2L2Nu_ext", "/ZZTo2L2Nu_13TeV_powheg_pythia8_ext1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 0.564)
# ZZTo2L2Q = kreator.makeMCComponent("ZZTo2L2Q", "/ZZTo2L2Q_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root",  3.28)
# ZZTo4L = kreator.makeMCComponent("ZZTo4L", "/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "CMS", ".*root", 1.256)
# ZZTo4L_ext = kreator.makeMCComponent("ZZTo4L_ext", "/ZZTo4L_13TeV_powheg_pythia8_ext1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 1.256)

# WZTo1L3Nu   = kreator.makeMCComponent("WZTo1L3Nu"  , "/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"  , "CMS", ".*root", (47.13)*(3*0.108)*(0.2) )
# WZTo1L1Nu2Q = kreator.makeMCComponent("WZTo1L1Nu2Q", "/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root",  10.71)
# WZTo2L2Q    = kreator.makeMCComponent("WZTo2L2Q"   , "/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2//MINIAODSIM"   , "CMS", ".*root",  5.60)
# WZTo3LNu    = kreator.makeMCComponent("WZTo3LNu"   , "/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "CMS", ".*root", 4.42965)
# WZTo3LNu_ext    = kreator.makeMCComponent("WZTo3LNu_ext"   , "/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM", "CMS", ".*root", 4.42965)

### ----------------------------- summary ----------------------------------------

#mcSamples = [ZJToMuMu_powhegMiNNLO_pythia8_testProd, ZJToMuMu_powhegMiNNLO_pythia8_photos_testProd ]
mcSamples = [ZJToMuMu_mWPilot_powhegMiNNLO_pythia8_photos] + WamcAtNLO + Top + DiBosons + [QCD_Mu15]

samples = mcSamples

### ---------------------------------------------------------------------

from CMGTools.TTHAnalysis.setup.Efficiencies import *
dataDir = "$CMSSW_BASE/src/CMGTools/TTHAnalysis/data"

#Define splitting
for comp in mcSamples:
    comp.isMC = True
    comp.isData = False
    comp.splitFactor = 250

if __name__ == "__main__":
    from CMGTools.RootTools.samples.tools import runMain
    runMain(samples,localobjs=locals())
