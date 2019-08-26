#!/bin/env python

import ROOT, os, sys, re, array

dryrun=0
doMuons=1
skipUnpack=1
skipMergeRoot=0
skipSingleCard=0
skipMergeCard=0   # disabled if fitting each charge (see below)
skipMergeCardFlavour=1 # requires both flavours, the electron cards should have all signal bins considered as signal (or be set up manually)
flavourCombinationOutdir = "muElCombination"

# REMEMBER TO ACTIVATE AGAIN THE FSR SYMMETRIZATION 
#
charges = ["plus", "minus"]
#charges = ["plus"]
# to deal with each single charge
fitSingleCharge = 0
singleChargesToFit = ["plus", "minus"]
if fitSingleCharge: 
    skipMergeCard = 1
    skipMergeCardFlavour = 1
#
useXsecWptWeights = 1  
allPtBinsSignal = 1   # usually 1 for muons or combination, 0 for electrons
distinguishNameSigAsBkg = 1 # mainly for electrons to prepare for combination, it gives a different name for pt bins that are treated as background
useBinUncEffStat = False
useBinEtaPtUncorrUncEffStat = False
uncorrelateFakesNuisancesByCharge = False # need to rerun the MergeRoot when changing this one
uncorrelateNuisancesByCharge = ".*CMS_Wmu_muscale.*" # use regular expression, otherwise keep ""
# note that there is always a part of the uncertainty that is charge-uncorrelated
freezePOIs = False
skipFitData = False
skipFitAsimov = False

if not skipMergeCardFlavour:
    allPtBinsSignal = 1

# el
#folder_el = "diffXsec_el_2019_05_13_eta0p2widthFrom1p3_last2p1to2p4/" # keep "/" at the end
#th3file_el = "cards/" + folder_el + "wel_eta0p2widthFrom1p3_last2p1to2p4_fixLepScale_uncorrPtScale.root"
#folder_el = "diffXsec_el_2019_06_21_zptReweight_allPtBinsAsSignal/" # keep "/" at the end
#th3file_el = "cards/" + folder_el + "wel_15June2019_zptReweight.root"
#folder_el = "diffXsec_el_2019_06_21_zptReweight_fixEffStat/" # keep "/" at the end
#folder_el = "diffXsec_el_2019_06_21_zptReweight_allPtBinsAsSignal/" # keep "/" at the end
#th3file_el = "cards/" + folder_el + "wel_13July2019_fixEffStat.root"
#folder_el = "diffXsec_el_2019_07_20_latestScaleFactor_AllIn_IDwithMConlyStat_allPtBinsAsSignal/" # keep "/" at the end
#folder_el = "diffXsec_el_2019_07_20_latestScaleFactor_AllIn_IDwithMConlyStat/"
folder_el = "diffXsec_el_2019_07_20_latestScaleFactor_AllIn_IDwithMConlyStat_NEWPROCESSNAME/"
#folder_el = "diffXsec_el_2019_07_24_testCoarserPt/" # keep "/" at the end
th3file_el = "cards/" + folder_el + "wel_05August_smoothSF_fsrNormGenXsec.root"
#th3file_el = "cards/" + folder_el + "wel_24July2019_testCoarserPt.root"
#folder_el = "diffXsec_el_2019_07_28_testPt2GeV/"
#th3file_el = "cards/" + folder_el + "wel_28July2019_testPt2GeV.root"
# mu
#folder_mu = "diffXsec_mu_2019_04_28_eta0p2widthFrom1p3_last2p1to2p4/" # keep "/" at the end
#th3file_mu = "cards/" + folder_mu + "wmu_eta0p2widthFrom1p3_last2p1to2p4_fixLepScale_uncorrPtScale_addBinUncEffStat.root"
#folder_mu = "diffXsec_mu_2019_05_09_recoEta0p1_recoPt1_genEta0p2from1p3_last2p1to2p4_genPt2/" # keep "/" at the end
#th3file_mu = "cards/" + folder_mu + "wmu_recoEta0p1_recoPt1_genEta0p2from1p3_last2p1to2p4_genPt2.root"
#folder_mu = "diffXsec_mu_2019_06_23_zptReweight_ptReco30/" # keep "/" at the end
#th3file_mu = "cards/" + folder_mu + "wmu_23June2019_zptReweight_ptReco30.root"
#folder_mu = "diffXsec_mu_2019_06_17_zptReweight_chargeUncorrQCDscales/" # keep "/" at the end
#folder_mu = "diffXsec_mu_2019_06_17_zptReweight_chargeUncorrQCDscales_EffStatOnlyStatUnc/" # keep "/" at the end
folder_mu = "diffXsec_mu_2019_06_17_zptReweight_chargeUncorrQCDscales_EffStatOnlyStatUncDataMC/" # keep "/" at the end
#folder_mu = "diffXsec_mu_2019_07_12_noSyst_onlyEtaMinus/"
#folder_mu = "diffXsec_mu_2019_06_17_zptReweight_chargeUncorrQCDscales_unfixedFSRcharge_testBinUncEffStat/" # keep "/" at the end
#folder_mu = "diffXsec_mu_2019_06_17_zptReweight/" # keep "/" at the end
#folder_mu = "diffXsec_mu_2019_08_02_testBinnedSFandUnc/"
#th3file_mu = "cards/" + folder_mu + "wmu_15June2019_zptReweight.root"
#th3file_mu = "cards/" + folder_mu + "wmu_07July2019_zptReweight_unfixedFSRcharge_EffStatOnlyStatUnc.root"
#th3file_mu = "cards/" + folder_mu + "wmu_22July2019_integratedFSR.root"
th3file_mu = "cards/" + folder_mu + "wmu_05August_smoothSF_fsrNormGenXsec.root"
#th3file_mu = "cards/" + folder_mu + "wmu_02August_testBinnedSFandUnc.root"
#th3file_mu = "cards/" + folder_mu + "wmu_04July2019_zptReweight_unfixedFSRcharge_testBinUncEffStat.root"
#th3file_mu = "cards/" + folder_mu + "wmu_12July2019_onlyRecoEtaMinus.root"

folder = folder_mu if doMuons else folder_el
th3file = th3file_mu if doMuons else th3file_el

#================================
# some more things are set below
excludeNuisRegexp = "CMS_DY,CMS_.*FR.*_slope,CMS_.*FR.*_continuous,CMS.*sig_lepeff"
if useBinUncEffStat: 
    excludeNuisRegexp = excludeNuisRegexp + "{comma}".format(comma="," if len(excludeNuisRegexp) else "") + ".*ErfPar\d+EffStat.*,.*BinEtaPtUncorrUncEffStat.*"
elif useBinEtaPtUncorrUncEffStat:
    excludeNuisRegexp = excludeNuisRegexp + "{comma}".format(comma="," if len(excludeNuisRegexp) else "") + ".*ErfPar\d+EffStat.*,.*BinUncEffStat.*"
else:
    excludeNuisRegexp = excludeNuisRegexp + "{comma}".format(comma="," if len(excludeNuisRegexp) else "") + ".*Bin.*UncEffStat.*"

nuisToSymmetrizeVsEta = " '.*fsr.*|.*CMS_W.*scale.*|.*mW.*' "
#nuisToSymmetrizeVsEta = ""

optionsForRootMerger = " --uncorrelate-QCDscales-by-charge --test-eff-syst --etaBordersForFakesUncorr " + ("0.5,1.0,1.5,2.1 " if doMuons else "0.5,1.0,1.5,2.1 ")
binnedSystOpt = " --WZ-testEffSyst-shape '0.0,1.0,1.5' --WZ-ptScaleSyst-shape '0.0,2.1' " if doMuons else " --WZ-testEffSyst-shape '0.0,1.0,1.479,2.0' --WZ-ptScaleSyst-shape '0.0,1.0,1.5,2.1' "
optionsForCardMaker = " --uncorrelate-QCDscales-by-charge --unbinned-QCDscale-Z --sig-out-bkg  "
#optionsForCardMaker += " --exclude-nuisances '.*' --keep-nuisances 'fsr'  "
optionsForCardMaker += " --exclude-nuisances '{expr}' ".format(expr=excludeNuisRegexp) 
optionsForCardMaker += binnedSystOpt 

if len(uncorrelateNuisancesByCharge):
    optionsForRootMerger += " --uncorrelate-nuisances-by-charge '{expr}' ".format(expr=uncorrelateNuisancesByCharge)
    optionsForCardMaker  += " --uncorrelate-nuisances-by-charge '{expr}' ".format(expr=uncorrelateNuisancesByCharge)

if useBinUncEffStat:
    optionsForRootMerger += " --useBinUncEffStat "
    optionsForCardMaker  += " --useBinUncEffStat "
if useBinEtaPtUncorrUncEffStat:
    optionsForRootMerger += " --useBinEtaPtUncorrUncEffStat "
    optionsForCardMaker  += " --useBinEtaPtUncorrUncEffStat "

if useXsecWptWeights: 
    optionsForCardMaker  += " --use-xsec-wpt "

if nuisToSymmetrizeVsEta:
    optionsForRootMerger += " --symmetrize-syst-ratio {regexp} ".format(regexp=nuisToSymmetrizeVsEta)

#--WZ-testEffSyst-LnN 0.012" 
# --wXsecLnN 0.038 # exclude ptslope for fakes, we use that one uncorrelated versus eta 
### --uncorrelate-fakes-by-charge   
# --fakesChargeLnN 0.03 --tauChargeLnN 0.03

postfixCardMaker = "_symFSRptScalemW"
optionsForCardMaker = optionsForCardMaker + " --postfix " + postfixCardMaker

postfixCardMakerMerger = ""
if doMuons: postfixCardMakerMerger = "finalFixes_symFSRptScalemW_ptScaleUncorrCharge"
else      : postfixCardMakerMerger = "finalFixes_sigBkgInAcc_symFSRptScalemW"
optionsForCardMakerMerger = " --postfix " + postfixCardMakerMerger + " --sig-out-bkg  " # --no-combinetf " #--useSciPyMinimizer  " 
if freezePOIs:  
    optionsForCardMakerMerger += " --freezePOIs "
    optionsForCardMaker       += " --freezePOIs "
if skipFitData: 
    optionsForCardMakerMerger += " --skip-fit-data "
    optionsForCardMaker       += " --skip-fit-data "
if skipFitAsimov:
    optionsForCardMakerMerger += " --skip-fit-asimov "
    optionsForCardMaker       += " --skip-fit-asimov "
# --no-correlate-xsec-stat

optionsForCardMakerMergerFlavour = " --postfix combinedLep_finalFixes_elePt01Bkg_bkgNotInGroupOrMaskedChan_scaleXsecPt01by2 --sig-out-bkg  "#--no-text2hdf5 --no-combinetf " # --useSciPyMinimizer "  

#================================

if uncorrelateFakesNuisancesByCharge:
    optionsForRootMerger += " --uncorrelate-fakes-by-charge "
    optionsForCardMaker += " --uncorrelate-fakes-by-charge "

flavour = "mu" if doMuons else "el"
# when creating cards for mu-el combination, all signal bins must be treated in the same way as for muons, so as signal
optionsForCardMaker = optionsForCardMaker  + (" --pt-range-bkg 25.9 30.1  " if not allPtBinsSignal else "") #--eta-range-bkg 1.39 1.61 "
if distinguishNameSigAsBkg:
    optionsForCardMaker = optionsForCardMaker  + " --distinguish-name-sig-as-bkg "

print ""

for charge in charges:

    unpack = "python makeTH1FromTH3.py {th3file} -o cards/{fd} -f {fl} -c {ch} --binfile cards/{fd}binningPtEta.txt".format(th3file=th3file, fd=folder, fl=flavour, ch=charge)
    if not skipUnpack:
        print unpack
        if not dryrun: os.system(unpack)

    mergeRoot = "python mergeRootComponentsDiffXsec.py -f {fl} -c {ch} --indir-bkg  cards/{fd}part0/  --indir-sig cards/{fd} -o cards/{fd} {opt}".format(fd=folder, fl=flavour, ch=charge, opt=optionsForRootMerger)
    if not skipMergeRoot:
        print mergeRoot
        if not dryrun: os.system(mergeRoot)


    makeCards = "python makeCardsFromTH1.py -i cards/{fd} -f {fl} -c {ch} -s w-helicity-13TeV/wmass_{lep}/systsEnv.txt {opt} ".format(fd=folder, fl=flavour, ch=charge, lep="e" if flavour == "el" else "mu", opt=optionsForCardMaker)
    if fitSingleCharge and charge in singleChargesToFit: 
         makeCards += " --fit-single-charge "
    if not skipSingleCard:
        print makeCards
        if not dryrun: os.system(makeCards)


combineCards = "python makeCardsFromTH1.py -i cards/{fd} -f {fl} --comb-charge {opt}".format(fd=folder, fl=flavour, opt=optionsForCardMakerMerger)
if not skipMergeCard:
    print combineCards
    if not dryrun: os.system(combineCards)


combineCardsFlavour = "python makeCardsFromTH1.py --indir-mu cards/{fdmu} --indir-el cards/{fdel} --comb-flavour".format(fdmu=folder_mu, fdel=folder_el) 
combineCardsFlavour += " {opt} -o cards/{combOut}/ ".format(opt=optionsForCardMakerMergerFlavour,combOut=flavourCombinationOutdir)
if not skipMergeCardFlavour:
    print combineCardsFlavour
    if not dryrun: os.system(combineCardsFlavour)

print ""
print "DONE"
print ""

