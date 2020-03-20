#!/bin/env python

import ROOT, os, sys, re, array

dryrun=0
doMuons=1
skipUnpack=1
skipMergeRoot=1
skipSingleCard=0
skipMergeCard=0 # disabled if fitting each charge (see below)
skipMergeCardFlavour=1 # requires both flavours, the electron cards should have all signal bins considered as signal (or be set up manually)
flavourCombinationOutdir = "muElCombination_allSig_nativeMCatNLOxsec"
#flavourCombinationOutdir = "muElCombination_allSig_nativeMCatNLOxsec_1sigBin_4fixedPOI"
#flavourCombinationOutdir = "muElCombination_allSig_nativeMCatNLOxsec_1sigBin_4fixedPOI_ptMax45"

doImpactsOnMW=True
# exclude some nuisances
# some more things are set below
excludeNuisRegexp = "CMS_DY,CMS_.*FR.*_slope,CMS_.*FR.*_continuous,CMS.*sig_lepeff" #,.*pdf.*|.*alphaS.*|.*mu(R|F).*|.*mW.*"


combineElePt01asBkg = 0
# if True, it requires existence of cards with distinguished name for POIs treated as background
# note that electron folder is modified below
if not skipMergeCardFlavour:
    if combineElePt01asBkg:
        flavourCombinationOutdir = "muElCombination_elePt01bkg_nativeMCatNLOxsec"
    if not skipMergeRoot or not skipSingleCard or not skipMergeCard:
        print "Warning: skipMergeCardFlavour is deactivated, as well as other options. Please check."
        quit()
else:
    combineElePt01asBkg = 0

flavour = "mu" if doMuons else "el"
print "="*30
print "flavour: %s" % flavour
print "="*30

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
#
# manage usage of pt scales
useNativeMCatNLOxsecW = 1 # use 60400 pb as xsec for W and Tau, instead of 3*20508.9 (signal should already be made with that, so it is done on tau scaling histograms in the merger, and in the fit the proper xsec file is chosen)
useSmoothPtScales = 0 # new smooth pt scales defined in mergeRootComponentsDiffXsec.py with addSmoothLepScaleSyst
addSmoothPtScalesWithStandardOnExtremePtBins = 1  # use option --add-smooth-ptscale-extremePtFromStandard for merger
# previous option requires useSmoothPtScales to be set to 0: basically this options add the smooth pt scales in the file, but still produces the old ones in the merger. Then, in the fits, it uses the new smooth ones
useAnalyticSmoothPtScales = 1 #if doMuons else 0 # suggested option: it disables the other ones for the pt scales (for muons, the Rochester corrections are used)
# for electrons, old functions behave better

clipSystVariations = 1.3  # use <=0 not to use it. To clip at XX% pass 1+XX
clipSystVariationsSignal = 1.15
uncorrelateFSRbyFlavour = 1
noPtScalesOutliers = 0  # 0 means using the scales on outliers
profileEleScales = 1
useExpNonProfiledErrs = 0 # for now only works for ele when not profiling scales, ineffective if profiling everything
useXsecWptWeights = 1 
allPtBinsSignal = 1   # usually 1 for muons or combination, 0 for electrons
distinguishNameSigAsBkg = 1 # mainly for electrons to prepare for combination, it gives a different name for pt bins that are treated as background (can stay true for muons, because in the merger the name is changed only if allPtBinsSignal = 0)
useBinUncEffStat = False
useBinEtaPtUncorrUncEffStat = False
uncorrelateFakesNuisancesByCharge = False # need to rerun the MergeRoot when changing this one
uncorrelatePtScalesByCharge = False if useAnalyticSmoothPtScales else True if doMuons else False
uncorrelateQCDscalesByCharge = True
uncorrelateNuisancesByCharge = ""  # use regular expression, otherwise keep ""
uncorrelatePtscaleByEtaside = False if useAnalyticSmoothPtScales else True
# note that there is always a part of the uncertainty that is charge-uncorrelated
freezePOIs = True  # no need for dedicated postfix, added automatically to the one given below
skipFitData = False
skipFitAsimov = False

if not skipMergeCardFlavour:
    allPtBinsSignal = 1

# el
#folder_el = "diffXsec_el_2019_09_22_nativeMCatNLOxsec/"
folder_el = "diffXsec_el_2019_09_22_nativeMCatNLOxsec_allPtBinsAsSignal/"
th3file_el = "cards/" + folder_el + "wel_20Sept2019_smoothSF_fsrNormGenXsec_WptNormGenXsec_nativeMCatNLOxsec.root"
#folder_el = "diffXsec_el_2019_09_22_nativeMCatNLOxsec_allPtBinsAsSignal_1sigBin_4fixedPOI/"
#th3file_el = "cards/" + folder_el + "wel_15oct2019_smoothSF_fsrNormGenXsec_WptNormGenXsec_nativeMCatNLOxsec_singleSignalBin.root"
#folder_el = "diffXsec_el_2019_09_22_onlyZandTau_recoPt29to46p5_1sigBin_4fixedPOI/"
#th3file_el = "cards/" + folder_el + "wel_15oct2019_smoothSF_fsrNormGenXsec_WptNormGenXsec_nativeMCatNLOxsec_singleSignalBin_recoPt29to46p5.root"
#folder_el = "diffXsec_el_2019_09_22_nativeMCatNLOxsec_allPtBinsAsSignal_1sigBin_4fixedPOI_ptMax45/"
#th3file_el = "cards/" + folder_el + "wel_15oct2019_smoothSF_fsrNormGenXsec_WptNormGenXsec_nativeMCatNLOxsec_singleSignalBin_ptMax45.root"

#folder_el = "diffXsec_el_2019_09_22_onlyZandTau_recoPt29to57_1sigBin_4fixedPOI/"
#th3file_el = "cards/" + folder_el + "wel_15oct2019_smoothSF_fsrNormGenXsec_WptNormGenXsec_nativeMCatNLOxsec_singleSignalBin_recoPt29to57.root"


# mu
folder_mu = "diffXsec_mu_2019_09_19_nativeMCatNLOxsec/"
th3file_mu = "cards/" + folder_mu + "wmu_19Sept2019_smoothSF_fsrNormGenXsec_WptNormGenXsec_nativeMCatNLOxsec.root"
#folder_mu = "diffXsec_mu_2019_09_19_nativeMCatNLOxsec_1sigBin_4fixedPOI/"
#th3file_mu = "cards/" + folder_mu + "wmu_15oct2019_smoothSF_fsrNormGenXsec_WptNormGenXsec_nativeMCatNLOxsec_singleSignalBin.root"
#folder_mu = "diffXsec_mu_2019_09_19_onlyZandTau_recoPt25to57_1sigBin_4fixedPOI/"
#th3file_mu = "cards/" + folder_mu + "wmu_15oct2019_smoothSF_fsrNormGenXsec_WptNormGenXsec_nativeMCatNLOxsec_singleSignalBin_recoPt25to57.root"
#folder_mu = "diffXsec_mu_2019_09_19_nativeMCatNLOxsec_1sigBin_4fixedPOI_ptMax45/"
#th3file_mu = "cards/" + folder_mu + "wmu_15oct2019_smoothSF_fsrNormGenXsec_WptNormGenXsec_nativeMCatNLOxsec_singleSignalBin_ptMax45.root"
#folder_mu = "diffXsec_mu_2019_09_19_onlyZandTau_recoPt25to46p5_1sigBin_4fixedPOI/"
#th3file_mu = "cards/" + folder_mu + "wmu_15oct2019_smoothSF_fsrNormGenXsec_WptNormGenXsec_nativeMCatNLOxsec_singleSignalBin_recoPt29to46p5.root"

if not skipMergeCardFlavour and combineElePt01asBkg:
    folder_el = "diffXsec_el_2019_09_22_nativeMCatNLOxsec_LowPtBkg4Comb/"
    folder_mu = "diffXsec_mu_2019_09_19_nativeMCatNLOxsec_4comb/"


folder = folder_mu if doMuons else folder_el
th3file = th3file_mu if doMuons else th3file_el

if addSmoothPtScalesWithStandardOnExtremePtBins and useSmoothPtScales:
    print "Warning: conflicting flags addSmoothPtScalesWithStandardOnExtremePtBins and useSmoothPtScales. Abort"
    quit()


#================================
if useSmoothPtScales or addSmoothPtScalesWithStandardOnExtremePtBins or useAnalyticSmoothPtScales:
    excludeNuisRegexp += ",.*CMS_W.*scale.*"
if useBinUncEffStat: 
    excludeNuisRegexp = excludeNuisRegexp + "{comma}".format(comma="," if len(excludeNuisRegexp) else "") + ".*ErfPar\d+EffStat.*,.*BinEtaPtUncorrUncEffStat.*"
elif useBinEtaPtUncorrUncEffStat:
    excludeNuisRegexp = excludeNuisRegexp + "{comma}".format(comma="," if len(excludeNuisRegexp) else "") + ".*ErfPar\d+EffStat.*,.*BinUncEffStat.*"
else:
    excludeNuisRegexp = excludeNuisRegexp + "{comma}".format(comma="," if len(excludeNuisRegexp) else "") + ".*Bin.*UncEffStat.*"

# pt scales are uncorrelated vs eta side. Is the symetrization necessary? Maybe not if they are the smoothed ones
nuisToSymmetrizeVsEta = ".*fsr.*|.*mW.*"
if useSmoothPtScales or useAnalyticSmoothPtScales:
    pass
    #nuisToSymmetrizeVsEta += "|.*smooth.*scale\d+.*"
else:
    nuisToSymmetrizeVsEta += "|.*CMS_W.*scale\d+.*"
#nuisToSymmetrizeVsEta = ""

optionsForRootMerger = " --test-eff-syst --etaBordersForFakesUncorr " + ("0.5,1.0,1.5,2.1 " if doMuons else "0.5,1.0,1.5,2.1 ")
#binnedSystOpt = " --WZ-testEffSyst-shape '0.0,1.0,1.5' --WZ-ptScaleSyst-shape '0.0,2.1' " if doMuons else " --WZ-testEffSyst-shape '0.0,1.0,1.5,2.1' --WZ-ptScaleSyst-shape '0.0,1.0,1.5,2.1' "
binnedSystOpt = " --WZ-testEffSyst-shape '0.0,1.0,1.5' " if doMuons else " --WZ-testEffSyst-shape '0.0,1.0,1.5,1.9,2.1' "

if not useAnalyticSmoothPtScales:
    if doMuons:
        binnedSystOpt += " --WZ-ptScaleSyst-shape '0.0,2.1' "
    else:
        binnedSystOpt += " --WZ-ptScaleSyst-shape '0.0,1.0,1.5,2.1' "


optionsForCardMaker = " --unbinned-QCDscale-Z --sig-out-bkg  "
#optionsForCardMaker += " --exclude-nuisances '.*' --keep-nuisances 'fsr'  "
optionsForCardMaker += " --exclude-nuisances '{expr}' ".format(expr=excludeNuisRegexp) 
optionsForCardMaker += binnedSystOpt 

#optionsForCardMaker += " --wAllXsecLnN '0.05/0.03' "  # justified from first principles, because our decorrelation scheme for the qcd scales otherwise reduces the overall xsec uncertainty unreasonably

if useNativeMCatNLOxsecW:
    optionsForRootMerger += " --use-native-MCatNLO-xsec "
    optionsForCardMaker  += " --use-native-MCatNLO-xsec "

if uncorrelateQCDscalesByCharge:
    optionsForRootMerger += " --uncorrelate-QCDscales-by-charge "
    optionsForCardMaker  += " --uncorrelate-QCDscales-by-charge "    

if len(uncorrelateNuisancesByCharge):
    optionsForRootMerger += " --uncorrelate-nuisances-by-charge '{expr}' ".format(expr=uncorrelateNuisancesByCharge)
    optionsForCardMaker  += " --uncorrelate-nuisances-by-charge '{expr}' ".format(expr=uncorrelateNuisancesByCharge)

if useAnalyticSmoothPtScales:
    optionsForRootMerger += " --use-analytic-smooth-pt-scales "
    if flavour == "mu":
        optionsForCardMaker  += " --use-analytic-smooth-pt-scales " 
    else:
        optionsForCardMaker  += " --use-analytic-smooth-pt-scales " 
        #optionsForCardMaker  += " --use-smooth-ptscale " 
else:
    if useSmoothPtScales:
        optionsForRootMerger += " --use-smooth-ptscale "
        optionsForCardMaker  += " --use-smooth-ptscale "
    if addSmoothPtScalesWithStandardOnExtremePtBins:
        optionsForRootMerger += " --add-smooth-ptscale-extremePtFromStandard " # here do not use --use-smooth-ptscale
        optionsForCardMaker  += " --use-smooth-ptscale " # note that here the smooth pt scales should be used in the datacard

if uncorrelatePtScalesByCharge:
    optionsForRootMerger += " --uncorrelate-ptscale-by-charge "
    optionsForCardMaker  += " --uncorrelate-ptscale-by-charge "
if uncorrelatePtscaleByEtaside:
    optionsForRootMerger += " --uncorrelate-ptscale-by-etaside "
    optionsForCardMaker  += " --uncorrelate-ptscale-by-etaside "

if useBinUncEffStat:
    optionsForRootMerger += " --useBinUncEffStat "
    optionsForCardMaker  += " --useBinUncEffStat "
if useBinEtaPtUncorrUncEffStat:
    optionsForRootMerger += " --useBinEtaPtUncorrUncEffStat "
    optionsForCardMaker  += " --useBinEtaPtUncorrUncEffStat "

if useXsecWptWeights: 
    optionsForCardMaker  += " --use-xsec-wpt "

if nuisToSymmetrizeVsEta:
    optionsForRootMerger += " --symmetrize-syst-ratio '{regexp}' ".format(regexp=nuisToSymmetrizeVsEta)

if uncorrelateFSRbyFlavour:
    optionsForRootMerger += " --uncorrelate-fsr-by-flavour"
    optionsForCardMaker  += " --uncorrelate-fsr-by-flavour"

if noPtScalesOutliers:
    optionsForCardMaker  += " --no-pt-scales-outliers "

#--WZ-testEffSyst-LnN 0.012" 
# --wXsecLnN 0.038 # exclude ptslope for fakes, we use that one uncorrelated versus eta 
### --uncorrelate-fakes-by-charge   
# --fakesChargeLnN 0.03 --tauChargeLnN 0.03

#optionsForCardMaker += " --wXsecLnN 0.3 "

postfixCardMaker = "_symFSRptScalemW" # for single-charge fit in single flavor
optionsForCardMaker = optionsForCardMaker + " --postfix " + postfixCardMaker

postfixCardMakerMerger = ""
if doMuons: postfixCardMakerMerger = "nativeMCatNLOxsecW_RochesterCorrUncert_cropNegBinNomi_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_FSRshapeOnly_addInclXsec_addImpactsOnMw"
else      : postfixCardMakerMerger = "nativeMCatNLOxsecW_profilePtScales_newSmoothUncorrScale_cropNegBinNomi_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_noSplitPtSystByPt_FSRshapeOnly_addInclXsec"

if noPtScalesOutliers:
    postfixCardMakerMerger += "_noPtScalesOutliers"

optionsForCardMakerMerger = " --postfix " + postfixCardMakerMerger + " --sig-out-bkg  " # --no-combinetf " #--useSciPyMinimizer  " 


#optionsForCardMakerMergerFlavour = " --postfix combinedLep_elePt01Bkg_bkgNotInGroupOrMaskedChan_symFSRmWptScale_smoothPtScaleUncorrEtaMuElUncorrChargeMuExtremePtFromOld2BinsForOut_LnN0p03Up0p05DownOnAllW --sig-out-bkg --skip-hadd-xsec "
#--no-text2hdf5 --no-combinetf " # " --useSciPyMinimizer " # " --skip-hadd-xsec --just-fit " 
optionsForCardMakerMergerFlavour = " --postfix combinedLep_allSig_nativeMCatNLOxsec_profileLepScale_cropNegBinNomi_uncorrFSRbyFlav_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_noSplitElePtSystByPt_FSRshapeOnly_addInclXsec --sig-out-bkg " #--no-text2hdf5 --no-combinetf " # --useSciPyMinimizer "  

if combineElePt01asBkg:
    optionsForCardMakerMergerFlavour = optionsForCardMakerMergerFlavour.replace("_allSig_","_elePt01bkg_")

if skipFitData: 
    optionsForCardMakerMerger += " --skip-fit-data "
    optionsForCardMaker       += " --skip-fit-data "
    optionsForCardMakerMergerFlavour += " --skip-fit-data "
if skipFitAsimov:
    optionsForCardMakerMerger += " --skip-fit-asimov "
    optionsForCardMaker       += " --skip-fit-asimov "
    optionsForCardMakerMergerFlavour += " --skip-fit-asimov "
# --no-correlate-xsec-stat
if doImpactsOnMW:
    optionsForCardMaker += " --impacts-mW "
    optionsForCardMakerMerger += " --impacts-mW "
    optionsForCardMakerMergerFlavour += " --impacts-mW "


if combineElePt01asBkg:
    optionsForCardMakerMergerFlavour += " --combineElePt01asBkg "

if flavour == "el" and not profileEleScales:
    optionsForCardMaker  += " --no-profile-pt-scales "
    if useExpNonProfiledErrs:
        optionsForCardMaker       += " --useExpNonProfiledErrs "
        optionsForCardMakerMerger += " --useExpNonProfiledErrs "
        optionsForCardMakerMergerFlavour += " --useExpNonProfiledErrs "
#================================

if uncorrelateFakesNuisancesByCharge:
    optionsForRootMerger += " --uncorrelate-fakes-by-charge "
    optionsForCardMaker += " --uncorrelate-fakes-by-charge "

# when creating cards for mu-el combination, all signal bins must be treated in the same way as for muons, so as signal
ptRangeBkg = " --pt-range-bkg 25.9 30.1  " if not allPtBinsSignal else ""
optionsForCardMaker = optionsForCardMaker  + ptRangeBkg #--eta-range-bkg 1.39 1.61 "
if distinguishNameSigAsBkg:
    optionsForRootMerger = optionsForRootMerger + " --distinguish-name-sig-as-bkg " + ptRangeBkg
    optionsForCardMaker  = optionsForCardMaker  + " --distinguish-name-sig-as-bkg "

if freezePOIs:  
    optionsForCardMakerMerger += " --freezePOIs "
    optionsForCardMaker       += " --freezePOIs "
    optionsForCardMakerMergerFlavour += " --freezePOIs "

if clipSystVariationsSignal > 0:
    optionsForCardMaker       += (" --clipSystVariationsSignal " + str(clipSystVariationsSignal))
    optionsForCardMakerMergerFlavour += (" --clipSystVariationsSignal " + str(clipSystVariationsSignal))
if clipSystVariations > 0:
    optionsForCardMaker       += (" --clipSystVariations " + str(clipSystVariations))
    optionsForCardMakerMergerFlavour += (" --clipSystVariations " + str(clipSystVariations))
    

print ""

for charge in charges:

    unpack = "python makeTH1FromTH3.py {th3file} -o cards/{fd} -f {fl} -c {ch} --binfile cards/{fd}binningPtEta.txt".format(th3file=th3file, fd=folder, fl=flavour, ch=charge)
    if not skipUnpack:
        print unpack
        if not dryrun: 
            ret = os.system(unpack)
            if ret:
                print "Error in unpacker"
                quit()
                
    mergeRoot = "python mergeRootComponentsDiffXsec.py -f {fl} -c {ch} --indir-bkg  cards/{fd}part0/  --indir-sig cards/{fd} -o cards/{fd} {opt}".format(fd=folder, fl=flavour, ch=charge, opt=optionsForRootMerger)
    if not skipMergeRoot:
        print mergeRoot
        if not dryrun: 
            ret = os.system(mergeRoot)
            if ret:
                print "Error in merger"
                quit()


    makeCards = "python makeCardsFromTH1.py -i cards/{fd} -f {fl} -c {ch} -s w-helicity-13TeV/wmass_{lep}/systsEnv.txt {opt} ".format(fd=folder, fl=flavour, ch=charge, lep="e" if flavour == "el" else "mu", opt=optionsForCardMaker)
    if fitSingleCharge and charge in singleChargesToFit: 
         makeCards += " --fit-single-charge "
    if not skipSingleCard:
        print makeCards
        if not dryrun: 
            ret = os.system(makeCards)
            if ret:
                print "Error in card maker"
                quit()


combineCards = "python makeCardsFromTH1.py -i cards/{fd} -f {fl} --comb-charge {opt}".format(fd=folder, fl=flavour, opt=optionsForCardMakerMerger)
if not skipMergeCard:
    print combineCards
    if not dryrun: 
        ret = os.system(combineCards)
        if ret:
            print "Error in charge-combined card maker"
            quit()


combineCardsFlavour = "python makeCardsFromTH1.py --indir-mu cards/{fdmu} --indir-el cards/{fdel} --comb-flavour".format(fdmu=folder_mu, fdel=folder_el) 
combineCardsFlavour += " {opt} -o cards/{combOut}/ ".format(opt=optionsForCardMakerMergerFlavour,combOut=flavourCombinationOutdir)
if not skipMergeCardFlavour:
    print combineCardsFlavour
    if not dryrun: 
        ret = os.system(combineCardsFlavour)
        if ret:
            print "Error in flavour-combined card maker"
            quit()


print ""
print "# DONE"
print ""

