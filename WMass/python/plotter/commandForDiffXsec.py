#!/bin/env python

import ROOT, os, sys, re, array

dryrun=0
doMuons=1
skipUnpack=0
skipMergeRoot=0
skipSingleCard=0
skipMergeCard=0
skipMergeCardFlavour=1 # requires both flavours, and the electron cards must have all signal bins considered as signal

allPtBinsSignalElectron = 1

# el
folder_el = "diffXsec_el_2019_03_14_ptMax56_dressed_FRpol2Above48GeV/" # keep "/" at the end
th3file_el = "cards/" + folder_el + "wel_pt56_L1prefire_fixEffStat.root"
# mu
folder_mu = "diffXsec_mu_2019_03_12_ptMax56_dressed/" # keep "/" at the end
th3file_mu = "cards/" + folder_mu + "wmu_pt56_L1prefire_fixEffStat.root"

folder = folder_mu if doMuons else folder_el
th3file = th3file_mu if doMuons else th3file_el

uncorrelateFakesNuisancesByCharge = False # need to rerun the MergeRoot when changing this one
# note that there is always a part of the uncertainty that is charge-uncorrelated

#================================
# some more things are set below

optionsForRootMerger = " --etaBordersForFakesUncorr " + ("0.5,1.0,1.4,1.6,2.0 " if doMuons else "0.5,1.0,1.4,1.6,2.0 ")

optionsForCardMaker = " --unbinned-QCDscale-Z --sig-out-bkg  --tauChargeLnN 0.03 --exclude-nuisances 'CMS_DY,CMS_.*FR.*_slope,CMS_.*FR.*_continuous'  " # --wXsecLnN 0.038 # exclude ptslope for fakes, we use that one uncorrelated versus eta ### --uncorrelate-fakes-by-charge   # .*FakesPtNormUncorrelated.* --fakesChargeLnN 0.03

optionsForCardMakerMerger = " --postfix combinedLep --sig-out-bkg --useSciPyMinimizer " #--no-text2hdf5 --no-combinetf " 

optionsForCardMakerMergerFlavour = " --postfix combinedLep --useSciPyMinimizer " # --no-text2hdf5 --no-combinetf "  

#================================

if uncorrelateFakesNuisancesByCharge:
    optionsForRootMerger += " --uncorrelate-fakes-by-charge "
    optionsForCardMaker += " --uncorrelate-fakes-by-charge "

flavour = "mu" if doMuons else "el"
if flavour == "el":
    # when creating cards for mu-el combination, all signal bins must be treated in the same way as for muons, so as signal
    optionsForCardMaker = optionsForCardMaker  + (" --pt-range-bkg 25.9 30.1  " if not allPtBinsSignalElectron else "") #--eta-range-bkg 1.39 1.61 "

charges = ["plus", "minus"]

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
    if not skipSingleCard:
        print makeCards
        if not dryrun: os.system(makeCards)


combineCards = "python makeCardsFromTH1.py -i cards/{fd} -f {fl} --comb-charge {opt}".format(fd=folder, fl=flavour, opt=optionsForCardMakerMerger)
if not skipMergeCard:
    print combineCards
    if not dryrun: os.system(combineCards)


combineCardsFlavour = "python makeCardsFromTH1.py --indir-mu cards/{fdmu} --indir-el cards/{fdel} --comb-flavour".format(fdmu=folder_mu, fdel=folder_el) 
combineCardsFlavour += " --comb-flavour {opt} -o cards/muElCombination/ ".format(opt=optionsForCardMakerMergerFlavour)
if not skipMergeCardFlavour:
    print combineCardsFlavour
    if not dryrun: os.system(combineCardsFlavour)

print ""
print "DONE"
print ""
