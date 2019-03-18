#!/bin/env python

import ROOT, os, sys, re, array

dryrun=0
skipUnpack=1
skipMergeRoot=1
skipSingleCard=0
skipMergeCard=0

folder = "diffXsec_el_2019_03_14_ptMax56_dressed_FRpol2Above48GeV/" # keep "/" at the end
th3file = "cards/" + folder + "wel_pt56_L1prefire.root"

#folder = "diffXsec_mu_2019_03_12_ptMax56_dressed/" # keep "/" at the end
#th3file = "cards/" + folder + "wmu_pt56_L1prefire.root"

optionsForRootMerger = " --etaBordersForFakesUncorr 0.5,1.0,1.6,2.0 " # use 0.5,1.0,1.5,2.0 for muons, where eta bins are 0.1 wide
optionsForCardMaker = " --unbinned-QCDscale-Z  --sig-out-bkg  --tauChargeLnN 0.03 --exclude-nuisances 'CMS_DY,CMS_.*FR.*_slope,CMS_.*FR.*_continuous,FakesPtNorm.*' " # --wXsecLnN 0.038 # exclude ptslope for fakes, we use that one uncorrelated versus eta
optionsForCardMakerMerger = " --postfix noFakesPtNorm "  # --no-text2hdf5

# folder = "diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3/"  # keep "/" at the end
# th3file = "cards/" + folder + "wmass_varhists_ele_pt1GeV_eta0p2From1p3.root"
# optionsForRootMerger = " --etaBordersForFakesUncorr 0.5,1.0,1.5,2.1 " # use 0.5,1.0,1.5,2.0 for muons, where eta bins are 0.1 wide
# optionsForCardMaker = " --unbinned-QCDscale-Z --sig-out-bkg --exclude-nuisances 'CMS_DY' "
# optionsForCardMakerMerger = " --postfix allSyst "  # --no-text2hdf5

flavour = "el" if "_el_" in folder else "mu"
if flavour == "el":
    optionsForCardMaker = optionsForCardMaker + " --pt-range-bkg 25.9 30.1  " #--eta-range-bkg 1.39 1.61 "


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

print ""
print "DONE"
print ""
