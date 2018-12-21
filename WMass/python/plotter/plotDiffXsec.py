#!/bin/env python

import ROOT, os, sys, re, array

# to run plots from Asimov fit and data. For toys need to adapt this script

dryrun = 0
skipData = 0
seed = 123456789
folder = "diffXsec_el_2018_12_18_onlyBkg_pt2GeV_last3GeV_eta0p2From2p0/"
postfix = "allSyst_bbb1_cxs1"
opts = "--lumi-norm 35900.0  -n --palette 55 --hessian"

flavour = "el" if "_el_" in folder else "mu"
lepton = "electron" if flavour == "el"  else "muon"

fits = ["Asimov", "Data"]

print ""

for fit in fits:

    if skipData and fit == "Data": continue

    typedir = "hessian" if fit == "Asimov" else "data"

    command = "python w-helicity-13TeV/plotDiffXsecChargeAsymmetry.py"
    command += " -i cards/{fd} -c {fl}".format(fd=folder, fl=flavour)
    command += " -o plots/diffXsec/chargeAsymmetry/{lep}/{fd}".format(lep=lepton,fd=folder)
    command += " -t cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
    command += " {opt} --suffix {fit}_{pf}".format(opt=opts,fit=fit,pf=postfix)

    if fit == "Data":
        command += " --fit-data --expected-toyfile cards/{fd}/fit/hessian/fitresults_{s}_Asimov_{pf}.root".format(fd=folder,s=seed,pf=postfix)

    print ""
    print command
    if not dryrun:
        os.system(command)




## example
# python w-helicity-13TeV/plotDiffXsecChargeAsymmetry.py -i cards/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3 -o plots/diffXsec/chargeAsymmetry/electron/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3/ -c el -t cards/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3/fit/data/fitresults_123456789_Data_bbb1_cxs1.root  --suffix Data_bbb1_cxs1 --lumi-norm 35900.0  -n --hessian --palette 55 --fit-data --expected-toyfile cards/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3/fit/hessian/fitresults_123456789_Asimov_bbb1_cxs1.root
