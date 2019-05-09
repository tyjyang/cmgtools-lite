#!/bin/env python

import ROOT, os, sys, re, array

# to run plots from Asimov fit and data. For toys need to adapt this script

doMuElComb = 0
dryrun = 0
skipData = 0
onlyData = 1

skipPlot = 1
skipTemplate = 1
skipDiffNuis = 0
skipPostfit = 1  # only for Data
skipCorr = 1
skipImpacts = 1

allPtBinsSignalElectron = 0

seed = 123456789

#folder = "diffXsec_el_2019_04_13_newSystAndWtau/"
#folder = "diffXsec_mu_2019_04_09_newSystAndWtau_fixTriSF/"
folder = "diffXsec_mu_2019_04_28_eta0p2widthFrom1p3_last2p1to2p4/"
#folder = "diffXsec_mu_2019_05_04_etaReco0p1_etaGen0p2from1p3_last2p1to2p4//"
if doMuElComb:
    folder = "muElCombination"

postfix = "testEffSyst_uncorrEta_fixLepScale_uncorrPtScale"
#postfix = "combinedLep"
if doMuElComb:
    postfix = "combinedLep"
postfix += "_bbb1_cxs1"
#postfix += "_bbb0"

flavour = "el" if "_el_" in folder else "mu"
lepton = "electron" if flavour == "el"  else "muon"
if doMuElComb:
    flavour = "lep"
    lepton = "lepton"

fits = ["Asimov", "Data"]

ptBinsSetting = " --pt-range-bkg 25.9 30.1 --pt-range '30,56' " if (flavour == "el" and not allPtBinsSignalElectron) else ""  # " --eta-range-bkg 1.39 1.61 "
ptMinForImpacts = " --pt-min-signal 30" if (flavour == "el" and not allPtBinsSignalElectron) else ""
optTemplate = " --norm-width --draw-selected-etaPt 2.25,39.5 --syst-ratio-range 'template' --palette 57 --do-signal-syst '.*scale1.*|.*lepeff.*' "  # --draw-selected-etaPt 0.45,38 --zmin 10 # kLightTemperature=87
ptMaxTemplate = "56"
ptMinTemplate = "30" if flavour == "el" else "26"

# do not ask Wplus.*_ieta_.*_mu$ to select signal strength rejecting pmasked, because otherwise you must change diffNuisances.py
# currently it uses GetFromHessian with keepGen=True, so _mu$ would create a problem (should implement the possibility to reject a regular expression)
# if you want mu rejecting pmasked do _mu_mu or _el_mu (for electrons _mu works because it doesn't induce ambiguities with the flavour)
diffNuisances_pois = [#"pdf.*|alphaS|mW", 
                      #"muR.*|muF.*", 
                      #"Fakes(Eta|Pt).*[0-9]+mu.*", 
                      #"Fakes(Eta|Pt).*[0-9]+el.*", 
                      #"ErfPar0EffStat.*", 
                      #"ErfPar1EffStat.*", 
                      #"ErfPar2EffStat.*", 
                      "CMS_.*|.*TestEffSyst.*|mW", 
                      #"Wplus.*_ieta_.*_mu",     
                      #"Wminus.*_ieta_.*_mu"
                      ]

# this is appended to nuis below
#correlationSigRegexp = {"Wplus_ieta6ipt6" : ".*Wplus.*_ieta_6_ipt_6_.*mu"  
#                        }
# need to specify which matrix, by default the one for mu is chosen and the bin label looks like Wplus_el_ieta_1_ipt_0_Wplus_el, without ending mu
correlationSigRegexp = {"Wplus_ieta6ipt8" : ".*Wplus_.*_ieta_6_ipt_8_.*"
                        }

correlationNuisRegexp = {# "allPDF"           : "pdf.*", 
                         "somePDFandAlphaS" : "^pdf([1-9]|1[0-9]|20)$,alphaS", 
                         "QCDscales"        : "muR.*|muF.*", 
                         # "muR"              : "^muR[1-9]+", 
                         # "muF"              : "^muF[1-9]+", 
                         # "muRmuF"           : "^muRmuF[1-9]+", 
                         "FakesEtaPtUncorr" : "Fakes(Eta|Pt).*[0-9]+mu.*",
                         "FakesEtaPtUncorr" : "Fakes(Eta|Pt).*[0-9]+el.*", 
                         "CMSsyst"          : "CMS_.*",
                         # "ErfPar0EffStat"   : "ErfPar0EffStat.*",
                         # "ErfPar1EffStat"   : "ErfPar1EffStat.*",
                         # "ErfPar2EffStat"   : "ErfPar2EffStat.*"
                         }

correlationMatrixTitle = {"allPDF"           : "all PDFs", 
                          "somePDFandAlphaS" : "some PDFs + #alpha_{S}", 
                          "QCDscales"        : "all QCD scales", 
                          "muR"              : "#mu_{R} QCD scales", 
                          "muF"              : "#mu_{F} QCD scales", 
                          "muRmuF"           : "#mu_{R}#mu_{F} QCD scales", 
                          "FakesEtaPtUncorr" : "fakes #eta normalizations and p_{T} slopes", 
                          "CMSsyst"          : "some experimental systematics",
                          "ErfPar0EffStat"   : "signal efficiency (uncorr. par 0)",
                          "ErfPar1EffStat"   : "signal efficiency (uncorr. par 1)",
                          "ErfPar2EffStat"   : "signal efficiency (uncorr. par 2)"
                         }

# for impacts
targets = [#"mu", 
           #"xsec", 
           "xsecnorm",
           #"etaptasym",
           #"etaxsec",
           #"etaxsecnorm",
           #"etaasym"
           ]

# impacts_nuis = [".*pdf.*", 
#                 ".*muR.*|.*muF.*", 
#                 ".*ErfPar0EffStat.*", 
#                 ".*ErfPar1EffStat.*", 
#                 ".*ErfPar2EffStat.*", 
#                 ".*CMS_.*",
#                 "GROUP"     # this will do groups, I can filter some of them, but they are few, so I will use --nuisgroups '.*'
#                 ]
impacts_nuis = ["GROUP"]     # this will do groups, I can filter some of them, but they are few, so I will use --nuisgroups '.*'
#groupnames = 'binByBinStat,stat,pdfs,wmodel,EffStat,scales,alphaS'
groupnames = 'binByBinStat,stat,luminosity,pdfs,QCDTheo,Fakes,OtherBkg,OtherExp,EffStat,EffSyst,lepScale'

# no longer used: for impacts in the form of graphs, use "W.*_ieta_.*" for pt-integrated stuff, and "W.*_ieta_.*_ipt_XX" for the rest, where XX is a given pt bin 
impacts_pois = [#"Wplus.*_ipt_2_.*" if flavour == "el" else "Wplus.*_ipt_0_.*",
                #"Wplus.*_ipt_8_.*",
                #"Wplus.*_ipt_11_.*",
                #"W.*_ieta_0_.*",
                #"W.*_ieta_.*",
                #"W.*_ieta_.*_ipt_2_.*",
                #"W.*_ieta_17_.*"
                ]


print ""

for charge in ["plus","minus"]:
    print ""
    print "="*30
    print "TEMPLATE ROLLING for charge " + charge
    print ""

    command = "python w-helicity-13TeV/templateRolling.py"
    command += " cards/{fd} -o plots/diffXsecAnalysis/{lep}/{fd}/templateRolling/{pfx}/ -c {fl}".format(fd=folder, lep=lepton, fl=flavour, pfx=postfix)
    command += " --plot-binned-signal -a diffXsec -C {ch} --pt-range '{ptmin},{ptmax}' ".format(ch=charge, ptmin=ptMinTemplate, ptmax=ptMaxTemplate)
    command += " {opt} ".format(opt=optTemplate)
    if not skipTemplate:
        print ""
        print command
        if not dryrun:
            os.system(command)


for fit in fits:

    if skipData and fit == "Data": continue
    if onlyData and fit != "Data": continue

    typedir = "hessian" if fit == "Asimov" else "data"

    print ""
    print ">"*30
    print ">>>> Running on",fit
    print ">"*30
    print ""

    print ""
    print "="*30
    print "DIFFERENTIAL CROSS SECTION AND CHARGE ASYMMETRY"
    print ""

    ## DIFFERENTIAL CROSS SECTION AND CHARGE ASYMMETRY
    command = "python w-helicity-13TeV/plotDiffXsecChargeAsymmetry.py"
    command += " -i cards/{fd} -c {fl}".format(fd=folder, fl=flavour)
    command += " -o plots/diffXsecAnalysis/{lep}/{fd}/plotDiffXsecChargeAsymmetry/".format(lep=lepton,fd=folder)
    command += " -t cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
    command += " --lumi-norm 35900.0  -n --palette 57 --hessian --suffix {fit}_{pf} {ptOpt}".format(fit=fit,pf=postfix, ptOpt=ptBinsSetting)
    if fit == "Data":
        command += " --fit-data --invert-ratio --expected-toyfile cards/{fd}/fit/hessian/fitresults_{s}_Asimov_{pf}.root".format(fd=folder,s=seed,pf=postfix)
    if not skipPlot:
        print ""
        print command
        if not dryrun:
            os.system(command)


    print ""
    print "="*30
    print "POSTFIT NUISANCES"
    print ""

    ## POSTFIT NUISANCES
    command = "python w-helicity-13TeV/diffNuisances.py"
    command += " --infile cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
    command += " --outdir plots/diffXsecAnalysis/{lep}/{fd}/diffNuisances/{pf}/".format(lep=lepton,fd=folder, pf=postfix)
    command += " -a --format html --type hessian  --suffix  {fit} ".format(fit=fit)
    for poi in diffNuisances_pois:
        tmpcommand = command + " --pois '{poi}' ".format(poi=poi)
        if not skipDiffNuis:
            print ""    
            print tmpcommand
            if not dryrun:
                os.system(tmpcommand)

    print ""
    print "="*30
    print "POSTFIT PLOTS"
    print ""

    ## POSTFIT PLOTS
    command = "python w-helicity-13TeV/postFitPlots_xsec.py"
    command += " cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root cards/{fd}/ ".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
    command += " -o plots/diffXsecAnalysis/{lep}/{fd}/postFitPlots_xsec/{pf}/".format(lep=lepton,fd=folder, pf=postfix)
    command += " --no2Dplot-signal-bin -n  ".format(fit=fit)  # -n
    if fit == "Data":
        if not skipPostfit:
            print ""    
            print command
            if not dryrun:
                os.system(command)


    print ""
    print "="*30
    print "CORRELATION MATRIX"
    print ""

    ## ADD CORRELATION
    command = "python w-helicity-13TeV/subMatrix.py cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
    command += " --outdir plots/diffXsecAnalysis/{lep}/{fd}/subMatrix/{pf}".format(lep=lepton,fd=folder, pf=postfix)
    command += " --nContours 51 --type hessian --suffix  {fit}".format(fit=fit)
    for poiRegexp in correlationSigRegexp:
        for nuisRegexp in correlationNuisRegexp:
            tmpcommand = command + " --params '{nuis},{poi}' ".format(nuis=correlationNuisRegexp[nuisRegexp],poi=correlationSigRegexp[poiRegexp])
            title = "correlation matrix for {t}".format(t=correlationMatrixTitle[nuisRegexp])
            tmpcommand += " --parNameCanvas {nuis}_{poi} --title '{t}' ".format(nuis=nuisRegexp, poi=poiRegexp, t=title)
            if not skipCorr:
                print ""
                print tmpcommand
                if not dryrun:
                    os.system(tmpcommand)

    print ""
    print "="*30
    print "IMPACTS"
    print ""

    ## IMPACTS
    command = "python w-helicity-13TeV/impactPlots.py cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
    command += " -o plots/diffXsecAnalysis/{lep}/{fd}/impactPlots/{pf}/  --suffix {fit} {pt} ".format(lep=lepton,fd=folder,fit=fit,pf=postfix, pt=ptMinForImpacts)
    command += " --abs-value --nContours 51 --margin '0.16,0.15,0.05,0.25' --canvasSize '1500,1200' --splitOutByTarget "
    # --palette 70 --invertPalette: kDarkBody from light blue to red
    # with following line will make impacts with graphs, not matrix                        
    command += " --etaptbinfile cards/{fd}/binningPtEta.txt ".format(fd=folder)
    if flavour != "lep":
        command += " -c {fl}".format(fl=flavour)
    for nuis in impacts_nuis:
        if nuis == "GROUP":
            varopt = " --nuisgroups '{ng}' ".format(ng=groupnames)
        else:
            varopt = " --nuis '{nuis_regexp}' ".format(nuis_regexp=nuis)

        for target in targets:
            if any(target == x for x in ["etaxsec", "etaxsecnorm", "etaasym"]):
                poi_regexp = ["W.*_ieta_.*"]
            else:
                poi_regexp = ["W.*_ieta_.*_ipt_%d_.*" % i for i in [2, 7, 16] ]

            for poi in poi_regexp:
                tmpcommand = command + " {vopt} --target {t} --pois '{poi_regexp}' ".format(vopt=varopt, t=target, poi_regexp=poi)  
                    
                if not skipImpacts:
                    print ""
                    print tmpcommand
                    if not dryrun:
                        os.system(tmpcommand)

print ""
print "THE END"
print ""            

## example
# python w-helicity-13TeV/plotDiffXsecChargeAsymmetry.py -i cards/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3 -o plots/diffXsec/chargeAsymmetry/electron/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3/ -c el -t cards/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3/fit/data/fitresults_123456789_Data_bbb1_cxs1.root  --suffix Data_bbb1_cxs1 --lumi-norm 35900.0  -n --hessian --palette 55 --fit-data --expected-toyfile cards/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3/fit/hessian/fitresults_123456789_Asimov_bbb1_cxs1.root


##
#python w-helicity-13TeV/diffNuisances.py --infile cards/diffXsec_el_2018_12_18_onlyBkg_pt2GeV_last3GeV_eta0p2From2p0/fit/data/fitresults_123456789_Data_allSyst_bbb1_cxs1.root --type hessian --pois ".*FakesEta.*" --outdir plots/diffXsec/diffNuisances/diffXsec_el_2018_12_18_onlyBkg_pt2GeV_last3GeV_eta0p2From2p0/Data_allSyst_bbb1_cxs1/ -a --format html

