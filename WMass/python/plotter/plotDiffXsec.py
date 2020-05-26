#!/bin/env python

import ROOT, os, sys, re, array

# to run plots from Asimov fit and data. For toys need to adapt this script

doMuElComb = 1
doMuon = 1 # if 0 do electrons, but doMuElComb overrides it if doMuElComb=1, which runs the combination
combineElePt01asBkg = 0
dryrun = 1
skipPreliminary = True # passed as an option to some scripts to print "Preliminary" in plot
skipData = 0
onlyData = 1
corrXsecStat = 1 # default should be 1, i.e. combinetf had option correlate-xsec-stat, else 0

skipInclusivePlot = 1
skipPlot = 1
skipTemplate = 1
skipDiffNuis = 0
skipPostfit = 1  # only for Data
skipCorr = 1
skipCorr1D = 1
skipCorrAll4HEPdata = 1
skipImpacts = 1
skipImpactsAll4HEPdata = 1
skipImpactsEtaPt = 1
skipMuElComparison = 1
#outFolderComparison = "test_nativeMCatNLOxsecW_profileLepScale_cropNegBinNomi_uncorrFSRbyFlav_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_noSplitElePtSystByPt_FSRshapeOnly" # update name here when using skipMuElComparison, or just use postfix

useXsecWptWeights = 0 # to plot the band better to keep the unweighted xsec (so keep 0)
allPtBinsSignal = 1
forceAllptbinsTheoryband = 1 # for electrons when making xsec plots, to use all pt bins to make theory band
#
# some script allow to plot a single charge
plotSingleCharge = 0
singleChargeToPlot = "minus" # "minus"

seed = 123456789

#folder = "diffXsec_el_2019_04_13_newSystAndWtau/"
#folder = "diffXsec_el_2019_05_13_eta0p2widthFrom1p3_last2p1to2p4/"
#folder = "diffXsec_mu_2019_04_09_newSystAndWtau_fixTriSF/"
#folder = "diffXsec_mu_2019_04_28_eta0p2widthFrom1p3_last2p1to2p4/"
#folder = "diffXsec_mu_2019_05_04_etaReco0p1_etaGen0p2from1p3_last2p1to2p4//"
#folder = "diffXsec_mu_2019_05_09_recoEta0p1_recoPt1_genEta0p2from1p3_last2p1to2p4_genPt2/"
#folder = "diffXsec_mu_2019_06_17_zptReweight/"
#folder = "diffXsec_el_2019_06_21_zptReweight/"
#folder = "diffXsec_mu_2019_06_17_zptReweight_chargeUncorrQCDscales_EffStatOnlyStatUncDataMC/"
#folder = "diffXsec_mu_2019_06_17_zptReweight_chargeUncorrQCDscales_fixFSRcharge/"
#folder = "diffXsec_mu_2019_06_17_zptReweight_chargeUncorrQCDscales_unfixedFSRcharge_testBinUncEffStat/"
#folder = "diffXsec_mu_2019_07_12_noSyst/"
#folder = "diffXsec_el_2019_06_21_zptReweight_fixEffStat/"
#folder = "diffXsec_el_2019_07_20_latestScaleFactor_AllIn_IDwithMConlyStat/"
#folder = "diffXsec_el_2019_07_20_latestScaleFactor_AllIn_IDwithMConlyStat_allPtBinsAsSignal/"
#folder = "diffXsec_el_2019_07_28_testPt2GeV/"
#folder = "diffXsec_mu_2019_08_02_testBinnedSFandUnc/"
#folder = "diffXsec_mu_2019_09_19_nativeMCatNLOxsec/"
#folder = "diffXsec_mu_2019_09_19_nativeMCatNLOxsec_1sigBin_4fixedPOI_ptMax45/"
#folder = "diffXsec_mu_2019_09_19_nativeMCatNLOxsec_1sigBin_4fixedPOI/"
#folder = "diffXsec_el_2019_09_22_nativeMCatNLOxsec/"
#folder = "diffXsec_el_2019_09_22_nativeMCatNLOxsec_testPrefirePtLess35/"

folder = ""
folder_mu = "diffXsec_mu_2019_09_19_nativeMCatNLOxsec/"
#folder_mu = "diffXsec_mu_2019_09_19_nativeMCatNLOxsec_1sigBin_4fixedPOI_ptMax45/"
#folder_mu = "diffXsec_mu_2020_04_10_nativeMCatNLOxsec_fixJetPrefire_1sigBin_4fixedPOI/"
folder_el = "diffXsec_el_2019_09_22_nativeMCatNLOxsec/"

flavour = ""

if doMuon:
    folder = folder_mu
    flavour = "mu"
else:
    folder = folder_el
    flavour = "el"

lepton = "electron" if flavour == "el"  else "muon"

if doMuElComb:
    allPtBinsSignal = 1
    folder = "muElCombination_allSig_nativeMCatNLOxsec"
    if combineElePt01asBkg:
        folder = "muElCombination_elePt01bkg_nativeMCatNLOxsec"
    #folder = "muElCombination_1Sept2019"
    skipTemplate = 1
    flavour = "lep"
    lepton = "lepton"
    if plotSingleCharge:
        print"Error: conflicting flags doMuElComb and plotSingleCharge. Abort"
        quit()

postfix_el = "nativeMCatNLOxsecW_profilePtScales_newSmoothUncorrScale_cropNegBinNomi_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_noSplitPtSystByPt_FSRshapeOnly_addInclXsec"
postfix_mu = "nativeMCatNLOxsecW_RochesterCorrUncert_cropNegBinNomi_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_FSRshapeOnly_addInclXsec"
#postfix_mu = "nativeMCatNLOxsecW_RochesterCorrUncert_cropNegBinNomi_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_FSRshapeOnly_addInclXsec_addImpactsOnMw_fixedPOI"

if flavour == "el":
    postfix = postfix_el
else:
    postfix = postfix_mu
if doMuElComb:
    postfix = "combinedLep_allSig_nativeMCatNLOxsec_profileLepScale_cropNegBinNomi_uncorrFSRbyFlav_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_noSplitElePtSystByPt_FSRshapeOnly_addInclXsec" 
    if combineElePt01asBkg:
        postfix = "combinedLep_elePt01bkg_nativeMCatNLOxsec_profileLepScale_outLnN30_cropNegBinNomi_uncorrFSRbyFlav_clipSyst1p3_clipSigSyst1p15"
               
if plotSingleCharge:
    postfix = "_symFSRptScalemW_singleCharge{ch}".format(ch=singleChargeToPlot)

postfix += "_bbb1_cxs%d" % corrXsecStat
postfix_el += "_bbb1_cxs%d" % corrXsecStat
postfix_mu += "_bbb1_cxs%d" % corrXsecStat
#postfix += "_bbb1_cxs0"
#postfix += "_bbb0"


fits = ["Asimov", "Data"]

ptBinsSetting = " --pt-range-bkg 25.9 30.1 --pt-range '30,56' " if (not allPtBinsSignal) else ""  # " --eta-range-bkg 1.39 1.61 "
ptMinForImpacts = " --pt-min-signal 30" if (not allPtBinsSignal) else ""
optTemplate = " --norm-width --draw-selected-etaPt 2.05,35.0 --syst-ratio-range 'template' --palette 57 --do-signal-syst '.*smooth.*scaleSyst.*|.*smooth.*scaleStat.*|.*mW.*|.*fsr.*' "  # --draw-selected-etaPt 0.45,38 --zmin 10 # kLightTemperature=87
ptMaxTemplate = "56"
ptMinTemplate = "30" if (flavour == "el" or not allPtBinsSignal) else "26"

# do not ask Wplus.*_ieta_.*_mu$ to select signal strength rejecting pmasked, because otherwise you must change diffNuisances.py
# currently it uses GetFromHessian with keepGen=True, so _mu$ would create a problem (should implement the possibility to reject a regular expression)
# if you want mu rejecting pmasked do _mu_mu or _el_mu (for electrons _mu works because it doesn't induce ambiguities with the flavour)
allSystNuisances = "CMS_.*|.*smooth.*scale.*|.*TestEffSyst.*|mW|fsr.*|L1Prefire.*|OutOfAccPrefire.*|ErfPar.*|Fakes(Eta|Pt).*[0-9]+(mu|el).*|pdf.*|alphaS|muR.*|muF.*"
diffNuisances_pois = [#"pdf.*|alphaS", 
                      #"muR.*|muF.*", 
                      ##"Fakes(Eta|Pt).*[0-9]+mu.*", 
                      ##"Fakes(Eta|Pt).*[0-9]+el.*",
                      #"Fakes(Eta|Pt).*[0-9]+(mu|el).*", 
                      ##"ErfPar0EffStat.*", 
                      ##"ErfPar1EffStat.*", 
                      ##"ErfPar2EffStat.*", 
                      #"CMS_.*|.*smooth.*scale.*|.*TestEffSyst.*|mW|fsr|L1Prefire.*|OutOfAccPrefire.*", 
                      ##"Wplus.*_ieta_.*_mu",     
                      ##"Wminus.*_ieta_.*_mu",
                      allSystNuisances, 
                      ]

# this is appended to nuis below
#correlationSigRegexp = {"Wplus_ieta6ipt6" : ".*Wplus.*_ieta_6_ipt_6_.*mu"  
#                        }
# need to specify which matrix, by default the one for mu is chosen and the bin label looks like Wplus_el_ieta_1_ipt_0_Wplus_el, without ending mu
#correlationSigRegexp = {"Wplus_ieta6ipt8" : ".*Wplus_.*_ieta_6_ipt_8_.*"
#                        }
correlationSigRegexp = {}

correlationNuisRegexp = {#"allPDF"           : "pdf.*,alphaS", 
                         #"someTheory"       : "^pdf([1-9]|1[0-9]|20)$,alphaS,fsr.*,mW", 
                         #"QCDscales"        : "muR.*|muF.*", 
                         # "muR"              : "^muR[1-9]+", 
                         # "muF"              : "^muF[1-9]+", 
                         # "muRmuF"           : "^muRmuF[1-9]+", 
                         "FakesEtaPtUncorr" : "Fakes(Eta|Pt).*[0-9]+.*",
                         #"CMSsyst"          : "CMS_.*|.*TestEffSyst.*|smooth.*scaleSyst.*",
                         #"pTscales"         : ".*smooth.*scale.*",
                         # "ErfPar0EffStat"   : "ErfPar0EffStat.*",
                         # "ErfPar1EffStat"   : "ErfPar1EffStat.*",
                         # "ErfPar2EffStat"   : "ErfPar2EffStat.*"
                         }

correlationMatrixTitle = {"allPDF"           : "all PDFs + #alpha_{S}", 
                          "someTheory"       : "some PDFs + #alpha_{S} + fsr + m_{W}", 
                          "QCDscales"        : "all QCD scales", 
                          "muR"              : "#mu_{R} QCD scales", 
                          "muF"              : "#mu_{F} QCD scales", 
                          "muRmuF"           : "#mu_{R}#mu_{F} QCD scales", 
                          "FakesEtaPtUncorr" : "fakes #eta normalizations and p_{T} slopes", 
                          "CMSsyst"          : "some experimental systematics",
                          "pTscales"          : "p_{T} scales",
                          "ErfPar0EffStat"   : "signal efficiency (uncorr. par 0)",
                          "ErfPar1EffStat"   : "signal efficiency (uncorr. par 1)",
                          "ErfPar2EffStat"   : "signal efficiency (uncorr. par 2)"
                         }

# for impacts
targets = [#"mu", 
           #"xsec", 
           #"xsecnorm",
           #"etaptasym",
           "etaxsec",
           "etaxsecnorm",
           "etaasym",
           "ptxsec",
           "ptxsecnorm",
           "ptasym"
           ]

# for impacts vs pT-eta
targetsPtEta = [#"mu", 
                "xsec", 
                "xsecnorm",
                "asym",
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
#impacts_nuis = ["muTestEffSyst0","muTestEffSyst1","muTestEffSyst2"] 
#groupnames = 'binByBinStat,stat,pdfs,wmodel,EffStat,scales,alphaS'
#groupnames = 'binByBinStat,stat,luminosity,pdfs,QCDTheo,Fakes,OtherBkg,OtherExp,EffStat,EffSyst,lepScale,QEDTheo'
groupnames = 'EffStat,Fakes,QCDTheo,pdfs,luminosity,stat,binByBinStat'
groupnames_4HEPdata = 'binByBinStat,stat,luminosity,pdfs,QCDTheo,Fakes,OtherBkg,OtherExp,EffStat,EffSyst,lepScale,QEDTheo'
#if flavour == "el" or doMuElComb:
#    groupnames += ',L1Prefire'
#groupnamesEtaPt = groupnames
#groupnamesEtaPt = "EffStat,Fakes,binByBinStat,stat"
groupnamesEtaPt = "EffSyst"

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
    command += " cards/{fd} -o plots/diffXsecAnalysis_new/{lep}/{fd}/templateRolling/{pfx}/ -c {fl}".format(fd=folder, lep=lepton, fl=flavour, pfx=postfix)
    command += " --plot-binned-signal -a diffXsec -C {ch} --pt-range '{ptmin},{ptmax}' ".format(ch=charge, ptmin=ptMinTemplate, ptmax=ptMaxTemplate)
    command += "  {opt} ".format(opt=optTemplate)
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
    print "INCLUSIVE CROSS SECTION AND CHARGE ASYMMETRY"
    print ""

    ## INCLUSIVE CROSS SECTION AND CHARGE ASYMMETRY
    command = "python w-helicity-13TeV/plotInclusiveXsec.py"
    command += " -i cards/{fd} -c {fl}".format(fd=folder, fl=flavour)
    command += " -o plots/diffXsecAnalysis_new/{lep}/{fd}/plotInclusiveXsec/".format(lep=lepton,fd=folder)
    command += " -t cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
    command += " --lumi-norm 35900.0 --hessian --suffix {fit}_{pf} {ptOpt}".format(fit=fit,pf=postfix, ptOpt=ptBinsSetting)
    #if plotSingleCharge:
    #    command += " -C {ch}".format(ch=singleChargeToPlot)
    #if combineElePt01asBkg:
    #    command += " --combineElePt01asBkg "
    if useXsecWptWeights:
        command += " --use-xsec-wpt "
    if skipPreliminary:
        command += " --skipPreliminary "
    if forceAllptbinsTheoryband:
        command += " --force-allptbins-theoryband "
    if fit == "Data":
        command += " --fit-data --expected-toyfile cards/{fd}/fit/hessian/fitresults_{s}_Asimov_{pf}.root".format(fd=folder,s=seed,pf=postfix)
        # --invert-ratio
    if not skipInclusivePlot:
        print ""
        print command
        if not dryrun:
            os.system(command)

    ###########
    print ""
    print "="*30
    print "DIFFERENTIAL CROSS SECTION AND CHARGE ASYMMETRY"
    print ""

    ## DIFFERENTIAL CROSS SECTION AND CHARGE ASYMMETRY
    command = "python w-helicity-13TeV/plotDiffXsecChargeAsymmetry.py"
    command += " -i cards/{fd} -c {fl}".format(fd=folder, fl=flavour)
    command += " -o plots/diffXsecAnalysis_new/{lep}/{fd}/plotDiffXsecChargeAsymmetry/".format(lep=lepton,fd=folder)
    command += " -t cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
    command += " --lumi-norm 35900.0  -n --palette -1 --hessian --suffix {fit}_{pf} {ptOpt}".format(fit=fit,pf=postfix, ptOpt=ptBinsSetting)
    if plotSingleCharge:
        command += " -C {ch}".format(ch=singleChargeToPlot)
    if combineElePt01asBkg:
        command += " --combineElePt01asBkg "
    if useXsecWptWeights:
        command += " --use-xsec-wpt "
    if skipPreliminary:
        command += " --skipPreliminary "
    if forceAllptbinsTheoryband:
        command += " --force-allptbins-theoryband "
    if fit == "Data":
        command += " --fit-data --expected-toyfile cards/{fd}/fit/hessian/fitresults_{s}_Asimov_{pf}.root".format(fd=folder,s=seed,pf=postfix)
        # --invert-ratio
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
    command += " --outdir plots/diffXsecAnalysis_new/{lep}/{fd}/diffNuisances/{pf}/".format(lep=lepton,fd=folder, pf=postfix)
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
    command += " -o plots/diffXsecAnalysis_new/{lep}/{fd}/postFitPlots_xsec/{pf}/".format(lep=lepton,fd=folder, pf=postfix)
    command += " --no2Dplot-signal-bin -n ".format(fit=fit)  # -n
    if fit == "Data":
        if not skipPostfit:
            if doMuElComb:
                for flav in ["mu", "el"]:
                    command2 = command + " --plot-flavour-from-combination {fl} ".format(fl=flav)
                    print ""    
                    print command2
                    if not dryrun:
                        os.system(command2)
            else:
                print ""    
                print command
                if not dryrun:
                    os.system(command)


    print ""
    print "="*30
    print "CORRELATION MATRIX"
    print ""

    ## ADD CORRELATION
    for matrixType in ["channelpmaskedexp", "channelpmaskedexpnorm"]:
        command = "python w-helicity-13TeV/subMatrix.py cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
        command += " --outdir plots/diffXsecAnalysis_new/{lep}/{fd}/subMatrix/{pf}/{mt}/".format(lep=lepton,fd=folder, pf=postfix,mt=matrixType)
        command += " --matrix-type {mt}  --nContours 51 --type hessian --suffix  {fit}".format(fit=fit,mt=matrixType)
        #for poiRegexp in correlationSigRegexp:
        for nuisRegexp in correlationNuisRegexp:
            #tmpcommand = command + " --params '{nuis},{poi}' ".format(nuis=correlationNuisRegexp[nuisRegexp],poi=correlationSigRegexp[poiRegexp])
            #title = "correlation matrix for {t}".format(t=correlationMatrixTitle[nuisRegexp])
            #tmpcommand += " --parNameCanvas {nuis}_{poi} --title '{t}' ".format(nuis=nuisRegexp, poi=poiRegexp, t=title)
            tmpcommand = command + " --params '{nuis}' ".format(nuis=correlationNuisRegexp[nuisRegexp])
            title = "0" # set to no title
            tmpcommand += " --parNameCanvas {nuis} --title '{t}' ".format(nuis=nuisRegexp, t=title)
            if not skipCorr:
                print ""
                print tmpcommand
                if not dryrun:
                    os.system(tmpcommand)

    ## ADD CORRELATION for pt-(eta) integrated cross section
    for matrixType in ["channelsumpois", "channelsumpoisnorm"]:
        command = "python w-helicity-13TeV/subMatrix.py cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
        command += " --outdir plots/diffXsecAnalysis_new/{lep}/{fd}/subMatrix/{pf}/{mt}".format(lep=lepton,fd=folder, pf=postfix, mt=matrixType)
        command += " --nContours 51 --type hessian --suffix  {fit} --matrix-type {mt} --margin 0.18,0.11,0.08,0.22".format(fit=fit, mt=matrixType)
        command += " --etaptbinfile cards/{fd}/binningPtEta.txt --canvasSize '1200,1000' ".format(fd=folder)
        for poiRegexp in [".*Wplus_.*ieta.*", ".*Wplus_.*ipt.*", ".*Wminus_.*ieta.*", ".*Wminus_.*ipt.*"]:            
            tmpcommand = command + " --params '{poi},fsr.*,mW,smooth.*scaleSyst.*,.*EffSyst.*,CMS_.*' ".format(poi=poiRegexp)
            title = "correlation matrix: {whichXsec} {var}-integrated cross section for W^{{{chs}}}".format(whichXsec="normalized" if "norm" in matrixType else "absolute",
                                                                                                            var="p_{T}" if "ieta" in poiRegexp else "|#eta|", 
                                                                                                            chs="+" if "plus" in poiRegexp else "-")
            tmpcommand += " --parNameCanvas {poi}_someNuis --title '{t}' -c {fl} ".format(poi=poiRegexp.replace(".*",""), t=title, fl=flavour)
            if not skipCorr1D:
                print ""
                print tmpcommand
                if not dryrun:
                    os.system(tmpcommand)

    ## ADD COVARIANCE with all nuisances to prepare hepdata
    ## for now only saving covariance matrix, not correlation
    filter_matrixType_poiRegexp = {"channelpmaskedexp"     : "W.*pmaskedexp",
                                   "channelpmaskedexpnorm" : "W.*pmaskedexpnorm",
                                   "channelsumpois"        : "W.*sumxsec",
                                   "channelsumpoisnorm"    : "W.*sumxsecnorm",
                                   "channelchargepois"     : "W.*chargeasym",
                                   "channelchargemetapois" : "W.*chargemetaasym",
                                   "channelratiometapois"  : "W.*ratiometaratio",
                               }
    for matrixType in filter_matrixType_poiRegexp.keys():
        poiRegexp=filter_matrixType_poiRegexp[matrixType]
        command = "python w-helicity-13TeV/subMatrix.py cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
        command += " --outdir plots/diffXsecAnalysis_new/{lep}/{fd}/subMatrix_4HEPdata/{pf}/{mt}".format(lep=lepton,fd=folder, pf=postfix, mt=matrixType)
        command += " --nContours 51 --type hessian --suffix  {fit} --matrix-type {mt} --which-matrix covariance --divide-covariance-by-bin-area --margin 0.18,0.11,0.08,0.22 ".format(fit=fit, mt=matrixType)
        command += " --etaptbinfile cards/{fd}/binningPtEta.txt --canvasSize '1200,1000' ".format(fd=folder)
        tmpcommand = command + " --params '{poiRegexp}' --show-all-nuisances ".format(poiRegexp=poiRegexp)
        tmpcommand += " --parNameCanvas AllPoisAndNuis_{poi} -c {fl} ".format(poi=poiRegexp.replace(".*",""), fl=flavour)
        if not skipCorrAll4HEPdata:
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
    command += " -o plots/diffXsecAnalysis_new/{lep}/{fd}/impactPlots/{pf}/  --suffix {fit} ".format(lep=lepton,fd=folder,fit=fit,pf=postfix)
    command += " {pt} ".format(pt=ptMinForImpacts)
    command += " --nContours 51 --margin '0.16,0.2,0.1,0.12' --canvasSize '1200,1000' --splitOutByTarget "
    # --palette 70 --invertPalette: kDarkBody from light blue to red
    # with following line will make impacts with graphs, not matrix                        
    command += " --etaptbinfile cards/{fd}/binningPtEta.txt ".format(fd=folder)
    if flavour != "lep":
        command += " -c {fl}".format(fl=flavour)
    if skipPreliminary:
        command += " --skipPreliminary "
    for nuis in impacts_nuis:

        for target in targets:
            if any(target == x for x in ["etaxsec", "etaxsecnorm", "etaasym"]):
                poi_regexp = ["W.*_ieta_.*"]
            elif any(target == x for x in ["ptxsec", "ptxsecnorm", "ptasym"]):
                poi_regexp = ["W.*_ipt_.*"]
            else:
                poi_regexp = ["W.*_ieta_.*_ipt_%d_.*" % i for i in  [2, 3, 7, 16] ]

            if nuis == "GROUP":
                if all(target != x for x in ["etaxsec", "ptxsec"]):
                    varopt = " --nuisgroups '{ng}' ".format(ng=groupnames.replace(',luminosity',''))
                else:
                    varopt = " --nuisgroups '{ng}' ".format(ng=groupnames)
            else:
                varopt = " --nuis '{nuis_regexp}' ".format(nuis_regexp=nuis)
            for poi in poi_regexp:
                tmpcommand = command + " {vopt} --target {t} --pois '{poi_regexp}' ".format(vopt=varopt, t=target, poi_regexp=poi)  
                    
                if not skipImpacts:
                    print ""
                    print tmpcommand
                    if not dryrun:
                        os.system(tmpcommand)

    print ""
    print "="*30
    print "IMPACTS 4 HEPdata"
    print ""

    ## IMPACTS for HEPdata (more lines than in paper)
    command = "python w-helicity-13TeV/impactPlots.py cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
    command += " -o plots/diffXsecAnalysis_new/{lep}/{fd}/impactPlots_4HEPdata/{pf}/  --suffix {fit} ".format(lep=lepton,fd=folder,fit=fit,pf=postfix)
    command += " {pt} ".format(pt=ptMinForImpacts)
    command += " --nContours 51 --margin '0.16,0.2,0.1,0.12' --canvasSize '1200,1000' --splitOutByTarget "
    # --palette 70 --invertPalette: kDarkBody from light blue to red
    # with following line will make impacts with graphs, not matrix                        
    command += " --etaptbinfile cards/{fd}/binningPtEta.txt ".format(fd=folder)
    if flavour != "lep":
        command += " -c {fl}".format(fl=flavour)
    if skipPreliminary:
        command += " --skipPreliminary "
    for nuis in impacts_nuis:

        for target in targets:
            if any(target == x for x in ["etaxsec", "etaxsecnorm", "etaasym"]):
                poi_regexp = ["W.*_ieta_.*"]
            elif any(target == x for x in ["ptxsec", "ptxsecnorm", "ptasym"]):
                poi_regexp = ["W.*_ipt_.*"]
            else:
                poi_regexp = ["W.*_ieta_.*_ipt_%d_.*" % i for i in  [2, 3, 7, 16] ]

            if nuis == "GROUP":
                varopt = " --nuisgroups '{ng}' ".format(ng=groupnames_4HEPdata)
            else:
                varopt = " --nuis '{nuis_regexp}' ".format(nuis_regexp=nuis)
            for poi in poi_regexp:
                tmpcommand = command + " {vopt} --target {t} --pois '{poi_regexp}' ".format(vopt=varopt, t=target, poi_regexp=poi)  
                    
                if not skipImpactsAll4HEPdata:
                    print ""
                    print tmpcommand
                    if not dryrun:
                        os.system(tmpcommand)


    print ""
    print "="*30
    print "IMPACTS VS PT-ETA"
    print ""
    ## IMPACTS versus pT-eta
    command = "python w-helicity-13TeV/impactPlots_singleNuis_vsPtEta.py" 
    command += " cards/{fd}/fit/{typedir}/fitresults_{s}_{fit}_{pf}.root".format(fd=folder,typedir=typedir,s=seed,fit=fit,pf=postfix)
    command += " -o plots/diffXsecAnalysis_new/{lep}/{fd}/impactPlots_singleNuis_vsPtEta/{pf}/  --suffix {fit} ".format(lep=lepton,fd=folder,
                                                                                                                        fit=fit,pf=postfix)
    command += " --nContours 51 --margin '0.16,0.2,0.1,0.12' --canvasSize '1500,1200' --splitOutByTarget --palette 109 "
    command += " --etaptbinfile cards/{fd}/binningPtEta.txt ".format(fd=folder)
    if flavour != "lep":
        command += " -c {fl}".format(fl=flavour)
    for nuis in impacts_nuis:
        if nuis == "GROUP":
            varopt = " --abs-value --nuisgroups '{ng}' ".format(ng=groupnamesEtaPt)
        else:
            varopt = " --zrange 'template' --nuis '{nuis_regexp}' ".format(nuis_regexp=nuis)

        for target in targetsPtEta:
            if any(target == x for x in ["etaxsec", "etaxsecnorm", "etaasym", "ptxsec", "ptxsecnorm", "ptasym"]):
                continue
            tmpcommand = command + " {vopt} --target {t}  ".format(vopt=varopt, t=target)                      
            if not skipImpactsEtaPt:
                print ""
                print tmpcommand
                if not dryrun:
                    os.system(tmpcommand)

                    
    print ""
    print "="*30
    print "MU-EL-COMB COMPARISON"
    print ""
    ## comparison with muon/electron/combination
    command = "python w-helicity-13TeV/compareMuElDiffXsec.py" 
    command += " --input-muon plots/diffXsecAnalysis_new/muon/{fm}/plotDiffXsecChargeAsymmetry/hessian_{fit}_{pfm}/plotDiffXsecChargeAsymmetry.root".format(fm=folder_mu,fit=fit,pfm=postfix_mu) 
    command += " --input-electron plots/diffXsecAnalysis_new/electron/{fe}/plotDiffXsecChargeAsymmetry/hessian_{fit}_{pfe}/plotDiffXsecChargeAsymmetry.root".format(fe=folder_el,fit=fit,pfe=postfix_el)
    command += " --input-combination plots/diffXsecAnalysis_new/lepton/{fc}/plotDiffXsecChargeAsymmetry/hessian_{fit}_{pfc}/plotDiffXsecChargeAsymmetry.root".format(fc=folder,fit=fit,pfc=postfix)
    command += " -o plots/diffXsecAnalysis_new/comparisons/{pfc}/".format(pfc=postfix)
    command += " -b cards/{fc}/binningPtEta.txt".format(fc=folder)
    if skipPreliminary:
        command += " --skipPreliminary "
 
    if doMuElComb and not skipMuElComparison:
        print ""
        print command
        if not dryrun:
            os.system(command)


print ""
print "THE END"
print ""            

## example
# python w-helicity-13TeV/plotDiffXsecChargeAsymmetry.py -i cards/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3 -o plots/diffXsec/chargeAsymmetry/electron/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3/ -c el -t cards/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3/fit/data/fitresults_123456789_Data_bbb1_cxs1.root  --suffix Data_bbb1_cxs1 --lumi-norm 35900.0  -n --hessian --palette 55 --fit-data --expected-toyfile cards/diffXsec_el_2018_12_16_onlyBkg_pt1GeV_eta0p2From1p3/fit/hessian/fitresults_123456789_Asimov_bbb1_cxs1.root


##
#python w-helicity-13TeV/diffNuisances.py --infile cards/diffXsec_el_2018_12_18_onlyBkg_pt2GeV_last3GeV_eta0p2From2p0/fit/data/fitresults_123456789_Data_allSyst_bbb1_cxs1.root --type hessian --pois ".*FakesEta.*" --outdir plots/diffXsec/diffNuisances/diffXsec_el_2018_12_18_onlyBkg_pt2GeV_last3GeV_eta0p2From2p0/Data_allSyst_bbb1_cxs1/ -a --format html


##
# python w-helicity-13TeV/compareMuElDiffXsec.py --input-muon plots/diffXsecAnalysis_new/muon/diffXsec_mu_2019_09_19_nativeMCatNLOxsec/plotDiffXsecChargeAsymmetry/hessian_Data_nativeMCatNLOxsecW_RochesterCorrUncert_cropNegBinNomi_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_FSRshapeOnly_bbb1_cxs1/plotDiffXsecChargeAsymmetry.root --input-electron plots/diffXsecAnalysis_new/electron/diffXsec_el_2019_09_22_nativeMCatNLOxsec/plotDiffXsecChargeAsymmetry/hessian_Data_nativeMCatNLOxsecW_profilePtScales_newSmoothUncorrScale_cropNegBinNomi_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_noSplitPtSystByPt_FSRshapeOnly_bbb1_cxs1/plotDiffXsecChargeAsymmetry.root --input-combination plots/diffXsecAnalysis_new/lepton/muElCombination_allSig_nativeMCatNLOxsec/plotDiffXsecChargeAsymmetry/hessian_Data_combinedLep_allSig_nativeMCatNLOxsec_profileLepScale_cropNegBinNomi_uncorrFSRbyFlav_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_noSplitElePtSystByPt_FSRshapeOnly_bbb1_cxs1/plotDiffXsecChargeAsymmetry.root -o plots/diffXsecAnalysis_new/comparisons/test_nativeMCatNLOxsecW_profileLepScale_cropNegBinNomi_uncorrFSRbyFlav_clipSyst1p3_clipSigSyst1p15_clipPtScale1p15_decorrPtScaleSystByEta_noSplitElePtSystByPt_FSRshapeOnly/ -b cards/muElCombination_allSig_nativeMCatNLOxsec/binningPtEta.txt
