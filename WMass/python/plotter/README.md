# README

## _running the cards if you need to_

### Make histograms for the fit

This part uses the same plotting functionalities provided by **mcPlots.py**.
First, one needs all the proper input txt files, in particular the one with the plot definitions.

#### Prepare txt file with plots to produce

Use the following command (for Wlike analysis)
```
python w-mass-13TeV/testingNano/cfg/makePlotsCFG_systTH3.py -o w-mass-13TeV/testingNano/cfg/plots_wlike_sysTH3.txt -a wlike
```
or the following for Wmass analysis (default, option _-a wmass_ can be skipped)
```
python w-mass-13TeV/testingNano/cfg/makePlotsCFG_systTH3.py -o w-mass-13TeV/testingNano/cfg/plots_wmass_sysTH3.txt -a wmass
```
This is set up to make a TH2 (eta-pt) for nominal histogram, and some TH3 (eta-pt and systematic index) for systematic uncertainties. It currently defines histograms for:
- QCD scale systematics, unbinned and in bins of boson pT (binned ones on signal and unbinned ones on anti-signal, e.g. Z for Wmass)
- PDFs and alphaS
- Efficiency statistical uncertainty (a nuisance for each bin of the tag-and-probe), although the current scale factor histograms might actually store the stat+syst uncertainty (**TO BE FIXED**)
- mass weights (actually the branch also include other EW stuff) 
Other systematics can be added in a similar way

#### Produce the histograms

Basically use **mcPlots.py** adding option _--skipPlot_ to avoid making plots (as we have many TH3), although we may at least let it plot the nominal ones for cross check. This means all options that are specific to plotting are not needed.

The command below is for positive charge (at some point both charges might be done together, using TH4 or something). For negative charge the command is the same, but you have to change some options accordingly.

Some notes:
- the command currently is set up to make all the plots in the plot txt file (i.e. _--sP ".*"_)
- some input files (e.g. the mca.*txt) have the same content for wmass and wlike, despite their different names
- can add _--allow-negative-results_ or _--neglist "<regular_expression>"_ to avoid cropping negative bins to 0 (can be done in next steps if needed, or by combinetf directly)
- can add _--updateRootFile_ to remake only some histograms (for instance selected with _--sP \<regexp\>_) and update an existing file
- can add _--filter-proc-files ".\*" ".\*\_1.root"_ to run on only some files (1 per process in this case), e.g. for testing
- for Wlike, need to select only odd (even) events for positive (negative) charge
- we are currently removing the mT cut (with _-X mtl1pf40_), until we agree on the MET to use (and also because for fakes we may use a simultaneous fit including the low mT region
- luminosity value passed to option _-l_ is used both as event weight and as the number appearing in plots. Thus, if the normalization already appears elsewhere (e.g. in the MCA file multiplying the cross section for preVFP and postVFP) the value passed to _-l_  should be offset by adding _1.0/luminosity_ as event weight with option _-W_ (or change the current behaviour)

Wlike (using odd events for charge plus)
```
python mcPlots.py w-mass-13TeV/testingNano/cfg/mca-wlike.txt w-mass-13TeV/testingNano/cfg/test/cuts_wlike.txt w-mass-13TeV/testingNano/cfg/plots_wlike_sysTH3.txt -P /data/shared/originalNANO/ -p "data,Zmumu,Wmunu,Ztautau,Wtaunu" --pg "data := data_preVFP,data_postVFP" --pg "Wmunu := Wmunu_plus_preVFP,Wmunu_plus_postVFP,Wmunu_minus_preVFP,Wmunu_minus_postVFP" --pg "Wtaunu := Wtaunu_plus_preVFP,Wtaunu_plus_postVFP,Wtaunu_minus_preVFP,Wtaunu_minus_postVFP" --pg "Zmumu := Zmumu_preVFP,Zmumu_postVFP" --pg "Ztautau := Ztautau_preVFP,Ztautau_postVFP" --sP ".*" --nanoaod-tree -X mtl1pf40  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_wlike.txt  --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 -f -l 36.3 -W "_get_fullMuonSF(Muon_pt[goodMuonsCharge][0],Muon_eta[goodMuonsCharge][0],Muon_charge[goodMuonsCharge][0],Muon_pt[goodMuonsOther][0],Muon_eta[goodMuonsOther][0],eraVFP)*_get_MuonPrefiringSF(Muon_eta,Muon_pt,Muon_looseId,eraVFP)*puw_2016UL_era(Pileup_nTrueInt,eraVFP)" --out cards/wlike/plus/wlike.root --skipPlot -A trigger oddevents "isOddEvent(event)"
```

Wmass (example for charge plus)
```
python mcPlots.py w-mass-13TeV/testingNano/cfg/mca-wmass.txt w-mass-13TeV/testingNano/cfg/test/cuts_wmass.txt w-mass-13TeV/testingNano/cfg/plots_wmass_sysTH3.txt -P /data/shared/originalNANO/ -p "data,Wmunu_plus,Wmunu_minus,Zmumu,Ztautau,Wtaunu_plus,Wtaunu_minus" --pg "data := data_preVFP,data_postVFP" --pg "Wmunu_plus := Wmunu_plus_preVFP,Wmunu_plus_postVFP" --pg "Wmunu_minus := Wmunu_minus_preVFP,Wmunu_minus_postVFP" --pg "Wtaunu_plus := Wtaunu_plus_preVFP,Wtaunu_plus_postVFP" --pg "Zmumu := Zmumu_preVFP,Zmumu_postVFP" --pg "Wtaunu_minus := Wtaunu_minus_preVFP,Wtaunu_minus_postVFP" --pg "Ztautau := Ztautau_preVFP,Ztautau_postVFP" --sP ".*" --nanoaod-tree -X mtl1pf40  --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_wlike.txt  --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 -f -l 36.3 -W "_get_fullMuonSF(Muon_pt[goodMuonsCharge][0],Muon_eta[goodMuonsCharge][0],Muon_charge[goodMuonsCharge][0],-1,-1,eraVFP)*_get_MuonPrefiringSF(Muon_eta,Muon_pt,Muon_looseId,eraVFP)*puw_2016UL_era(Pileup_nTrueInt,eraVFP)" --out cards/wmass/plus/wmass.root --skipPlot -A onemuon chargeplus "Muon_charge[goodMuons][0] > 0"
```

**IMPORTANT**: We can actually now make all histograms including fakes alltogether, which is currently accomplished by actually computing yields simultaneously in the four regions defined by mT/isolation. Check the section __Making QCD template for Wmass__ for details.

#### Unpack the histograms into TH2 (eta-pt) for combinetf

Currently using **makeHistogramsWMass.py** as an independent script, but it might be merged inside the card maker script. It takes as input the root file produced in previous step, and also produce the alternate histograms for some systematics by mirroring the alternative template with respect to nominal one (e.g. for PDFs, which do not have Up and Down by default). It also write the histograms with a proper name following the conventions used by combinetf.

Some nuisances might be decorrelated by charge, in which case the histograms should be named accordingly. This can be done automatically with _--decorrelate-by-charge_ as shown below (Plus or Minus is appended to the original syst name)

Wlike
```
python makeHistogramsWMass.py -i cards/wlike/plus/wlike.root -o cards/wlike/Zmmu_plus_shapes.root -c plus --wlike --decorrelate-by-charge ".*effStatTnP|.*muR\d+|.*muF\d+" [--crop-negative-bin]
```

Wmass
```
python makeHistogramsWMass.py -i cards/wmass/plus/wmass.root -o cards/wmass/Wmunu_plus_shapes.root -c plus --decorrelate-by-charge ".*effStatTnP|.*muR\d+|.*muF\d+" [--crop-negative-bin]
```

Once you have these final root files, you can inspect their content with **printRootFileContent.py**, it also has some options to exclude or select some histogram names using regular expressions.

E.g. the following will print histograms whose name contains _pdfs_ or _muRmuF3_, but excluding those for process Wtaunu
```
python printRootFileContent.py cards/wmass/Wmunu_plus_shapes.root -r ".*pdf.*|.*muRmuF3" -x "x_Wtaunu.*"
```

### Check some alternative histograms for some systematics and processes

Example command
```
python w-mass-13TeV/makeSystRatios.py cards/wmass/Wmunu_plus_shapes.root plots/testNanoAOD/WmassPlots/fakeRateRegion_postVFP_plus_systTH3/postprocessing/testHistoFits/ -p "Zmumu,Wmunu_plus,Top,data_fakes" -s ".*muonL1Prefire(1|3)Up|.*pdf12|.*muRmuF9PlusUp|.*muRUp|.*effStatTnP(1|112|500)PlusUp|.*lumi.*"
```

### Make the cards and run the fit with the following command

**IMPORTANT**: The command to run the fit needs to be used from a cmssw-cc7 environment, within a release (e.g. __CMSSW_10_6_19_patch2__ should work) and running cmsenv beforehand. Currently, __w-mass-13TeV/cardMaker.py__ is run using python (not python3), even though the script is written using python3 (but it seems to be fine, except that printed messages are displayed as ntuples).

The following command is just a simplified example to produce the datacard for a single charge (can use _-c plus,minus_ to make cards for both charges). The script allows one to make datacards for a single charge or both, possibly combining them, and to execute the commands to actually run the fit. It has some options to customize the datacard content (for instance, to exclude some nuisances on the fly) and to configure the text2hdf5 or combinetf commands. By default the fit is not run, to do it some options are needed (under testing at the moment)

Wlike
```
python w-mass-13TeV/cardMaker.py -i cards/wlike/  -f mu -c plus --wlike
```

Wmass (default)
```
python w-mass-13TeV/cardMaker.py -i cards/wmass/  -f mu -c plus
```

The systematics to be applied on a specific process are configured in __w-mass-13TeV/wmass_mu/systsFit.txt__ (for wmass) and __w-mass-13TeV/wlike_mu/systsFit.txt__ (for wlike). One may also propagate all theoretical and experimental uncertainties on fakes (assuming the corresponding histograms exist) by manipulating these files.

Among the main general options:
- _--comb_: combine the datacards for the two charges
- _--exclude-nuisances_: exclude some nuisances when writing the card, using regular expressions
- _--keep-nuisances_: keep these nuisances (overriding those excluded by _--exclude-nuisances_)
- _--mass-nuis_: use only this nuisance parameters for mass shift, neglecting all the others

Some options to customize fit (check the script for more)
- _--fit-single-charge_: fit single charge
- _--postfix_: enable the same option of text2hdf5.py, adding a postfix to output file (to distinguish different files without having to change output folder). This is propagated automatically to combinetf.py


Examples for some realistic fits (Wmass)
- fit each single charge independently, freezing POIs (to measure mW), selecting mass shift of 100 MeV, and skipping fit to data (thus only doing Asimov)
```
python w-mass-13TeV/cardMaker.py -i cards/wmass/  -f mu -c "plus,minus" --fit-single-charge --freezePOIs --mass-nuis massShift100MeV --impacts-mW --skip-fit-data --all-proc-background
```
- fit combination
```
python w-mass-13TeV/cardMaker.py -i cards/wmass/  -f mu -c "plus,minus" --comb --freezePOIs --mass-nuis massShift100MeV --impacts-mW --skip-fit-data --all-proc-background
```
### Plot impacts on mW

```
python w-mass-13TeV/makeImpactsOnMW.py cards/wmass/fit/hessian/fitresults_123456789_Asimov_bbb1_cxs1.root  -o plots/testNanoAOD/fits/test_impacts/ --nuisgroups ALL --prefitUncertainty 100 --scaleToMeV --showTotal
```

## Prepare scale factors for the analysis

Once the scale factors are computed with the tag-and-probe and the output root file is available somewhere, you can plot them all and make their product using the following command

```
python w-mass-13TeV/plotSF.py testMuonSF/2021-05-31_allSFs_nodz_dxybs.root plots/testNanoAOD/testSF/SFeta0p1_31May2021_nodz_dxybs/globalAndPerEra/  -e "BtoF,GtoH,B,C,D,E,F,G,H" -n 'trigger,idip,iso,antiiso,isonotrig,antiisonotrig' --makePreOverPost
```
This will plot the SF passed to _-n_ for pre and postVPF eras (and possibly single eras if available). It also plots the absolute and relative uncertainties to provide a full picture of how they look like. Option _--makePreOverPost_ will make the script call  _w-mass-13TeV/makeEffRatioPrePostVFP.py_ to compute the scale factors for preVFP/postVFP in data and MC.
The script also makes and plots the products of the SF (e.g. for trigger with either charge, isolation, and idip), whose usage is more convenient at analysis level. Tracking and reco SF are compatible with 1 (efficiencies are close to 100%, so better to neglect them)

The previous command assumes that the input root file contains histograms for the preVFP or postVFP era. For tests, we made scale factors also split by era. Currently, they can be manipulated using the following script (not worth merging it with the other one, as we don't expect to use scale factors for each single era).
The following command will make a lot of plots, including data/data and MC/MC scale factors, comparison of inclusive efficiencies with luminosity-weighted average of the different eras, and many more things. 
```
python w-mass-13TeV/compareSFperEra.py testMuonSF/2021-05-31_allSFs_nodz_dxybs.root plots/testNanoAOD/testSF/SFeta0p1_31May2021_nodz_dxybs/ -e BtoF
```
A subfolder named "compareSFperEra" will be automatically created inside the provided output folder, and all plots will be created there.

## Making QCD template for Wmass

### Prepare MCA file
The idea is to make histograms in four regions defined by isolation and transverse mass, where the low mT region also has the requirement of at least one jet not overlapping with the muon.
Instead of using the general formula to apply the fake rate method, which relies on both the fake and prompt rates, we can use directly the simplified formula given by f/(1-f) * N(QCD)_highIso, where f denotes the fake rate and N(QCD)_highIso is the number of QCD background events in the kinematic region with high muon isolation and high mT (so orthogonal to the signal region).
Since f=N_passIso/N_tot (evaluated in the fake rate computation region, with low mT and 1 jet), the factor f/(1-f) reduces to N_passIso/N_failIso. At the same time, N(QCD)_highIso can be estimated as data-MC, and we can then directly multiply this template to obtain the QCD template in the signal region. 
To make these histograms, first you need to use the following command to prepare the MCA file
```
python w-mass-13TeV/testingNano/cfg/makePlotsCFG_systTH3.py -o w-mass-13TeV/testingNano/cfg/plots_fakerate_systTH3.txt --a wmass -b 48,-2.4,2.4,116,26,142 --ptVar "Muon_pt[goodMuons][0]+29.0*regionIsoMt(Muon_pfRelIso04_all[goodMuons][0]<0.15,transverseMass<40)"
```
This also makes the histograms for the systematic variations on the prompt lepton templates.

### Produce all the histograms for data and MC processes.

The MC histograms will make the prompt component to be subtracted from data.
One can make all the histograms with the following command (separately for preVFP and postVFP era for now, they may be gathered at a later stage).
One can make the histograms separately for each charge by running twice with option __-c (plus|minus)__. The charge dependence may arise from the prompt subtraction, which is largely dominated by W.
```
python runFakeRate.py -o plots/testNanoAOD/WmassPlots/ -e postVFP --variables ".*" --plot-file "plots_fakerate_systTH3.txt" --options " --skipPlot " -c plus -s
python runFakeRate.py -o plots/testNanoAOD/WmassPlots/ -e postVFP --variables ".*" --plot-file "plots_fakerate_systTH3.txt" --options " --skipPlot " -c minus -s
```
One can also make the histograms (and thus the analysis) for a specific sub era adding option __--sub-era D__ (for Run2016D). The proper scale factors, pileup, and muon prefiring weights will be customized accordingly inside the script (please always check that the correct commands are built before really running, you can use option __-d__ to print without running).

Once the histograms are available, one has to manipulate them to get the QCD prediction, according to the formula described above. This can be done using this command, where the input file is just the output of the previous command
```
python w-mass-13TeV/plotFakesTemplate.py plots/testNanoAOD/WmassPlots/fakeRateRegion_postVFP_plus_systTH3/plots_fakerate_systTH3.root plots/testNanoAOD/WmassPlots/fakeRateRegion_postVFP_plus_systTH3/postprocessing/ -b "29,26,55"
python w-mass-13TeV/plotFakesTemplate.py plots/testNanoAOD/WmassPlots/fakeRateRegion_postVFP_minus_systTH3/plots_fakerate_systTH3.root plots/testNanoAOD/WmassPlots/fakeRateRegion_postVFP_minus_systTH3/postprocessing/ -b "29,26,55"
```
Note that the pt binning passed with __-b__ must be consistent with the one used before (which had a range 4 times larger because including the 4 iso/mT regions).

The command above will also propagate a default 1.2% uncertainty for luminosity on the prompt component, which is also propagated to the fake templates through the subtractions of prompt events from data.

At this point one can continue by running **makeHistogramsWMass.py** as explained in section __Unpack the histograms into TH2 (eta-pt) for combinetf_

#### Make plots in iso/mT regions for checks

The same __runFakeRate.py__ script can be used without option __-s__ and removing __--skipPlot__ to plot histograms of given variables in all the four iso/mT regions. An example command is the following (default variables should be already adequate):
```
python runFakeRate.py -e postVFP --plot-file "plots_fakerate.txt" [-c plus|minus]
```

### Vertex study
Prepare the commands to make histograms for this study with the utility script __runVertexStudy.py__
```
python runVertexStudy.py
```
Once the histograms are available, the actual plots can be made with a command like the following (see other options therein)
```
python w-mass-13TeV/makeVertexStudy.py plots/testNanoAOD/vertexStudy/Wboson_noPUorPrefire/ -p Wmunu_plus_preVFP -v wpt --hname "mueta_dzGenRecoVtx" -e -2.39 -2.21 -e 2.21 2.39 --postfix endcap2p2to2p4
```

### MC truth efficiency
Prepare the commands to make histograms for this study with the utility script __runMCtruthEfficiency.py__ (you may need to customize the command inside, or the eras you want to process)
```
python runMCtruthEfficiency.py
```
Then, you can make all sort of plots with a command like the following (check available options to customize binning, ranges, or plots to make)
```
python w-mass-13TeV/makeMCtruthEffStudy.py plots/testNanoAOD/MCtruthEfficiency/W_perEra_noRecoAccept_finePt/plus/ --palette 87 --plot-ratio-era H --tnp-file testMuonSF/2021-05-31_allSFs_nodz_dxybs.root --rebin-pt 4 -w 'trackerOrGlobal,tracker,standalone,global,accept,idip,trig,trigNoBit,iso,idipANDtrig,idipANDisonotrig,idipANDtrigANDiso,idipANDtrigNoBit'
```


## Details about making plots

### cut file
Currently still using txt file for cuts, as in previous CMGTools versions. Expressions can use standard C++ syntax. Some special character needs to be escaped when they are also used as field separators. For instance, std::abs need to be used as std\\:\\:abs, because ':' is the separator between cut name and expression

Can also add a third field separated from the rest by ';', where a comma-separated list of several options can be passed (format is key=value or just value for bools, the latter case sets the options to True).
**NEW**: can use BEFOREDEFINE=True to make this filter be defined before other RDataframe Defines are declared. This should speed up things a little bit, but of course those cuts should not require any newly defined column, except for some special ones that are internally defined in tree2yield.py (see inside getManyPlotsRaw() function in the TreeToYield class). Typically it can be used on the trigger bits or json file selection.

### plot file
Added option to customize event weight per histogram.
- Use **AddWeight="expression"** (just as it was already possible in MCA file per process) to assign an additional weight for those specific histograms (commas in expression needs to be escaped as '\,'). This can be a constant or a function or even a product of functions.
- use **ReplaceWeight="oldexpr->newexpr"** to replace oldexpr in the nominal weight with newexpr (they are separated by '->'). The format is the same as for AddWeight (e.g. commas need to be escaped). This can be used to replace a part of an expression, e.g. to change name of a function or even its arguments, although it should be used with caution to avoid bugs (regular expressions are not supported, it is simply used within python in the replace() method for strings).

### built-in features of mcPlots.py
By default, the script relies on classes from mcAnalysis and tree2yields. For MC processes, the gen event weight column is defined internally (in mcAnalysis.py). When looping on events with RDF, the filters are created separately for each line in the CUT file, to speed things up a little bit. This splitting is also necessary to print yields for the cutflow, which is the default behaviour now (among other outputs, a *_yields.txt file is created).
The yields printed with the cutflow currently do not take into account gen weights for MC, this will be added at some point.

Inside tree2yield.py, the main function managing the histograms is getManyPlotsRaw() from the TreeToYields class. Here, some RDF columns are created to ease the usage of different processes dynamically. Currently they are:
- **eraVFP**: set to BToF (GToH) for preVFP (postVFP), and BToH otherwise. Whether a process is pre or post VFP depends on its name as defined in the MCA file. The values are enum types defined in ccFiles/defines.h
- **isData**: with similar logic to previous item, it is set to 'Data' when the process name starts with 'data', and MC otherwise. Again, Data and MC are enum types in ccFiles/defines.h