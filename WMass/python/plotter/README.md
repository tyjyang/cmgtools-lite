# README


## _running the cards if you need to_

### run the cards with the following command

```
python w-mass-13TeV/wmass_mu/make_cards_mu.py --suffix SUFFIX --syst (--dry-run)
```
by default it:
- decorrelates the QCD scale systematics by charge and in bins of W-pT
- reweights the boson-pT for Ws and Zs

### once you think all have run, check with the DatacardChecker.py script:
```
python w-mass-13TeV/DatacardsChecker.py -c cards/wmass_20XX_YY_ZZ_SUFFIX/ -g 5
```

this will put a new condor\_submit file into a retry\_X sub-directory. it will not automatically
submit the file though, so you have to do that. you can also change the grouping of the jobs for
this step


### FOR W-LIKE:

trees with all relevant gen-level variables for W and Z:
```
/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/ntuplesRecoil/TREES_prefiring_muons_fullTrees/
```
