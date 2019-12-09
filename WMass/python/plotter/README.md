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

### making friends
First to an interactive test with few events:
```
python postproc_batch.py <mainTreeDir> <friendsDir> --friend --log <logDir> -d <dastasetName> -c 0 -N 20000
```
Then flood condor (4hrs / 250k events)
```
python postproc_batch.py <mainTreeDir> <friendsDir> --friend --log <logDir> -N 250000  --submit  --runtime 240
```

### resubmit friends
```
./scripts/friendChunkResub.py <friendsDir> <mainTreeDir> --run-checker -N 250000
```
Then merge. First create a dir where to put the chunks. Make on the same filesystem where the chunks are to make the mv fast (no copy of the files)
The scripts expects the directory "Chunks", so link it
```
mkdir <mainTreeDir>/friends_chunks
ln -s <mainTreeDir>/friends_chunks Chunks
./scripts/friendChunkAdd.sh -c <mainTreeDir>/friends
```
then if everything went ok you can safely remove the chunks in "<mainTreeDir>/friends_chunks".
Then do a final check:
```
python scripts/checkMergedFriends.py <mainTreeDir>
```
