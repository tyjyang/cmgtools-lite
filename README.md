# Short recipe for CMGTools 

For the general recipe, [follow these instructions](https://twiki.cern.ch/twiki/bin/view/CMS/CMGToolsReleasesExperimental).

--------------
```
# setup CMSSW and the base git (do on a lxplus6 to build it)
cmsrel CMSSW_9_4_12
cd CMSSW_9_4_12/src
cmsenv
git cms-init

# add the wmass cmg-cmssw repository to get the Heppy heppy_94X branch
git remote add wmass-central  git@github.com:WMass/cmg-cmssw.git  -f  -t heppy_94X

# configure the sparse checkout, and get the base heppy packages
cp /afs/cern.ch/user/e/emanuele/public/wmass/sparse-checkout_94X_heppy .git/info/sparse-checkout 
git checkout -b heppy_94X wmass-central/heppy_94X

# now get the CMGTools subsystem from the cmgtools-lite repository
git clone -o wmass-central git@github.com:WMass/cmgtools-lite.git  -b wmass94X CMGTools

#compile
cd $CMSSW_BASE/src && scram b -j 8

#run
cd CMGTools/WMass/cfg
heppy Trash run_wmass_cfg.py -N 5000 --option test=1
```
