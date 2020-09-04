#!/bin/bash
cd /afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/rel_slc7/CMSSW_9_4_12/src/CMGTools/WMass/python/postprocessing
eval `scramv1 runtime -sh` 
python $*
