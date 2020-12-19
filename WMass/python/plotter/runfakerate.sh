#! /bin/bash
# this script is a wrapper for the commands used to compute and pack the fake rate
#
plotterPath="${CMSSW_BASE}/src/CMGTools/WMass/python/plotter"
treePath="/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/postNANO/dec2020/"
dataPeriod="all" # all, preVFP, postVFP
#
######################
# options to set
######################
#--------------------------
dryRun="n"
reweightZpt="n" # use W and Z with reweighted pt
useSignedEta="y" # distinguish bins of positive and negative rapidity (if passing binning with just positive values below, it will just skip the negative, so you are actually using half statistics)
charge=""  # "p", "n", or "" for positive, negative or both leptons
useLeptonScaleFactors="y" # to use weight for lepton scale factors (since they are obtained in a different phase space, one should do it with and without and compare)
#--------------------------
ptDefinition="pt"  # pt variable in w-mass-13TeV/make_fake_rates_xvars.txt 
#-------------------------
#today=`date +"%d_%m_%Y"`
#outdir="fr_${today}_eta_${ptDefinition}_mT40_${lumi/./p}fb_signedEta_jetPt30"
outdir="fr_eta_${ptDefinition}_${dataPeriod}"
lumi="35.9" 
if [[ "${dataPeriod}" == "preVFP" ]]; then
    lumi="19.3" 
elif [[ "${dataPeriod}" == "postVFP" ]]; then
    lumi="16.6" 
else:
    lumi="35.9" 
fi
######################
######################
# additional options to be passed to w-mass-13TeV/make_fake_rates_data.py
# can pass a new cut as you would do with mcPlots.py 
# for now passing options to run on nanoAOD (as implemented in mcAnalysis.py)
#addOption=" -A alwaystrue pfmet 'MET_T1_pt<30' "
addOption=" --nanoaod-tree --max-genWeight-procs 'W.*|Z.*' '50000.0' --clip-genWeight-toMax "

if [[ "${useLeptonScaleFactors}" != "y" ]]; then
    outdir="${outdir}_noLepSF"
fi

chargeCmd=""
if [[ "${charge}" == "p" ]]; then
    chargeCmd=" --charge \"p\" "
    outdir="${outdir}_plus"
elif [[ "${charge}" == "n" ]]; then
    chargeCmd=" --charge \"n\" "
    outdir="${outdir}_minus"
fi

outdir="${plotterPath}/plots/fake-rate/wmassUL2016/${outdir}/"
echo "Creating ${outdir}"
mkdir -p ${outdir}

cmdComputeFR="python ${plotterPath}/w-mass-13TeV/make_fake_rates_data.py --qcdmc  --pt ${ptDefinition} --lumi ${lumi} --tree-path ${treePath} --outdir ${outdir} ${chargeCmd} "
    
if [[ "${useLeptonScaleFactors}" != "y" ]]; then
    cmdComputeFR="${cmdComputeFR} --no-scaleFactors"
fi

if [[ "${useSignedEta}" == "y" ]]; then
    cmdComputeFR="${cmdComputeFR} --useSignedEta "
fi

if [[ "${reweightZpt}" == "y" ]]; then
    cmdComputeFR="${cmdComputeFR} --reweightZpt "
fi

if [[ "X${addOption}" != "X" ]]; then
    cmdComputeFR="${cmdComputeFR} --addOpts \"${addOption}\" "
fi

if [[ "${dryRun}" == "y" ]]; then
    cmdComputeFR="${cmdComputeFR} --dry-run "
    echo "Running this command:"
    echo ""
    echo "${cmdComputeFR} > commands4fakeRate.sh"
    echo ""
    echo "${cmdComputeFR} > commands4fakeRate.sh" | bash
    echo "The commands used for fake-rate are stored in commands4fakeRate.sh"
    echo "Use the following command to really run things"
    echo ""
    echo "cat commands4fakeRate.sh | bash"  # here we really run the commands saved in commands4fakeRate.sh
    echo ""
else
    echo "Running this command:"
    echo ""
    echo "${cmdComputeFR}"
    echo ""
    echo "${cmdComputeFR}" | bash    
fi

