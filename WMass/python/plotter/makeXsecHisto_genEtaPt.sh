#!/bin/bash

# can run as:
# bsub -q cmscaf1nw -oo /afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/plotter/makeXsecHisto_genEtaPt.log  "source /afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/plotter/makeXsecHisto_genEtaPt.sh"

cd /afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/plotter
eval `scramv1 runtime -sh`

# should first run the script to create the MCA:
#
# ELECTRON
# python w-helicity-13TeV/printMCAforXsec.py -o w-helicity-13TeV/wmass_e/mca-includes/ -n mca-80X-wenu-sigInclCharge_gen_eta_pt_4xsec.txt -c el -l preFSR -w "WJetsToLNu_NLO*"
#
# MUON
# python w-helicity-13TeV/printMCAforXsec.py -o w-helicity-13TeV/wmass_mu/mca-includes/ -n mca-80X-wmunu-sigInclCharge_gen_eta_pt_4xsec.txt -c mu -l preFSR -w "WJetsToLNu_ext2v5_*"  
# do not need all MC, depends on binning

plotterDir="${CMSSW_BASE}/src/CMGTools/WMass/python/plotter/"

# set these carefully
wDir="${plotterDir}w-helicity-13TeV/wmass_mu/"
mcafile="${wDir}mca-includes/mca-80X-wmunu-sigInclCharge_gen_eta_pt_4xsec.txt"
dirName="test/xsection_genAbsEtaPt_dressed_mu_pt0p5_eta0p1_etaGap_yields_v2/"
selectPlot=" --sP gen_ptl1_absetal1_dressed "  # gen_ptl1_absetal1_dressed 


#dirName="test/utilityDistributions/"
outDir="${plotterDir}plots/${dirName}"
#selectPlot=" --sP ptl1_gen_reco,etal1_gen_reco"  # use "" for all
#excludeProcesses="--xp data,Z_LO,W_LO,data_fakes,TauDecaysW,WFlips,Top,DiBosons"
#excludeProcesses=""
#selectProcesses=" -p Wplus_el_central,Wminus_el_central" 
#selectProcesses=" -p Wplus_mu_central,Wminus_mu_central" 
#selectProcesses="" 
ncores="8"  # use 2 or 4 if running on batch, can be 8 otherwise
#lumi=0.001   # N = L*sigma, I want plot of sigma in pb, so I compute N for 1/pb (lumi is defined in 1/fb in mcPlots.py)
lumi=35.9  # to plot absolute yields in 35.9/fb


# do not change the folder here
plotfile="${plotterDir}w-helicity-13TeV/wmass_e/tmp_plots.txt"
cutfile="${plotterDir}w-helicity-13TeV/wmass_e/alwaystrue.txt"
# do not use weights here, this is gen only and the weights induce reco cuts
weight=""

#treepath="/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/trees/"
#treedir="TREES_electrons_1l_V6_TINY/"
#treepath="/eos/user/m/mdunser/w-helicity-13TeV/trees"
#treedir="TREES_latest_1muskim/"
#treepath="/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/ntuplesRecoil/"
#treedir="TREES_SIGNAL_1l_recoil_fullTrees/"
treepath="/afs/cern.ch/work/m/mdunser/public/wmassTrees/"
treedir="SKIMS_muons_latest/"


#addCommand=" --contentAxisTitle 'cross section [pb]'"
addCommand=""
#addCommand=" --max-entries 1000 --drawBox 0,2.5,30,45"
#addCommand=" -X fiducial -X eleKin -X json -A onelep chargeMatch 'LepGood1_mcMatchId*LepGood1_charge!=-24' "

command="python ${plotterDir}mcPlots.py -f -l ${lumi} --s2v --tree treeProducerWMass --obj tree --lspam '#bf{CMS} #it{Preliminary}' --legendWidth 0.25 --legendFontSize 0.05 ${mcafile} ${cutfile} ${plotfile}  -P ${treepath}/${treedir} -F Friends ${treepath}/${treedir}friends/tree_Friend_{cname}.root ${selectPlot} --noCms ${addCommand} -j ${ncores} --pdir ${outDir} ${weight}"
# ${selectProcesses}   ${excludeProcesses}

echo "#-----------------------------------------------------"
echo "#-----------------------------------------------------"


    thiscommand="${command} "
    echo "${thiscommand}"
    #echo "${thiscommand}" | bash
    echo ""
        
echo "#-----------------------------------------------------"

echo ""


