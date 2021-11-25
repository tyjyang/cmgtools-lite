#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export CMSSW_BASE=$(dirname $(dirname $dir))
plotdir=$dir/WMass/python/plotter/plots

if [ -L $plotdir ] && [ ! -e ${plotdir} ]; then
    echo "Removing spurious link $plotdir"
    echo "It will be recreated"
    rm ${plotdir}
fi

if [ ! -L ${plotdir} ]; then
    echo "Creating symlink to plots directory at" $plotsdir
    wwwdir=/eos/user/${USER:0:1}/${USER}/www/WMassAnalysis
    ln -sv $wwwdir $plotdir
fi

loadSingularity=false
usePython2=false
defaultSingularity="/data/shared/singularity/pythonrootdevarch.sif"
python2Singularity=${defaultSingularity/pythonroot/python2root}


for i in "$@" ; do
    if [[ $i == "-s" ]] ; then
        loadSingularity=true
    fi
    if [[ $i == "--p2" ]] ; then
        usePython2=true
    fi
done

if [[ $HOSTNAME == "cmswmass2.cern.ch" ]]; then
    if [[ "$loadSingularity" = "true" ]]; then
	if [[ "$usePython2" = "true" ]]; then
	    echo "Setting up singularity for python 2!"
	    ${python2Singularity}
	else
	    echo "Setting up singularity for python 3!"
	    ${defaultSingularity}
	fi
    else
	echo "Singularity not set! Can call it with any of the following"
	echo "${defaultSingularity} (for python3, recommended)"
	echo "${python2Singularity} (for python2)"
    fi
fi

