#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export CMSSW_BASE=$(dirname $(dirname $dir))
plotdir=$CMSSW_BASE/src/CMGTools/WMass/python/plotter/plots

loadSingularity=false
usePython2=false

defaultSingularity="/data/shared/singularity/pythonrootdevf32.sif" # for python 3
python2Singularity=${defaultSingularity/pythonroot/python2root}


for i in "$@" ; do
    if [[ $i == "-s" ]] ; then
        loadSingularity=true
    fi
    if [[ $i == "--p2" ]] ; then
        usePython2=true
    fi
done


if [ -e ${plotdir} ]; then
    # this part should not be entered, unless the old link was committed
    echo "Removing spurious link $plotdir"
    echo "It will be recreated"
    rm ${plotdir}
fi
if [ ! -L ${plotdir} ]; then
    # Could also be ~/www
    echo "Creating link to www folder for plot storage"
    wwwdir=/eos/user/${USER:0:1}/${USER}/www/WMassAnalysis
    ln -sv $wwwdir $plotdir
fi

if [[ "$loadSingularity" = "true" ]]; then
    if [[ "$usePython2" = "true" ]]; then
	echo "Setting singularity for python 2!"
	${python2Singularity}
    else
	echo "Setting singularity for python 3!"
	${defaultSingularity}
    fi
else
    echo "Singularity not set! Can call it with any of the following"
    echo "${defaultSingularity} (for python3, recommended)"
    echo "${python2Singularity} (for python2)"
fi
