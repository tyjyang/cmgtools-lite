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
defaultSingularity="/scratch/singularity/pythonrootarchdev"

for i in "$@" ; do
    if [[ $i == "-s" ]] ; then
        loadSingularity=true
    fi
done

if [[ "$loadSingularity" = "true" ]]; then
    echo "Setting up singularity for python 3!"
    singularity run ${defaultSingularity}
fi
