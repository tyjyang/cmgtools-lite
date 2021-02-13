#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export CMSSW_BASE=$(dirname $(dirname $dir))
plotdir=$dir/WMass/python/plotter/plots

if [ -L $plotdir ] && [ ! -e ${plotdir} ]; then
    rm ${plotdir}
fi

if [ ! -L ${plotdir} ]; then
    echo "Creating symlink to plots directory at" $plotsdir
    # Could also be ~/www
    wwwdir=/eos/user/${USER:0:1}/${USER}/www/WMassAnalysis
    ln -s $wwwdir $plotdir
fi

p3=1
if [[ $HOSTNAME == "cmswmass2.cern.ch" ]]; then
    if [ $p3 -eq 1 ]; then
        echo "Setting up singularity for python3"
        /data/shared/singularity/pythonrootdevf32.sif
    else
        echo "Setting up singularity for python2"
        /data/shared/singularity/python2rootdevf32.sif
    fi
fi
