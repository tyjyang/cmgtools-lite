#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export CMSSW_BASE=$(dirname $(dirname $dir))
plotsdir=$dir/WMass/python/plotter/plots

if [ ! -e ${plotdir} ]; then
    rm ${plotdir}
fi
if [ ! -L ${plotdir} ]; then
    echo "Creating symlink to plots directory at" $plotsdir
    # Could also be ~/www
    wwwdir=/eos/user/${USER:0:1}/${USER}/www/WMassAnalysis
    ln -s $wwwdir $plotdir
fi

if [[ $HOSTNAME == "cmswmass2.cern.ch" ]]; then
    echo "Setting up singularity"
    #/data/shared/singularity/python2rootdevf32.sif
    /data/shared/singularity/pythonrootdevf32.sif
fi
