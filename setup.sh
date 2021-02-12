#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export CMSSW_BASE=$(dirname $(dirname $dir))
plotsdir=$CMSSW_BASE/WMass/python/plotter/plots

if [ ! -e ${plotdir} ]; then
    rm ${plotdir}
fi
if [ ! -L ${plotdir} ]; then
    # Could also be ~/www
    wwwdir=/eos/user/${USER:0:1}/${USER}/www/WMassAnalysis
    ln -s wwwdir $plotdir
fi

/data/shared/singularity/python2rootdevf32.sif
