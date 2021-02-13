# Setup recipe

This version of CMGTools strips away all CMSSW dependence, it's just a set of python scripts and macros currently

Still, several paths are expected to be in a path $CMSSW_BASE/src

```bash
gituser=$(git config --global user.github)
git clone git@github.com:${gituser}/cmgtools-lite.git src/CMGTools
cd src/CMGTools
. ./setup.sh
```

In general you will run commands from the director WMass/python/plotter
