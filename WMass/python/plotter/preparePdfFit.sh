pdf=$1
set -e

runFull=0
if [ $# -gt 1 ]; then
    runFull=$2
fi

if [ $runFull -gt 0 ]; then
    python w-mass-13TeV/testingNano/cfg/makePlotsCFG_systTHn.py -o w-mass-13TeV/testingNano/cfg/plots_fakerate_systTHn_${pdf}.txt --a wmass --pdf-weights $pdf
    python runWmassTHn.py -o plots/testNanoAOD/WmassPlots/W${pdf}_testTHn/ -e postVFP --variables ".*" --plot-file "plots_fakerate_systTHn_${pdf}.txt" --options " --skipPlot --no-rdf-runGraphs "
fi

python w-mass-13TeV/plotFakesTemplate.py plots/testNanoAOD/WmassPlots/W${pdf}_testTHn/fakeRateRegion_systTHn/postVFP/allTHn/plots_fakerate_systTHn_${pdf}.root
exit 1

for charge in plus minus; do
    python makeHistogramsWMass.py -i plots/testNanoAOD/WmassPlots/W${pdf}_testTHn/fakeRateRegion_systTHn/postVFP/allTHn/postprocessing/wmass_shapes.root --outdir cards/wmass/W${pdf}_testTHn/ -c plus --decorrelate-by-charge ".*effStatTnP|.*muR\d+|.*muF\d+"
done
python w-mass-13TeV/cardMaker.py -i cards/wmass/W${pdf}_noPDFonZ_noDzCut/  -f mu -c "plus,minus" --comb --freezePOIs --mass-nuis massShift100MeV --impacts-mW --skip-fit-data --all-proc-background 

