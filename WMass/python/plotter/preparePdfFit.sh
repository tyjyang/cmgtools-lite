set -e

pdf=nnpdf31
if [ $# -gt 0 ]; then
    pdf=$1
fi 

firstStep=2
if [ $# -gt 1 ]; then
    firstStep=$2
fi

echo "first step is $firstStep"

if [ $firstStep -le 2 ]; then
    python w-mass-13TeV/testingNano/cfg/makePlotsCFG_systTHn.py -o w-mass-13TeV/testingNano/cfg/plots_fakerate_systTHn_${pdf}.txt --a wmass --pdf-weights $pdf
    python runWmassTHn.py -o plots/testNanoAOD/WmassPlots/W${pdf}_testTHn/ -e postVFP --variables ".*" --plot-file "plots_fakerate_systTHn_${pdf}.txt" --options " --skipPlot " --pdf $pdf
fi

if [ $firstStep -le 3 ]; then
    python w-mass-13TeV/plotFakesTemplate.py plots/testNanoAOD/WmassPlots/W${pdf}_testTHn/fakeRateRegion_systTHn/postVFP/allTHn/plots_fakerate_systTHn_${pdf}.root --skip-plots
fi

if [ $firstStep -le 4 ]; then
    for charge in plus minus; do
        python makeHistogramsWMass.py -i plots/testNanoAOD/WmassPlots/W${pdf}_testTHn/fakeRateRegion_systTHn/postVFP/allTHn/postprocessing/wmass_shapes.root --outdir cards/wmass/W${pdf}_testTHn/ -c $charge --decorrelate-by-charge ".*effStatTnP|.*muR\d+|.*muF\d+"
    done
fi

python w-mass-13TeV/cardMaker.py -i cards/wmass/W${pdf}_testTHn/ -f mu -c "plus,minus" --comb --freezePOIs --mass-nuis massShift100MeV --impacts-mW --skip-fit-data --all-proc-background --dry-run
