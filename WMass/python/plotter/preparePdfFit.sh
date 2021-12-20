set -e

pdf=nnpdf31
if [ $# -gt 0 ]; then
    pdf=$1
fi 

firstStep=2
if [ $# -gt 1 ]; then
    firstStep=$2
fi

fitLabel="default"
if [ $# -gt 2 ]; then
    fitLabel=$3
fi

extraArgs=""
if [ $# -gt 3 ]; then
    extraArgs="--cf $4"
fi

#plots/testNanoAOD/WmassPlots/W${pdf}/postVFP/${label}/wmass_shapes.root 
fakeFile=plots/testNanoAOD/WmassPlots/W${pdf}_testTHn/fakeRateRegion_systTHn/postVFP/allTHn/plots_fakerate_systTHn_${pdf}.root 
histFile=plots/testNanoAOD/WmassPlots/W${pdf}_testTHn/fakeRateRegion_systTHn/postVFP/allTHn/postprocessing/wmass_shapes.root
cardInputs=cards/wmass/W${pdf}_${fitLabel} 
echo "first step is $firstStep"

if [ $firstStep -le 2 ]; then
    python w-mass-13TeV/testingNano/cfg/makePlotsCFG_systTHn.py -o w-mass-13TeV/testingNano/cfg/plots_fakerate_systTHn_${pdf}.txt --a wmass --pdf-weights $pdf --muonScaleUnc dummy1Bin
    python runWmassTHn.py -o plots/testNanoAOD/WmassPlots/W${pdf}/ -e postVFP --variables ".*" --plot-file "plots_fakerate_systTHn_${pdf}.txt" --options " --skipPlot " --pdf $pdf
fi

if [ $firstStep -le 3 ]; then
    python w-mass-13TeV/plotFakesTemplate.py $fakeFile --skip-plots
fi

if [ $firstStep -le 4 ]; then
    for charge in plus minus; do
        python makeHistogramsWMass.py -i $histFile --outdir $cardInputs -c $charge --decorrelate-by-charge ".*effStatTnP|.*muR\d+|.*muF\d+"
    done
fi

python w-mass-13TeV/cardMaker.py -i $cardInputs -f mu -c "plus,minus" --comb --mass-nuis massShift100MeV --all-proc-background --dry-run $extraArgs
