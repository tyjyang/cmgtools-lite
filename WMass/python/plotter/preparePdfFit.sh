pdf=$1
set -e

python w-mass-13TeV/testingNano/cfg/makePlotsCFG_systTH3.py -o w-mass-13TeV/testingNano/cfg/plots_fakerate_systTH3_${pdf}.txt --a wmass -b 48,-2.4,2.4,116,26,142 --ptVar "Muon_pt[goodMuons][0]+29.0*regionIsoMt(Muon_pfRelIso04_all[goodMuons][0]<0.15,transverseMass<40)" --pdf-weights $pdf --ptVarScaleTest "customPtTest+29.0*regionIsoMt(Muon_pfRelIso04_all[goodMuons][0]<0.15,transverseMass<40)"
python w-mass-13TeV/testingNano/cfg/makePlotsCFG_systTH3.py -o w-mass-13TeV/testingNano/cfg/plots_fakerate_systTH3_nnpdf31.txt --a wmass -b 48,-2.4,2.4,116,26,142 --ptVar "Muon_pt[goodMuons][0]+29.0*regionIsoMt(Muon_pfRelIso04_all[goodMuons][0]<0.15,transverseMass<40)" --pdf-weights nnpdf31 --ptVarScaleTest "customPtTest+29.0*regionIsoMt(Muon_pfRelIso04_all[goodMuons][0]<0.15,transverseMass<40)"

for charge in plus minus; do
    python runFakeRate.py -o plots/testNanoAOD/WmassPlots/W${pdf}_noPDFonZ_noDzCut/ -e postVFP --variables ".*" --plot-file "plots_fakerate_systTH3_${pdf}.txt" --options " --skipPlot " -c ${charge} -s
    python w-mass-13TeV/plotFakesTemplate.py plots/testNanoAOD/WmassPlots/W${pdf}_noPDFonZ_noDzCut/fakeRateRegion_systTH3/postVFP/${charge}/plots_fakerate_systTH3_${pdf}.root -b "29,26,55"
    python makeHistogramsWMass.py -i plots/testNanoAOD/WmassPlots/W${pdf}_noPDFonZ_noDzCut/fakeRateRegion_systTH3/postVFP/${charge}/postprocessing/wmass_shapes.root -o cards/wmass/W${pdf}_noPDFonZ_noDzCut/Wmunu_${charge}_shapes.root -c ${charge} --decorrelate-by-charge ".*effStatTnP|.*muR\d+|.*muF\d+" 
done
python w-mass-13TeV/cardMaker.py -i cards/wmass/W${pdf}_noPDFonZ_noDzCut/  -f mu -c "plus,minus" --comb --freezePOIs --mass-nuis massShift100MeV --impacts-mW --skip-fit-data --all-proc-background 

