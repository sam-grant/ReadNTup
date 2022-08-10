files=""
for i in {0..6}; do
	files=$files" "dMu/plots_trackTruthTrees.${i}.root
	hadd -f count_trackTruth_WORLD_250MeV_BQ_noVertCorr.${i}.root $files
done