
config="dMu" #Bz" #"Bz"

if [ ! -d ${config} ]; then
    mkdir ${config}
fi

# x4 /pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/1.8e-19
# x4 /pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/5.4e-18 

# X edmPlots_trackReco_WORLD_250MeV_AQ.root
# X edmPlots_trackReco_AAR_250MeV_AQ.root
# X edmPlots_trackTruth_AAR_250MeV_AQ.root
# X edmPlots_trackTruth_WORLD_250MeV_AQ.root
# X edmPlots_trackReco_WORLD_250MeV_BQ.root
# X edmPlots_trackReco_AAR_250MeV_BQ.root
# X edmPlots_trackTruth_AAR_250MeV_BQ.root
# X dmPlots_trackTruth_WORLD_250MeV_BQ.root

dir="/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/5.4e-18"

for tree in `ls $dir | sort -V`; do

    id=${tree##*truthTrees.}
    id=${id%%.root*}
    ./Plotter.exe ${dir}/${tree} ${config}/plots_${id}.root

done

rm -f plots_${config}.root
hadd -f plots_${config}.root ${config}/plots_*.root
rm -f ${config}/plots_*.root
