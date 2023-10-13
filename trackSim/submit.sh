# Submit trackTruth/Reco plotting jobs 

nCores=$1 
truthStr=$2
# alignDir=$3
# alignStr=$4

if [ $# -ne 2 ]; then
  echo "$# arguments supplied, two required (N cores and Truth/Reco)"
  return
fi

outDir="output" 

if [ ! -d ${outDir} ]; then
    mkdir ${outDir}
fi

rm -f ${config}/*.root

dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/5.4e-18

# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/$alignDir
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Plus1mm
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Minus1mm
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Plus0.1deg
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Minus0.1deg

for tree in `ls $dir | sort -V`; do

    id=${tree##*truthTrees.}
    id=${id%%.root*}

    echo "${dir}/${tree}" "${outDir}/plots_${id}.root" "$truthStr"

done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

hadd -f edmPlots_track${truthStr}_LAB_250MeV_BQ_randCorr.root ${outDir}/plots_*.root

# rm -f thetaYvsMomentum_track${truth}_WORLD_250MeV_BQ_noVertCorr_full.root && hadd -f thetaYvsMomentum_track${truth}_WORLD_250MeV_BQ_noVertCorr_full.root ${config}/plots_*.root
# rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_full.root ${config}/plots_*.root
# rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_${alignStr}_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_${alignStr}_full.root ${config}/plots_*.root
#rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_plus1mm_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_plus1mm_full.root ${config}/plots_*.root
# rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_minus1mm_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_minus1mm_full.root ${config}/plots_*.root
# rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_plus0.1deg_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_plus0.1deg_full.root ${config}/plots_*.root
# rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_minus0.1deg_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_minus0.1deg_full.root ${config}/plots_*.root
