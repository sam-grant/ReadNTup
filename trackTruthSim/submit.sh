nCores=$1 
truth=$2
# alignDir=$3
# alignStr=$4


if [ $# -ne 2 ]; then
  echo "$# arguments supplied, two required (N cores and Truth/Reco)"
  return
fi

config="dMu" #Bz" #"Bz"

if [ ! -d ${config} ]; then
    mkdir ${config}
fi

rm -f ${config}/*.root

dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/5.4e-18

# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/$alignDir
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Plus1mm
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Minus1mm
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Plus0.1deg
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Minus0.1deg

echo "Track vertices, quality track vertices"

for tree in `ls $dir | sort -V`; do

    if [[ $tree == attic ]]; then
        continue
    fi

    # Skip random file from dominika, not sure what this is
    if [[ $tree == "trackSimTreesfull.root" ]]; then
#	echo $tree
	continue
    fi

    id=${tree##*truthTrees.}
    id=${id%%.root*}

    echo "${dir}/${tree}" "${config}/plots_${id}.root" "$truth"

     # break

done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

rm -f thetaYvsMomentum_track${truth}_WORLD_250MeV_BQ_noVertCorr_full.root && hadd -f thetaYvsMomentum_track${truth}_WORLD_250MeV_BQ_noVertCorr_full.root ${config}/plots_*.root

# rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_full.root ${config}/plots_*.root

# rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_${alignStr}_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_${alignStr}_full.root ${config}/plots_*.root

#rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_plus1mm_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_plus1mm_full.root ${config}/plots_*.root

# rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_minus1mm_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_minus1mm_full.root ${config}/plots_*.root
# rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_plus0.1deg_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_plus0.1deg_full.root ${config}/plots_*.root
# rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_minus0.1deg_full.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_minus0.1deg_full.root ${config}/plots_*.root
