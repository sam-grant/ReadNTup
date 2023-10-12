# Submit histogramming jobs, read from ROOT ntuples on /pnfs
nCores=$1
stn=$2

if [ $# -ne 2 ]; then
  echo "$# arguments supplied, two required"
  return
fi

config="output" 

if [ ! -d ${config} ]; then
    mkdir ${config}
fi

rm -f ${config}/plots_*.root

dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/5.4e-18

for tree in `ls $dir | sort -V`; do

    id=${tree##*Trees_}
    id=${id%%.root*}

    if [ -f ${config}/plots_${id}_${stn}.root ]; then 
        echo "${config}/plots_${id}_${stn}.root exists, skipping..."
        continue 
    fi

    echo "${dir}/${tree}" "${config}/plots_${id}_${stn}.root" "$stn"

done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

# --- From Plotter.C ---
# Output format: edmPlots_<refFrame>_<recoType>_<binWidth>_<qualityInfo>_<corrInfo>.root
rm -f edmPlots_allDecays_LAB_250MeV_noQ.root && hadd -f edmPlots_allDecays_LAB_250MeV_noQ.root ${config}/plots_*.root
# rm -f edmPlots_allDecays_LAB_250MeV_noQ_accWeight${stn}.root && hadd -f edmPlots_allDecays_LAB_250MeV_noQ_accWeight${stn}.root ${config}/plots_*.root

# rm -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_accWeight${stn}_full.root && hadd -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_accWeight${stn}_full.root ${config}/plots_*.root

# rm -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_full.root && hadd -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_full.root ${config}/plots_*.root

# rm -f sanityPlots_allDecays_WORLD_250MeV_AQ_noVertCorr.root && hadd -f sanityPlots_allDecays_WORLD_250MeV_AQ_noVertCorr.root ${config}/plots_*.root
# rm -f thetaYvsMomentum_allDecays_WORLD_250MeV_AQ_noVertCorr_accWeight${stn}.root && hadd -f thetaYvsMomentums_allDecays_WORLD_250MeV_AQ_noVertCorr.root ${config}/plots_*.root

# rm -f thetaYvsMomentum_allDecays_WORLD_250MeV_AQ_noVertCorr.root && hadd -f thetaYvsMomentums_allDecays_WORLD_250MeV_AQ_noVertCorr.root ${config}/plots_*.root


# rm -f plots_${config}_${dataset}.root
# hadd -f plots_${config}_${stn}.root ${config}/plots_*_${stn}.root
# rm -f ${config}/plots_*_${dataset}.root

# -------> edmPlots_trackTruth_WORLD_250MeV_AQ_noVertCorr_minus1mm_full.root

# rm -f sanityPlots_allDecays_MainSample_BQ.root && hadd -f sanityPlots_allDecays_MainSample_BQ.root ${config}/plots_*.root
# rm -f sanityPlots_allDecays_Summer2022_${user}_BQ.root && hadd -f sanityPlots_allDecays_Summer2022_${user}_BQ.root ${config}/plots_*.root
# rm -f edmPlots_allDecays_Summer2022_${user}_BQ.root && hadd -f sanityPlots_allDecays_Summer2022_${user}_BQ.root ${config}/plots_*.root
# rm -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_Summer2022_${user}.root && hadd -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_Summer2022_${user}.root ${config}/plots_*.root
# 

