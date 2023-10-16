# Submit histogramming jobs, read from ROOT ntuples on /pnfs
nCores=$1
stn=$2

if [ $# -ne 2 ]; then
  echo "$# arguments supplied, two required"
  return
fi

outDir="output" 

if [ ! -d ${outDir} ]; then
    mkdir ${outDir}
fi

rm -f ${outDir}/plots_*.root

dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/5.4e-18

for tree in `ls $dir | sort -V`; do

    id=${tree##*Trees_}
    id=${id%%.root*}

    if [ -f ${outDir}/plots_${id}_${stn}.root ]; then 
        echo "${outDir}/plots_${id}_${stn}.root exists, skipping..."
        continue 
    fi

    echo "${dir}/${tree}" "${outDir}/plots_${id}_${stn}.root" "$stn"

done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

# --- From Plotter.C ---
# Output format: edmPlots_<refFrame>_<recoType>_<binWidth>_<qualityInfo>_<corrInfo>.root
rm -f edmPlots_allDecays_LAB_250MeV_noQ_randCorr.root && hadd -f edmPlots_allDecays_LAB_250MeV_noQ_randCorr.root ${outDir}/plots_*.root


# rm -f edmPlots_allDecays_LAB_250MeV_noQ_accWeight${stn}.root && hadd -f edmPlots_allDecays_LAB_250MeV_noQ_accWeight${stn}.root ${outDir}/plots_*.root
# rm -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_accWeight${stn}_full.root && hadd -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_accWeight${stn}_full.root ${outDir}/plots_*.root
# rm -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_full.root && hadd -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_full.root ${outDir}/plots_*.root
# rm -f sanityPlots_allDecays_WORLD_250MeV_AQ_noVertCorr.root && hadd -f sanityPlots_allDecays_WORLD_250MeV_AQ_noVertCorr.root ${outDir}/plots_*.root
# rm -f thetaYvsMomentum_allDecays_WORLD_250MeV_AQ_noVertCorr_accWeight${stn}.root && hadd -f thetaYvsMomentums_allDecays_WORLD_250MeV_AQ_noVertCorr.root ${outDir}/plots_*.root
# rm -f thetaYvsMomentum_allDecays_WORLD_250MeV_AQ_noVertCorr.root && hadd -f thetaYvsMomentums_allDecays_WORLD_250MeV_AQ_noVertCorr.root ${outDir}/plots_*.root
# rm -f plots_${outDir}_${dataset}.root
# hadd -f plots_${outDir}_${stn}.root ${outDir}/plots_*_${stn}.root
# rm -f ${outDir}/plots_*_${dataset}.root
# -------> edmPlots_trackTruth_WORLD_250MeV_AQ_noVertCorr_minus1mm_full.root
# rm -f sanityPlots_allDecays_MainSample_BQ.root && hadd -f sanityPlots_allDecays_MainSample_BQ.root ${outDir}/plots_*.root
# rm -f sanityPlots_allDecays_Summer2022_${user}_BQ.root && hadd -f sanityPlots_allDecays_Summer2022_${user}_BQ.root ${outDir}/plots_*.root
# rm -f edmPlots_allDecays_Summer2022_${user}_BQ.root && hadd -f sanityPlots_allDecays_Summer2022_${user}_BQ.root ${outDir}/plots_*.root
# rm -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_Summer2022_${user}.root && hadd -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_Summer2022_${user}.root ${outDir}/plots_*.root