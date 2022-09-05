nCores=$1 
truth=$2
dataset=$3

if [ $# -ne 3 ]; then
  echo "$# arguments supplied, three required"
  return
fi

config="dMu" #Bz" #"Bz"

if [ ! -d ${config} ]; then
    mkdir ${config}
fi

rm -f ${config}/*.root

# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/$align 

# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/5.4e-18 
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/1.8e-18
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AcceptedDecays/5.4e-18
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/ControlSample/AcceptedDecays
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/5.4e-18
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/1.8e-19
# Plus1mm

dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/5.4e-18

# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Plus1mm
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Minus1mm
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Plus0.1deg
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Minus0.1deg

for tree in `ls $dir | sort -V`; do

    if [[ $tree == attic ]]; then
        continue
    fi

    id=${tree##*truthTrees.}
    id=${id%%.root*}

    echo "${dir}/${tree}" "${config}/plots_${id}.root" "$truth" "$dataset"

     # break

done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

# rm -f thetaYvsMomentum_trackReco_WORLD_250MeV_BQ_noVertCorr.root && hadd -f thetaYvsMomentum_trackReco_WORLD_250MeV_BQ_noVertCorr.root ${config}/plots_*.root
# rm -f thetaYvsMomentum_trackReco_WORLD_250MeV_BQ_noVertCorr.root && hadd -f thetaYvsMomentum_trackReco_WORLD_250MeV_BQ_noVertCorr.root ${config}/plots_*.root
# 
rm -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_full_${dataset}_accWeight.root && hadd -f edmPlots_track${truth}_WORLD_250MeV_BQ_noVertCorr_full_${dataset}_accWeight.root ${config}/plots_*.root
# 
# ./Plots/MC/dMu/5.4e-18/Plots/sanityPlots_truth_full_plus1mm_BQ.root

# rm -f plots_${truth}.root && 

# hadd -f sanityPlots_${truth}_MainSample_plus1mm_BQ.root ${config}/plots_*.root
# hadd -f sanityPlots_${truth}_MainSample_minus1mm_BQ.root ${config}/plots_*.root
# hadd -f sanityPlots_${truth}_MainSample_plus0.1deg_BQ.root ${config}/plots_*.root
# hadd -f sanityPlots_${truth}_MainSample_minus0.1deg_BQ.root ${config}/plots_*.root
# hadd -f sanityPlots_${truth}_MainSample_minus0.1deg_BQ.root ${config}/plots_*.root
# rm -f sanityPlots_${truth}_MainSample_plus1mm_BQ.root && hadd -f sanityPlots_${truth}_MainSample_plus1mm_BQ.root ${config}/plots_*.root
# rm -f sanityPlots_${truth}_MainSample_minus1mm_BQ.root && hadd -f sanityPlots_${truth}_MainSample_minus1mm_BQ.root ${config}/plots_*.root
# 
# rm -f edmPlots_${truth}_MainSample_BQ_noVertCorr_${align}.root && hadd -f edmPlots_${truth}_MainSample_BQ_noVertCorr_${align}.root ${config}/plots_*.root
# rm -f edmPlots_${truth}_Summer2022_${user}_BQ_noVertCorr_${align}.root && hadd -f edmPlots_${truth}_Summer2022_${user}_BQ_noVertCorr_${align}.root ${config}/plots_*.root
# rm -f edmPlots_${truth}_${user}_Summer2022_full_BQ_noVertCorr.root && hadd -f edmPlots_${truth}_${user}_Summer2022_BQ_noVertCorr.root ${config}/plots_*.root

# rm -f edmPlots_${truth}_${user}_Summer2022_full_BQ_noVertCorr_plus1mm.root && hadd -f edmPlots_${truth}_${user}_Summer2022_full_BQ_noVertCorr_plus1mm.root ${config}/plots_*.root
# rm -f edmPlots_${truth}_${user}_Summer2022_full_BQ_noVertCorr_minus1mm.root && hadd -f edmPlots_${truth}_${user}_Summer2022_full_BQ_noVertCorr_minus1mm.root ${config}/plots_*.root
# rm -f edmPlots_${truth}_${user}_Summer2022_full_BQ_noVertCorr_plus0.1deg.root && hadd -f edmPlots_${truth}_${user}_Summer2022_full_BQ_noVertCorr_plus0.1deg.root ${config}/plots_*.root
# rm -f edmPlots_${truth}_${user}_Summer2022_full_BQ_noVertCorr_minus0.1deg.root && hadd -f edmPlots_${truth}_${user}_Summer2022_full_BQ_noVertCorr_minus0.1deg.root ${config}/plots_*.root

# rm -f plots_${truth}.root && hadd -f sanityPlots_${truth}_Plus1mmTestSample_BQ.root ${config}/plots_*.root
# rm -f plots_${truth}.root && hadd -f sanityPlots_${truth}_Minus1mmTestSample_BQ.root ${config}/plots_*.root
# rm -f plots_${truth}.root && hadd -f sanityPlots_${truth}_Plus0.1degTestSample_BQ.root ${config}/plots_*.root
#rm -f plots_${truth}.root && hadd -f sanityPlots_${truth}_Minus0.1degTestSample_BQ.root ${config}/plots_*.root
