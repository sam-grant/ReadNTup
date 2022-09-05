# stn=$1
# nCores=$2

# if [ $# -ne 2 ]; then
#   echo "$# arguments supplied, station and cores, e.g. 'S12' 2."
#   return
# fi

nCores=$1
stn=$2

if [ $# -ne 2 ]; then
  echo "$# arguments supplied, two required"
  return
fi

config="dMu" #Bz" #"Bz"

if [ ! -d ${config} ]; then
    mkdir ${config}
fi

rm -f ${config}/plots_*.root

dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/5.4e-18

# dir=  ${user}/AllDecaysNtuple 

# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Summer2022Campaign/${user}/AllDecaysNtuple 
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Summer2022Campaign/${user}/AllDecaysNtuple 
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Summer2022Campaign/${user}/AllDecaysNtuple 
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Summer2022Campaign/${user}/AllDecaysNtuple 

# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/1.8e-18
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AcceptedDecays/5.4e-18
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/ControlSample/AcceptedDecays

for tree in `ls $dir | sort -V`; do

    id=${tree##*Trees_}
    id=${id%%.root*}

    if [ -f ${config}/plots_${id}_${stn}.root ]; then 
        echo "${config}/plots_${id}_${stn}.root exists, skipping..."
        continue 
    fi

    echo "${dir}/${tree}" "${config}/plots_${id}_${stn}.root" "$stn"
    # break
done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

rm -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_accWeight${stn}_full.root && hadd -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_accWeight${stn}_full.root ${config}/plots_*.root

# rm -f sanityPlots_allDecays_WORLD_250MeV_AQ_noVertCorr.root && hadd -f sanityPlots_allDecays_WORLD_250MeV_AQ_noVertCorr.root ${config}/plots_*.root
# rm -f thetaYvsMomentum_allDecays_WORLD_250MeV_AQ_noVertCorr_accWeight${stn}.root && hadd -f thetaYvsMomentums_allDecays_WORLD_250MeV_AQ_noVertCorr.root ${config}/plots_*.root


# rm -f plots_${config}_${dataset}.root
# hadd -f plots_${config}_${stn}.root ${config}/plots_*_${stn}.root
# rm -f ${config}/plots_*_${dataset}.root

# -------> edmPlots_trackTruth_WORLD_250MeV_AQ_noVertCorr_minus1mm_full.root

# rm -f sanityPlots_allDecays_MainSample_BQ.root && hadd -f sanityPlots_allDecays_MainSample_BQ.root ${config}/plots_*.root
# rm -f sanityPlots_allDecays_Summer2022_${user}_BQ.root && hadd -f sanityPlots_allDecays_Summer2022_${user}_BQ.root ${config}/plots_*.root
# rm -f edmPlots_allDecays_Summer2022_${user}_BQ.root && hadd -f sanityPlots_allDecays_Summer2022_${user}_BQ.root ${config}/plots_*.root
# rm -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_Summer2022_${user}.root && hadd -f edmPlots_allDecays_WORLD_250MeV_AQ_noVertCorr_Summer2022_${user}.root ${config}/plots_*.root
# 

