# Get arguments (sam dataset and number of cores to use)
if [ "$#" -ne 1 ]; then
  echo "Number of cores required as argument."
  return
fi

nCores=$1

config="plots" 

rm -rf $config && mkdir $config

dir=/pnfs/GM2/persistent/EDM/MC/TrackerAcceptance/jprice_mc_digit_trackAcceptance
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/5.4e-18
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/0mm
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/0mm

for tree in `ls ${dir} | sort -V`; do

    id=${tree##*trackerAcceptanceTrees.}
    id=${id%%.root*}

    if [ ! -f ${dir}/${tree} ]; then 
        continue
    fi

    echo "${dir}/${tree}" "${config}/trackerAcceptancePlots.${id}.root"
    # ./Plotter.exe ${dir}/${tree} ${config}/trackerAcceptancePlots.${id}.root
    #./PlotterMisalignment.exe ${dir}/${tree} ${config}/trackerAcceptancePlots.${id}.root
    # break
done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

rm -f trackerAcceptancePlots.root && hadd -f trackerAcceptancePlots.root ${config}/trackerAcceptancePlots.*.root

# dir="/gm2/data/g2be/Production/Trees/Run1" 

# for run in `cat txt/${datasetName}.txt`; do
  
#   echo "${dir}/trackRecoTrees_${run}.root" "${datasetName}/trackRecoPlots_${run}.root" $dataset
  
# done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"
