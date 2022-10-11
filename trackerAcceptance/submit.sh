# Get arguments (sam dataset and number of cores to use)
if [ "$#" -ne 2 ]; then
  echo "Number of cores and truth required as argument."
  return
fi

nCores=$1
truth=$2

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

    echo "${dir}/${tree}" "${config}/trackerAcceptancePlots.${truth}.${id}.root" "$truth" 

done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

rm -f trackerAcceptancePlots.${truth}.fine.root && hadd -f trackerAcceptancePlots.${truth}.fine.root ${config}/trackerAcceptancePlots.${truth}.*.root