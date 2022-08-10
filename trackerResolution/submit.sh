# Get arguments (sam dataset and number of cores to use)
if [ "$#" -ne 1 ]; then
  echo "Number of cores required as argument."
  return
fi

nCores=$1

config="plots" 

rm -rf $config && mkdir $config

dir=/pnfs/GM2/persistent/EDM/MC/TrackerAcceptance/jprice_mc_digit_trackAcceptance

for tree in `ls ${dir} | sort -V`; do

    id=${tree##*trackerAcceptanceTrees.}
    id=${id%%.root*}

    if [ ! -f ${dir}/${tree} ]; then 
        continue
    fi

    echo "${dir}/${tree}" "${config}/trackerResolutionPlots.${id}.root"

done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

rm -f trackerResolutionPlots.root && hadd -f trackerResolutionPlots.root ${config}/trackerResolutionPlots.*.root

