# Get arguments (sam dataset and number of cores to use)
if [ "$#" -ne 2 ]; then
  echo "Number of cores and truth required as argument."
  return
fi

nCores=$1
truth=$2 # Default is to use truth for acceptance studies

output="output" 

if [ ! -d $output ]; then
  mkdir $output
fi

rm -f $output/*

# Use Joe's acceptance sample, since this has the same number of events between decays and tracks
dir=/pnfs/GM2/persistent/EDM/MC/TrackerAcceptance/jprice_mc_digit_trackAcceptance

for tree in `ls ${dir} | sort -V`; do

    id=${tree##*trackerAcceptanceTrees.}
    id=${id%%.root*}

    if [ ! -f ${dir}/${tree} ]; then 
        continue
    fi

    echo "${dir}/${tree}" "${output}/trackerAcceptancePlots.${truth}.${id}.root" "$truth" 

done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

rm -f trackerAcceptancePlots.${truth}.root && hadd -f trackerAcceptancePlots.${truth}.root ${output}/trackerAcceptancePlots.${truth}.*.root