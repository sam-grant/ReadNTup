# DO NOT USE IT FAILS TO INCLUDE THE PART.1. INTERMEDIATE FILES FOR SOME REASON

dataset=$1
nCores=$2

if [[ $dataset == "Run-1a" ]]; then
  datasetName="gm2pro_daq_full_run1_60h_5039A_GLdocDB16021-v2"
elif [[ $dataset == "Run-1b" ]]; then
  datasetName="gm2pro_daq_full_run1_HighKick_5042B_GLdocDB20949-v3"
elif [[ $dataset == "Run-1c" ]]; then
  datasetName="gm2pro_daq_full_run1_9d_5040A_GLdocDB17018-v3"
elif [[ $dataset == "Run-1d" ]]; then
  datasetName="gm2pro_daq_full_run1_EndGame_5042B_GLdocDB20839-v1"
fi

# make main file list
rm -f fileList_${dataset} && touch fileList_${dataset}
for file in `ls ../$datasetName`; do
	echo ../${datasetName}/${file} >> fileList_${dataset}
done

# split files in blocks of 50
rm -f splitFileList_${dataset}_*
split fileList_${dataset} -l 50 -a 3 -d splitFileList_${dataset}_

for list in `ls splitFileList_${dataset}_*`; do
	# echo dataset
	echo $dataset ${list##*_}
done | xargs -i --max-procs=$nCores bash -c ". haddFileList.sh {}"

# finish up
files=""
for file in `ls plots_${dataset}.*.root`; do
  	files=$files" "$file
done

hadd -f plots_${dataset}.root $files

# hadd trackTruthTrees.3.root splitFileList0*/*.root
