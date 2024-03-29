# Hadd files in FileList_12345.txt (made by haddSubRunsSAMList.sh)

# Get dataset from argument
dataset=$1

datasetName=""
if [[ $dataset == "Run-1a" ]]; then
  datasetName="gm2pro_daq_full_run1_60h_5039A_GLdocDB16021-v2"
elif [[ $dataset == "Run-1b" ]]; then
  datasetName="gm2pro_daq_full_run1_HighKick_5042B_GLdocDB20949-v3"
elif [[ $dataset == "Run-1c" ]]; then
  datasetName="gm2pro_daq_full_run1_9d_5040A_GLdocDB17018-v3"
elif [[ $dataset == "Run-1d" ]]; then
  datasetName="gm2pro_daq_full_run1_EndGame_5042B_GLdocDB20839-v1"
fi

files=""

for file in `ls $datasetName`; do

  files=$files" "${datasetName}/$file

done

# vaecho $files
# Hadd these files (using Joe's script that does a few at 
hadd.sh edmPlots_${dataset}.root $files
