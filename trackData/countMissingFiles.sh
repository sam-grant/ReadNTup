dataset=$1

if [[ $dataset == "Run-1a" ]]; then
  datasetName="gm2pro_daq_full_run1_60h_5039A_GLdocDB16021-v2"
elif [[ $dataset == "Run-1b" ]]; then
  datasetName="gm2pro_daq_full_run1_HighKick_5042B_GLdocDB20949-v3"
elif [[ $dataset == "Run-1c" ]]; then
  datasetName="gm2pro_daq_full_run1_9d_5040A_GLdocDB17018-v3"
elif [[ $dataset == "Run-1d" ]]; then
  datasetName="gm2pro_daq_full_run1_EndGame_5042B_GLdocDB20839-v1"
fi

a=`cat txt/${datasetName}.txt | wc -l`
b=`ls $datasetName | wc -l`

echo "${b} out of ${a} files written"
