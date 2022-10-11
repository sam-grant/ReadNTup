# Get arguments (sam dataset and number of cores to use)
if [ "$#" -ne 2 ]; then
  echo "Dataset and number of and cores required as arguments. EG: . submit.sh Run-1a 2 "
  return
fi

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

if [ ! -d ${datasetName} ]; then
  mkdir ${datasetName}

##########################################################
# Comment if re-submitting and you don't want to overwrite 
else 
rm -f ${datasetName}/*.root
##########################################################

fi

# echo $datasetName

dir="/gm2/data/g2be/Production/Trees/Run1" 

echo "Run, Subruns, All fills, Laser fills, Tracks, Quality tracks, Vertices, Quality vertices"

for run in `cat txt/${datasetName}.txt`; do
  
  echo "${dir}/trackRecoTrees_${run}.root" "${datasetName}/trackRecoPlots_${run}.root" $dataset
  
done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

sleep 2

rm -f edmPlots_${dataset}_250MeV_1000_2500MeV_randomised_BQ.root && hadd -f rm -f edmPlots_${dataset}_250MeV_1000_2500MeV_randomised_BQ.root ${datasetName}/trackRecoPlots_*.root

# hadd -f thetaYvsMomentum_${dataset}_BQ_noVertCorr.root ${datasetName}/trackRecoPlots_*.root

# hadd -f plots_${dataset}.root ${datasetName}/trackRecoPlots_*.root

# hadd -f count_${dataset}.root ${datasetName}/trackRecoPlots_*.root