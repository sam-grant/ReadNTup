# Get arguments (sam dataset and number of cores to use)
if [ "$#" -ne 1 ]; then
  echo "Number of cores to use required as an argument"
  return
fi

nCores=$1

########

dataset="gm2pro_daq_full_run1_60h_5039A_GLdocDB16021-v2" #Bz" #"Bz"

if [ ! -d ${dataset} ]; then
    mkdir ${dataset}
fi

dir="/gm2/data/g2be/Production/Trees/Run1" 

########

for run in `cat txt/${dataset}.txt`; do
  
  # ./Plotter.exe 
  echo ${dir}/trackRecoTrees_${run}.root ${dataset}/trackRecoPlots_${run}.root
  
done | xargs -i --max-procs=$nCores bash -c ". RunJob.sh {}"

########
