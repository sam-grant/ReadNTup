source /cvmfs/gm2.opensciencegrid.org/prod/g-2/setup && setup root v6_12_04e -q e15:prof
source /grid/fermiapp/products/common/etc/setups.sh
setup fife_utils

# Get arguments (sam dataset and number of cores to use)
if [ "$#" -ne 2 ]; then
  echo "SAM dataset & number of cores to use required as arguments"
  return
fi

dataset=$1
nCores=$2

# Make sub-directory (so we can do more than one)
if [ -d $dataset ]; then
  echo "Output directory $dataset already exists."
  return
fi

mkdir $dataset
cd $dataset
pwd

# Get list of sub-runs in this dataset
samweb -e gm2 list-definition-files $dataset > tmp.txt
sort -t_ -k4 tmp.txt > SAMFileList.txt
rm -f tmp.txt

# Check for jobs that haven't been processed already
currJobRun=-1
currRunPath=""
rm -f FileList.txt && touch FileList.txt

for line in `cat SAMFileList.txt`; do

  # Get run and sub-run  
  run=${line##*_}
  run=${run%%.*}
  subrun=${line##*${run}.}
  subrun=${subrun%%.*}

  # If we're in a new run, submit the old lot and clear FileList
  if [ $run -ne $currJobRun ]; then
    if [ -s FileList.txt ]; then
      mv FileList.txt FileList_${currJobRun}.txt
      echo "Stored FileList_${currJobRun}.txt"
      touch FileList.txt
    fi
    fullFilePath=`samweb -e gm2 locate-file $line`
    currRunPath=/pnfs${fullFilePath##*/pnfs}
    currRunPath=${currRunPath%%(*}
    currJobRun=$run
  fi

  # Set path and currJobRun if first time through
  if [ $currJobRun -lt 0 ]; then 
    fullFilePath=`samweb -e gm2 locate-file $line`
    currRunPath=/pnfs${fullFilePath##*/pnfs}
    currRunPath=${currRunPath%%(*}
    currJobRun=$run
  fi

  # Check output exists before adding - ignore those where the tracker had problems

  if [ -f /pnfs/GM2/scratch/users/sgrant/BrAna/CaloBeamPos_NewCuts/${dataset}/${run}/caloBeamPositionPlots_${run}.${subrun}.root ]; then
    echo /pnfs/GM2/scratch/users/sgrant/BrAna/CaloBeamPos_NewCuts/${dataset}/${run}/caloBeamPositionPlots_${run}.${subrun}.root >> FileList.txt
  else 
    echo "File \"/pnfs/GM2/scratch/users/sgrant/BrAna/CaloBeamPos_NewCuts/${dataset}/${run}/caloBeamPositionPlots_${run}.${subrun}.root\" does not exist, but run/sub-run are in dataset $dataset."
    rm -f FileList*.txt
    return;
  fi
  
done

# Catch final file list (as we've left loop)
if [ -s FileList.txt ]; then
 mv FileList.txt FileList_${currJobRun}.txt
 echo "Stored FileList_${currJobRun}.txt"
fi

# Loop over all the file lists we've made and submit them for hadd
for file in `ls FileList_*.txt`; do
  run=${file##*FileList_}
  run=${run%%.*}
  echo $run $dataset
done | xargs -i --max-procs=$nCores bash -c ". ../haddFileList.sh {}"

cd ../