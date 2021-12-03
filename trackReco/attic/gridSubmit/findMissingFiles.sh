# Check that samweb and root are installed
command -v samweb >/dev/null 2>&1 || { echo >&2 "samweb not found. You want:"; echo "source /grid/fermiapp/products/common/etc/setups.sh"; echo "setup fife_utils"; return; }

# Get sam dataset number
if [ "$#" -ne 1 ]; then
  echo "SAM dataset (only) required as argument"
  return
fi
dataset=$1
pnfsOutDir=/pnfs/GM2/scratch/users/${USER}/BrAna/CaloBeamPos_NewCuts/${dataset}
#echo "Checking output files for $dataset in $pnfsOutDir"

# Get list of sub-runs in this dataset
samweb -e gm2 list-definition-files $dataset > tmp.txt
sort -t_ -k4 tmp.txt > SAMFileList.txt
rm -f tmp.txt

# Check for jobs that haven't been processed already
for line in `cat SAMFileList.txt`; do

  # Get run and sub-run  
  run=${line##*_}
  run=${run%%.*}
  subrun=${line##*${run}.}
  subrun=${subrun%%.*}

  if [ ! -f ${pnfsOutDir}/${run}/caloBeamPositionPlots_${run}.${subrun}.root ]; then
      echo ${pnfsOutDir}/${run}/caloBeamPositionPlots_${run}.${subrun}.root
  fi


done

