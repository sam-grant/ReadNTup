dataset=$1
run=$2

# Get arguments (sam dataset and number of cores to use)
if [ "$#" -ne 2 ]; then
  echo "SAM dataset & run required (eg Run1) to use required as arguments"
  return
fi

echo "Copying ROOT files in /pnfs/GM2/scratch/users/sgrant/BrAna/CaloBeamPos_NewCuts/${dataset}/ to /pnfs/GM2/persistent/EDM/Data/RadialField/CaloBeamPos_NewCuts/${run}"

cp /pnfs/GM2/scratch/users/sgrant/BrAna/CaloBeamPos_NewCuts/${dataset}/*.root /pnfs/GM2/persistent/EDM/Data/RadialField/CaloBeamPos_NewCuts/${run}
