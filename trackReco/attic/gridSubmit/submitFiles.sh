exeFile="Plotter.exe"

# One argument allowed (the dataset) 
if [ $# -ne 1 ]; then
  echo "$# arguments supplied, one allowed."
  return
fi

# Get dataset name
dataset=$1

if [[ $dataset == "Test" ]]; then
  fileList="txt/testFileList.txt"
elif [[ $dataset == "Run-1a" ]]; then
  fileList="txt/gm2pro_daq_full_run1_60h_5039A_GLdocDB16021-v2.txt"
elif [[ $dataset == "Run-1b" ]]; then
  datasetName="txt/gm2pro_daq_full_run1_HighKick_5042B_GLdocDB20949-v3.txt"
elif [[ $dataset == "Run-1c" ]]; then
  datasetName="txt/gm2pro_daq_full_run1_9d_5040A_GLdocDB17018-v3.txt"
elif [[ $dataset == "Run-1d" ]]; then
  datasetName="txt/gm2pro_daq_full_run1_EndGame_5042B_GLdocDB20839-v1.txt"
fi

#Make sure run list exists
if [ ! -f $fileList ] || [ `cat $fileList | wc -l` -eq 0 ]; then
  echo "Input file list \"${fileList}\" is required.  This was not found or has 0 entries"
  return
fi

rm -f fileList.txt && cp $fileList fileList.txt

# Make output directory
pnfsOutDir=/pnfs/GM2/scratch/users/${USER}/Run1EDMPlots/${dataset}

if [ ! -d $pnfsOutDir ]; then
  mkdir -p $pnfsOutDir
  chmod -R g+w $pnfsOutDir
  echo "Output files will appear in $pnfsOutDir"
fi

rm -f newFileList.txt && touch newFileList.txt 
# Make sure that there's not already an output file that we'd overwrite, if there is, remove it 
while read run; do 

  if [ -f ${pnfsOutDir}/${dataset}/trackRecoPlots_${run}.root ]; then 
    echo "${pnfsOutDir}/${dataset}/trackRecoPlots_${run}.root already exists and would be overwritten. Removing file from list. "
  else
    echo $file >> newFileList.txt
  fi

done < fileList.txt

# A bit clumsy here but it just checks to see if the lists are the same length
a=`cat fileList.txt | wc -l`
b=`cat newFileList.txt | wc -l`

if [ "$a" -eq "$b" ]; then
  echo "No possible overwrites found for this sub-run."
  rm -r newFileList.txt
else
  echo "Overwrites found for this sub-run! Using modified list."
  mv newFileList.txt fileList.txt
fi

# Copy .exe file to pnfs so we can get it from grid jobs
if [ -f ${pnfsOutDir}/${exeFile} ]; then
  rm -f ${pnfsOutDir}/${exeFile} 
fi

sleep 5

ifdh cp ${exeFile} ${pnfsOutDir}/${exeFile}

echo "Copied ${exeFile} to $pnfsOutDir"

# Split into lists with 50 each (and submit one of these each time)
rm -f SplitFileList*
split FileList.txt -l 50 -a 3 -d SplitFileList

# what does this do?
# for file in `ls SplitFileList*`; do 
#   mv $file ${file}.txt
# done

for list in `ls SplitFileList*`; do 

  rm -f submission${list}.txt && touch submission${list}.txt

  while read run; do
    echo /gm2/data/g2be/Production/Trees/Run1/trackRecoTrees_${run}.root >> submission${list}.txt
    # echo root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/${longFileName} >> xrootdFileList${fileNum}.txt
  done < $list

  # Copy filelist, fcl & .so file to pnfs so we can get it from grid jobs
  if [ -f ${pnfsOutDir}/submission${list}.txt ]; then
    rm -f ${pnfsOutDir}/submission${list}.txt
  fi

  ifdh cp submission${list}.txt ${pnfsOutDir}/submission${list}.txt
  echo "Copied submission${list}.txt to $pnfsOutDir"

# Make script that we'll want to run on the grid
cat <<EOF > run${fileList}.sh

# Setup stuff to copy files back and forth
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup
setup ifdhc

# Copy submission${list}.txt to grid node 
if [ ! -z "\`ifdh ls ${pnfsOutDir}/submission${list}.txt\`" ]; then
  ifdh cp ${pnfsOutDir}/submission${list}.txt submission${list}.txt
else 
  echo "${pnfsOutDir}/submission${list}.txtdoes not exist"
fi

# Copy exe file
if [ ! -z "\`ifdh ls ${pnfsOutDir}/${exeFile}\`" ]; then
  ifdh cp ${pnfsOutDir}/${exeFile} ./${exeFile}
else 
  echo "${pnfsOutDir}/${exeFile}does not exist"
fi

# Setup gm2
source /cvmfs/gm2.opensciencegrid.org/prod/g-2/setup
setup gm2 v9_71_00 -q prof

# Add current directory to LD_LIBRARY_PATH for so file we copied
export LD_LIBRARY_PATH=\`pwd\`:\$LD_LIBRARY_PATH

# Some nodes have newer kernels
case `uname -r` in
  3.*) export UPS_OVERRIDE="-H Linux64bit+2.6-2.12";;
  4.*) export UPS_OVERRIDE="-H Linux64bit+2.6-2.12";;
esac

# Run jobs
while read run; do 

  echo "Processing \$run"

  ./${exeFile} \${run} 

  # Only copy back if job completed successfully (output from last command was 0)
  if [ \$? -eq 0 ]; then
    if [ -f trackRecoPlots_\${run}.root ]; then
      ifdh mv trackRecoPlots_\${run}.root ${pnfsOutDir}/trackRecoPlots_\${run}.root
    else 
      echo trackRecoPlots_\${run}.root does not exist.  ls reads:"
      ls
    fi
  fi
done < submission${list}.txt # xrootdFileList${fileNum}.txt

EOF

  if [ -f ${pnfsOutDir}/run${fileList}.sh ]; then
    rm -f ${pnfsOutDir}/run${fileList}.sh
  fi

  ifdh cp run${fileList}.sh ${pnfsOutDir}/run${fileList}.sh

  while [ ! -f ${pnfsOutDir}/run${fileList}.sh ]; do
    echo "${pnfsOutDir}/run${fileList}.sh not found after ifdh cp.  Sleeping 5..."
    sleep 5
  done

  echo "Running submission"

  jobsub_submit -N 1 -G gm2 --OS=SL7 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --expected-lifetime=30m --role=Analysis file://${pnfsOutDir}/runJob${fileNum}.sh
   
  rm -f runJob${fileNum}.sh
  rm -f xrootdFileList${fileNum}.txt
  rm -f SplitFileList${fileNum}.txt

done
