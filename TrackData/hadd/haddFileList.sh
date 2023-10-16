# DO NOT USE IT FAILS TO INCLUDE THE PART.1. INTERMEDIATE FILES FOR SOME REASON

dataset=$1
fileListID=$2

files=""
for file in `cat splitFileList_${dataset}_${fileListID}`; do
  	files=$files" "$file
done

hadd.sh plots_${dataset}.${fileListID}.root $files
