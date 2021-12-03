dataset=$1
fileListID=$2

files=""
for file in `cat splitFileList_${dataset}_${fileListID}`; do
  	files=$files" "$file
done

hadd.sh plots_${dataset}.${fileListID}.root $files
