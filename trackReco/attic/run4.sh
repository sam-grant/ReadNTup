dataset="gm2pro_daq_full_run1_EndGame_5042B_GLdocDB20839-v1"

if [ ! -d ${dataset} ]; then
	mkdir ${dataset}
else 
	rm -f ${dataset}/*.root
fi

dir="/gm2/data/g2be/Production/Trees/Run1" 

for run in `cat txt/${dataset}.txt`; do

	echo "Dataset ${dataset}"

	ls ${dir}/trackRecoTrees_${run}.root
	
	./VerticalOffset.exe "${dir}/trackRecoTrees_${run}.root" "${dataset}/trackRecoPlots_${run}.root" "Run-1d"
	
done