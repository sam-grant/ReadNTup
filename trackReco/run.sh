input=$1
output=$2
dataset=$3 

echo "Dataset ${dataset}"
echo "Input ${input}"
echo "Output ${output}"

# echo "Command: ./Plotter.exe $input $output $dataset"
# ./Plotter.exe $input $output $dataset
echo "Command: ./MomentumDist.exe $input $output $dataset"
./MomentumDist.exe $input $output $dataset
##########################################################
# Uncomment if re-submitting and you don't want to overwrite 
# if [[ ! -f $output ]]; then
#  	./Plotter.exe $input $output $dataset
# fi
##########################################################
