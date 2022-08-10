input=$1
output=$2
dataset=$3 

echo "Dataset ${dataset}"
echo "Input ${input}"
echo "Output ${output}"

echo "Command: ./OrthogonalPlotter.exe $input $output $dataset"
# ./Plotter.exe $input $output $dataset
# echo "Command: ./Acceptance.exe $input $output $dataset"
#echo "Command: ./VerticalOffset.exe $input $output $dataset"
#./Plotter.exe $input $output $dataset
./OrthogonalPlotter.exe $input $output $dataset
# ./VerticalOffsetVsTime.exe $input $output $dataset
# ./PlotVerticalPosition.exe $input $output $dataset
##########################################################
# Uncomment if re-submitting and you don't want to overwrite 

#if [[ ! -f $output ]]; then
#    ./Plotter.exe $input $output $dataset
#fi
##########################################################
