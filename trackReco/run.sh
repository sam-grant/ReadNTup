input=$1
output=$2
truth=$3
dataset=$4 


echo "Dataset ${dataset}"
echo "Input ${input}"
echo "Output ${output}"

# echo "Command: ./ThetaYvsMomentum.exe $input $output $dataset"
# ./ThetaYvsMomentum.exe $input $output $dataset
echo "Command: ./Plotter.exe $input $output $truth $dataset"
./Plotter.exe $input $output $truth $dataset
# echo "Command: ./Acceptance.exe $input $output $dataset"
#echo "Command: ./VerticalOffset.exe $input $output $dataset"
#./Plotter.exe $input $output $dataset
# ./OrthogonalPlotter.exe $input $output $dataset
# ./OrthogonalPlotter.exe $input $output $dataset
# ./VerticalOffsetVsTime.exe $input $output $dataset
# ./PlotVerticalPosition.exe $input $output $dataset
##########################################################
# Uncomment if re-submitting and you don't want to overwrite 

#if [[ ! -f $output ]]; then
#    ./Plotter.exe $input $output $dataset
#fi
##########################################################
