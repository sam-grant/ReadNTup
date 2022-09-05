input=$1
output=$2
truth=$3
dataset=$4 

# echo "Command: ./ThetaYvsMomentum.exe $input $output $dataset"
# ./ThetaYvsMomentum.exe $input $output $dataset
echo "Command: ./Plotter.exe $input $output $truth $dataset"
./Plotter.exe $input $output $truth $dataset