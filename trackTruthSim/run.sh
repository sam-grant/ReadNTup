input=$1
output=$2
truth=$3

echo "Command: ./ThetaYvsMomentum.exe $input $output $truth"
./ThetaYvsMomentum.exe $input $output $truth # $dataset

# echo "Command: ./Plotter.exe $input $output $truth"
# ./Plotter.exe $input $output $truth 

# echo "Command: ./CountEverything.exe $input"
# ./CountEverything.exe $input 
