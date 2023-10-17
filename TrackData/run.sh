input=$1
output=$2
dataset=$3 

##########################################################
# Uncomment if re-submitting and you don't want to overwrite 

#if [[ ! -f $output ]]; then
#    ./Plotter.exe $input $output $dataset
#fi
##########################################################

echo "Command: ./Plotter.out $input $output $dataset"
./Plotter.out $input $output $dataset