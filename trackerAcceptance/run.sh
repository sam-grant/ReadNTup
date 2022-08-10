<<<<<<< HEAD
input=$1
output=$2

./Plotter.exe $input $output
=======
config="plots" 

rm -rf $config && mkdir $config

dir=/pnfs/GM2/persistent/EDM/MC/TrackerAcceptance/jprice_mc_digit_trackAcceptance

for tree in `ls $dir | sort -V`; do
    
    id=${tree##*trackerAcceptanceTrees.}
    id=${id%%.root*}

# # 1567278068
#    if [ $id -eq 1567278068 ] || [ $id -eq 1567292467 ] || [ $id -eq 1567364467 ] || [ $id -eq 1567418467 ]; then
    # if [ $id -eq 1567292467 ] || [ $id -eq 1567364467 ] || [ $id -eq 1567418467 ]; then
#        echo ""
#        echo "${dir}/${tree} is corrupt, skipping..."
#        echo ""
#        continue
#    else
       echo "./Plotter.exe ${dir}/${tree} ${config}/trackerAcceptancePlots.${id}.root"
        ./Plotter.exe ${dir}/${tree} ${config}/trackerAcceptancePlots.${id}.root
#    fi

    # break
    
done

rm -f trackerAcceptancePlots.root && hadd -f trackerAcceptancePlots.root ${config}/trackerAcceptancePlots.*.root
>>>>>>> 61e4c6dfc811565a0487218eb445720b32e0a372
