input=$1
output=$2
stn=$3

echo "Command: ./ThetaYvsMomentum.exe $input $output $stn"
# ./Plotter.exe $input $output $stn
# echo "./CountEverything.exe $input"
# ./CountEverything.exe $input
# ./Count.exe $input $output $stn
./ThetaYvsMomentum.exe $input $output #$stn

# if [ ! -d ${config} ]; then
#     mkdir ${config}
# fi

# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/5.4e-18 
# #dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/1.8e-18
# # dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AcceptedDecays/5.4e-18
# # dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/ControlSample/AcceptedDecays

# for tree in `ls $dir | sort -V`; do
#     id=${tree##*truthTrees.}
#     id=${id%%.root*}
#     ./Plotter.exe ${dir}/${tree} ${config}/plots_${id}.root
#     #echo "./Plotter.exe ${dir}/${tree} ${config}/plots_${id}.root"
#     #./Acceptance.exe ${dir}/${tree} ${config}/plots_${id}.root
#     # ./SmallAngleApprox.exe ${dir}/${tree} ${config}/plots_${id}.root
#     # break
# done

# rm -f plots_${config}.root
# hadd -f plots_${config}.root ${config}/plots_*.root
# rm -f ${config}/plots_*.root

