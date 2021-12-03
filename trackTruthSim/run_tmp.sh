config="dMu" #Bz" #"Bz"

if [ ! -d ${config} ]; then
    mkdir ${config}
fi


dir="/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackTruth" 

for tree in `ls $dir | sort -V`; do
    id=${tree##*truthTrees.}
    id=${id%%.root*}
    ./Plotter.exe ${dir}/${tree} ${config}/plots_${id}.root
done

rm -f plots_${config}.root
hadd -f plots_${config}.root ${config}/plots_*.root
rm -f ${config}/plots_*.root
