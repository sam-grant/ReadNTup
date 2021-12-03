# file=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackReco/simTree.dMu.5.4e-18.root
# x2 /pnfs/GM2/persistent/EDM/MC/dMu/Trees/ControlSample/TrackReco/
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackReco
# ./Plotter.exe $file plots_${config}.root

config="dMu" 

if [ ! -d ${config} ]; then
    mkdir ${config}
fi

dir="/pnfs/GM2/persistent/EDM/MC/dMu/Trees/ControlSample/TrackReco"

for tree in `ls $dir | sort -V`; do
    id=${tree##*truthTrees.}
    id=${id%%.root*}
    ./Plotter.exe ${dir}/${tree} ${config}/plots_${id}.root
done

rm -f plots_${config}.root
hadd -f plots_${config}.root ${config}/plots_*.root
rm -f ${config}/plots_*.root

