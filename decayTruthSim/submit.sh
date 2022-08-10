stn=$1
nCores=$2

if [ $# -ne 2 ]; then
  echo "$# arguments supplied, station and cores, e.g. 'S12' 2."
  return
fi

config="dMu" #Bz" #"Bz"

if [ ! -d ${config} ]; then
    mkdir ${config}
fi

rm -f ${config}/*.root

dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/5.4e-18 
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/1.8e-18
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AcceptedDecays/5.4e-18
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/ControlSample/AcceptedDecays

for tree in `ls $dir | sort -V`; do
    id=${tree##*truthTrees.}
    id=${id%%.root*}
    echo "${dir}/${tree}" "${config}/plots_${id}_${stn}.root" "$stn"
    # ./Plotter.exe ${dir}/${tree} ${config}/plots_${id}_${stn}.root $stn
    # ./Acceptance.exe ${dir}/${tree} ${config}/plots_${id}.root
    # ./SmallAngleApprox.exe ${dir}/${tree} ${config}/plots_${id}.root
    # break
done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

# rm -f plots_${config}_${dataset}.root
# hadd -f plots_${config}_${stn}.root ${config}/plots_*_${stn}.root
# rm -f ${config}/plots_*_${dataset}.root

rm -f count_allDecays.root && hadd -f count_allDecays.root ${config}/plots_*.root
