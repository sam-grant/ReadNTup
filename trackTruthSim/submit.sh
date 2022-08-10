nCores=$1 
truth=$2

if [ $# -ne 2 ]; then
  echo "$# arguments supplied, two required (cores and truth, e.g. . submit.sh 3 truth)"
  return
fi

config="dMu" #Bz" #"Bz"

if [ ! -d ${config} ]; then
    mkdir ${config}
fi

rm -f ${config}/*.root

# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/5.4e-18 
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/1.8e-18
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/AcceptedDecays/5.4e-18
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/ControlSample/AcceptedDecays
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/5.4e-18
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/1.8e-19
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/0mm/
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/1mm
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Plus1mm
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Minus1mm
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Plus0.1deg
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/Minus0.1deg

# compare these two!

dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Summer2022Campaign/sgrant/MainNtuple
# dir=/pnfs/GM2/persistent/EDM/MC/dMu/Trees/TrackRecoAndTruth/5.4e-18


for tree in `ls $dir | sort -V`; do

    if [[ $tree == attic ]]; then
        continue
    fi

    id=${tree##*truthTrees.}
    id=${id%%.root*}

    echo "${dir}/${tree}" "${config}/plots_${id}.root" "$truth"

done | xargs -i --max-procs=$nCores bash -c ". run.sh {}"

rm -f plots_${truth}.root && hadd -f plots_${truth}.root ${config}/plots_*.root