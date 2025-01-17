#!/bin/bash

run_fraction_of_trajectory() {
tag=$1     # string tag for output files. Use the empty string ('' or "") for no tag
f0=$2      # fraction in [0,1]
f1=$3      # fraction in [0,1]
period=$4  # period, in timestep units

traj="./dump.lammpstrj"
MDTRAJ_PATH="/home/flavio/programmi/mdtraj/mixture"

MDTRAJ="${MDTRAJ_PATH}/bin/mdtraj -lammpstrj $traj -logtime ./schedule.times -fskip $f0 $f1 "
[[ $tag ]] && MDTRAJ="${MDTRAJ} -tag $tag" # if tag is not empty, add it as an argument
[[ $tag ]] && tagdot=".${tag}" || tagdot="" # if tag is not empty, prepend a '.' for output files' names

#echo "|||| $tag ||||  g(r), S(q) ..."
$MDTRAJ -rdf 0.01 -1 -sq 2 100 1 || exit 1

python ${MDTRAJ_PATH}/python/find_sq_local_maxima.py sq${tagdot}.ave qmax${tagdot}.dat hanning 5 0.1 1
L=$(grep BOX -A 1 $traj | head -n 2 | tail -n 1 | awk '{printf "%.7f\n",$2-$1}')
qpeak=$(head -n 1 qmax${tagdot}.dat | awk -v L=$L 'BEGIN{pi=atan2(0,-1); dq=pi/L}{printf "%.0f", $1/dq-1}')

echo "|||| $tag ||||  S(q,t) and MSD(t) and Q(t) ..."
$MDTRAJ -sqt 2 $((${qpeak}+5)) 1 -msd || exit 1

python plot_Qoverlap.py --files Qoverlap${tagdot}.ave --dt 0.002

MDTRAJ_PY=${MDTRAJ_PATH}/python
python ${MDTRAJ_PY}/plot_msd_average.py --file msd${tagdot}.ave --dt 0.002 --fitD True --inlabels labels${tagdot}.dat
python ${MDTRAJ_PY}/plot_sqt.py --inavg sqt${tagdot}.ave --dt 0.002 --fmt .- --normalize 1 --select_q $((${qpeak}-2)) 
[[ $tag ]] && (
mv msd.png msd${tagdot}.png
mv msd_D.dat msd_D${tagdot}.dat
mv msd_D.png msd_D${tagdot}.png
mv sqt.png sqt${tagdot}.png
)
echo "|||| $tag |||| Clearing PDFs and heaviest output files ..."
rm *.pdf log$tagdot 2> tmp
rm tmp
}

run_fraction_of_trajectory "" 0.0 0.0


