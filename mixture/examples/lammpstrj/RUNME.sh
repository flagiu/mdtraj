#!/bin/bash

run_fraction_of_trajectory() {
tag=$1     # string tag for output files. Use the empty string ('' or "") for no tag
f0=$2      # fraction in [0,1]
f1=$3      # fraction in [0,1]
period=$4  # period, in timestep units

# change this if needed #
MDTRAJ_PATH="/home/flavio/programmi/mdtraj/mixture"
traj="../dump.lammpstrj"
#-----------------------#

MDTRAJ="${MDTRAJ_PATH}/bin/mdtraj -lammpstrj ${traj} -fskip $f0 $f1"
[[ $tag ]] && MDTRAJ="${MDTRAJ} -tag $tag" # if tag is not empty, add it as an argument
[[ $tag ]] && tagdot=".${tag}" || tagdot="" # if tag is not empty, prepend a '.' for output files' names

echo "|||||||| $tag : g(r) and S(q) ..."
$MDTRAJ -rdf 0.02 -1 -sq 2 100 1

python ${MDTRAJ_PATH}/python/find_rdf_local_minima.py rdf${tagdot}.ave rcut${tagdot}.dat hanning 5 0.8 1
python ${MDTRAJ_PATH}/python/find_sq_local_maxima.py sq${tagdot}.ave qmax${tagdot}.dat hanning 5 0.1 1
L=$(grep BOX -A 1 $traj | head -n 2 | tail -n 1 | awk '{printf "%.7f\n",$2-$1}')
qpeak=$(head -n 1 qmax${tagdot}.dat | awk -v L=$L 'BEGIN{pi=atan2(0,-1); dq=pi/L}{printf "%.0f", $1/dq-1}')

echo "|||| $tag ||||  S(q,t) and MSD(t) ..."
$MDTRAJ -sqt 2 $((${qpeak}+5)) 1 -msd -period $period
echo "|||||||| $tag : Prob(angle), q_Errington-Debenetetti ..."
$MDTRAJ -rcut rcut${tagdot}.dat -adf 0.01 -edq
echo "|||||||| $tag : Q4 and detailed Coordination-Number ..."
$MDTRAJ -rcut rcut${tagdot}.dat -bo -l 4 -cn -v
# higher cutoff for ALTBC and NND
awk '(NR==2){print $0}' rcut${tagdot}.dat > tmp
echo "|||||||| $tag : ALTBC and nearest-neighbours-distances ..."
$MDTRAJ -rcut tmp -altbc 0.02 2.5 25 -nnd 7 -v
rm tmp

MDTRAJ_PY=${MDTRAJ_PATH}/python
python ${MDTRAJ_PY}/plot_msd_average.py --file msd${tagdot}.ave --dt 0.002 --fitD True --inlabels labels${tagdot}.dat
python ${MDTRAJ_PY}/plot_sqt.py --inavg sqt${tagdot}.ave --dt 0.002 --fmt .- --normalize 1 --select_q $((${qpeak}-2))
python ${MDTRAJ_PY}/plot_adf_average.py --inavg adf${tagdot}.ave --inlabels labels${tagdot}.dat
python ${MDTRAJ_PY}/plot_altbc.py --inavg altbc${tagdot}_type0.ave
python ${MDTRAJ_PY}/plot_coordnum_histogram.py --indat coordnum${tagdot}.dat --inlabels labels${tagdot}.dat
python ${MDTRAJ_PY}/plot_nnd.py --indat nnd${tagdot}.dat --inlabels labels${tagdot}.dat --xlim 2.6 4.6
python ${MDTRAJ_PY}/plot_ed_q_histogram.py --indat ed_q${tagdot}.dat --inlabels labels${tagdot}.dat
[[ $tag ]] && (
mv msd.png msd${tagdot}.png
mv msd_D.dat msd_D${tagdot}.dat
mv msd_D.png msd_D${tagdot}.png
mv adf.png adf${tagdot}.png
mv altbc.png altbc${tagdot}.png
mv coordnum_hist.png coordnum_hist${tagdot}.png
mv nnd.png nnd${tagdot}.png
mv ed_q_hist.png ed_q_hist${tagdot}.png
)
echo "|||||||| $tag : Clearing PDFs and heaviest output files ..."
rm *.pdf log$tagdot nnd${tagdot}.dat traj.jmd coordnum${tagdot}.dat ed_q${tagdot}.dat
rm boo*${tagdot}.dat boc*${tagdot}.dat boo*${tagdot}.local_ave boc*${tagdot}.local_ave
}

run_fraction_of_trajectory "" 0.1 0.1 10000
