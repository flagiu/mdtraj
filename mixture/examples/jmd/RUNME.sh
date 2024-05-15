#!/bin/bash

run_fraction_of_trajectory() {
tag=$1 # string tag for output files. Use the empty string ('' or "") for no tag
f0=$2  # fraction in [0,1]
f1=$3  # fraction in [0,1]

MDTRAJ_PATH="/home/flavio/programmi/mdtraj/mixture"

[[ -e traj.jmd ]] && rm traj.jmd
ls ../configurations/pos_* | sort -V | while read el; do cat $el >> traj.jmd; done
MDTRAJ="${MDTRAJ_PATH}/bin/mdtraj -jmd traj.jmd -fskip $f0 $f1"
[[ $tag ]] && MDTRAJ="${MDTRAJ} -tag $tag" # if tag is not empty, add it as an argument
[[ $tag ]] && tagdot=".${tag}" || tagdot="" # if tag is not empty, prepend a '.' for output files' names

echo "|||||||| $tag : g(r) and S(q) ..."
$MDTRAJ -rdf 0.02 -1 -sq 2 100 1 -out_xyz -pbc_out

python ${MDTRAJ_PATH}/python/find_rdf_local_minima.py rdf${tagdot}.ave rcut${tagdot}.dat hanning 5 0.8 0

echo "|||||||| $tag : Prob(angle), q_Errington-Debenetetti ..."
$MDTRAJ -rcut rcut${tagdot}.dat -adf 0.01 -edq
echo "|||||||| $tag : Q4 and detailed Coordination-Number ..."
$MDTRAJ -rcut rcut${tagdot}.dat -bo -l 4 -cn -v
echo "|||||||| $tag : Q6 ..."
$MDTRAJ -rcut rcut${tagdot}.dat -bo -l 6 -v
# higher cutoff for ALTBC and NND
awk '(NR==2){print $0}' rcut${tagdot}.dat > tmp
echo "|||||||| $tag : ALTBC and nearest-neighbours-distances ..."
$MDTRAJ -rcut tmp -altbc 0.02 2.5 25 -nnd 7
rm tmp

awk '(NR>1){print $2}' boo.l4${tagdot}.ave > tmp1
awk '(NR>1){print $2}' boo.l6${tagdot}.ave > tmp2
paste tmp1 tmp2 > q4q6${tagdot}.ave
rm tmp1 tmp2

awk '(NR>1){print $2}' boo_ave.l4${tagdot}.ave > tmp1
awk '(NR>1){print $2}' boo_ave.l6${tagdot}.ave > tmp2
paste tmp1 tmp2 > q4q6_ave${tagdot}.ave
rm tmp1 tmp2

MDTRAJ_PY=${MDTRAJ_PATH}/python
python ${MDTRAJ_PY}/plot_adf_average.py --inavg adf${tagdot}.ave --inlabels labels${tagdot}.dat
python ${MDTRAJ_PY}/plot_altbc.py --inavg altbc${tagdot}.ave
python ${MDTRAJ_PY}/plot_coordnum_histogram.py --indat coordnum${tagdot}.dat --inlabels labels${tagdot}.dat
python ${MDTRAJ_PY}/plot_nnd.py --indat nnd${tagdot}.dat --inlabels labels${tagdot}.dat --xlim 2.6 4.6
python ${MDTRAJ_PY}/plot_ed_q_histogram.py --indat ed_q${tagdot}.dat --inlabels labels${tagdot}.dat
[[ $tag ]] && (
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

run_fraction_of_trajectory "" 0.01 0.0
