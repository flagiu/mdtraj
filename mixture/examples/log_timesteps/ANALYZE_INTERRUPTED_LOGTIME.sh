#!/bin/bash

run_fraction_of_trajectory() {
tag=$1     # string tag for output files. Use the empty string ('' or "") for no tag
f0=$2      # fraction in [0,1]
f1=$3      # fraction in [0,1]

traj="../dump.lammpstrj"
MDTRAJ_PATH="/home/flavio/programmi/mdtraj/mixture"

MDTRAJ="${MDTRAJ_PATH}/bin/mdtraj -lammpstrj $traj -nodynamics -fskip $f0 $f1 "
[[ $tag ]] && MDTRAJ="${MDTRAJ} -tag $tag" # if tag is not empty, add it as an argument
[[ $tag ]] && tagdot=".${tag}" || tagdot="" # if tag is not empty, prepend a '.' for output files' names

cat > rcut.dat << EOF
3.58
5.34
5.34
EOF
echo 3.70 > rcut_clusters.dat
echo "|||| $tag |||| q4 ..."
$MDTRAJ -rcut rcut.dat -clusters rcut_clusters.dat -cn -bo -l 4 -qdot_th 0.65

echo "|||| $tag |||| Clearing PDFs and heaviest output files ..."
rm *.pdf log$tagdot nnd${tagdot}.dat coordnum${tagdot}.dat ed_q${tagdot}.dat 2> tmp
rm boo*${tagdot}.dat boc*${tagdot}.dat boo*${tagdot}.local_ave boc*${tagdot}.local_ave 2> tmp
rm tmp
}

run_fraction_of_trajectory "" 0.0 0.0

