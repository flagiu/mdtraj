#!/bin/bash

traj="../dump.lammpstrj"
MDTRAJ_PATH="/home/flavio/programmi/mdtraj/mixture"
MDTRAJ_PY=${MDTRAJ_PATH}/python
MDTRAJ="${MDTRAJ_PATH}/bin/mdtraj -lammpstrj $traj -nodynamics"

cat > rcut.dat << EOF
3.58
5.34
5.34
EOF
echo 3.70 > rcut_clusters.dat
echo "|||| q4 ..."
$MDTRAJ -rcut rcut.dat -clusters rcut_clusters.dat -bo -l 4 -qdot_th 0.65 -out_lammpsdump -pbc_out

bash clip_type_lmp-dump.sh clusters.l4.dump 2 > clipped.dump

MDTRAJ="${MDTRAJ_PATH}/bin/mdtraj -lammpstrj clipped.dump -nodynamics -fskip 0.76 0 -tag A -dynamic_types"
echo 4.5 > tmp
$MDTRAJ -rcut tmp -altbc 0.02 2.5 25 -v

for i in $(seq 0 1)
do
python ${MDTRAJ_PY}/plot_altbc.py --inavg altbc.A_type${i}.ave --outname altbc.A_type${i}
done
rm tmp

