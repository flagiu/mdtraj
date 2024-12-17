#!/bin/bash
root_path="/home/flavio/programmi/mdtraj/mixture"

#########################################################################################
# Analyze a jmd trajectory to distinguish crystal / non crystal and sc/a7  @ 1GPa, 400K #
#########################################################################################
# Type assignment for LAMMPS:
# 1 liquid
# 2 xtal
# 3 sc
# 4 a7
# Type assignment for the scikit_learn classifier
# 0 sc
# 1 a7

# paste frames into a single trajectory
[[ -e traj.jmd ]] && rm traj.jmd
ls ../configurations/pos_* | sort -V | while read el; do cat $el >> traj.jmd; done

# prepare cutoffs
cat > rcut.dat << EOF
3.770
5.210
5.210
EOF

echo 3.90 > rcut_clusters.dat

MDTRAJ="${root_path}/bin/mdtraj -jmd traj.jmd -rcut rcut.dat"

echo "# Analyzing crystalline clusters with q4dot ..."
$MDTRAJ -bo -l 4 -v -qdot_th 0.65 -clusters rcut_clusters.dat -out_lammpsdump -pbc_out

echo "# Analyzing q_l for l=4,...,8 ..."
$MDTRAJ -bo -l 4 5 6 7 8 -v

# collect all q's into a single file
for l in {4..8}; do
awk '(NR>1){print $3}' boo.l${l}.dat > tmp$l
done
paste tmp{4..8} > q.dat
rm tmp{4..8}

# collect all locally-averaged q's into a single file
for l in {4..8}; do
awk '(NR>1){print $3}' boo_ave.l${l}.dat > tmp$l
done
paste tmp{4..8} > q_ave.dat
rm tmp{4..8}

# collect q4dot
awk '(NR>1){print $3}' boc.l4.dat > tmp4

echo "# Identifiying sc vs. a7 ..."
python classify.py classifier_sc_a7_l45678.pkl q.dat tmp4 clusters.l4.dump sc_a7.dat > classified.dump
python classify.py classifier_sc_a7_l45678_ave.pkl q_ave.dat tmp4 clusters.l4.dump sc_a7_ave.dat > classified_ave.dump

#rm tmp4
#rm bo*.l{4..22}.dat bo*.l{4..22}.ave bo*.l{4..22}.local_ave

