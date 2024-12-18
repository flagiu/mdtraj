#!/bin/bash
root_path="/home/flavio/programmi/mdtraj/mixture"

#########################################################################################
# Analyze a jmd trajectory to distinguish crystal / non crystal and sc/a7  @ 1GPa, 400K #
#########################################################################################

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

echo "# Analyzing q_l for l=5,...,8 ..."
$MDTRAJ -bo -l 5 6 7 8 -v

echo "# Analyzing q_4 and crystalline clusters ..."
$MDTRAJ -bo -l 4 -v -qdot_th 0.65 -clusters rcut_clusters.dat -out_lammpsdump -pbc_out

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
echo "#   labeled trajectories: classified.dump classified_ave.dump"
echo "#   full remote paths:"
echo "#   sftp://$(whoami)@$(hostname -a | cut -d " " -f 1)/$(pwd)/clusters.l4.dump"
echo "#   sftp://$(whoami)@$(hostname -a | cut -d " " -f 1)/$(pwd)/classified.dump"
echo "#   sftp://$(whoami)@$(hostname -a | cut -d " " -f 1)/$(pwd)/classified_ave.dump"

rm tmp4
rm bo*.l*.dat bo*.l*.ave bo*.l*.local_ave

echo "# Plotting results ..."
python3 plot_nc_sc_a7.py
echo "#   plotted: sc_a7.png"
echo "#   quick visualization: xmgrace nc.l4.dat -nxy sc_a7.dat -nxy sc_a7_ave.dat"

