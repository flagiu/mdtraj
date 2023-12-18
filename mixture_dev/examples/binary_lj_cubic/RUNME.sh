#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/mixture_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
echo "RUNME.sh: Compiling:"
make
echo
cd $back

echo "RUNME.sh: Unzipping trajectory:"
gunzip dump.lammpstrj.gz
#mpirun lmp -in in.lmp > out.lmp

echo "RUNME.sh: Test liquid g(r), coordination number and bond order at equilibrium:"
../../bin/mdtraj -lammpstrj dump.lammpstrj -fskip 0.25 0.0 -rdf 0.1 -cn -rcut rcut.save -edq -bo -l 4
../../bin/mdtraj -lammpstrj dump.lammpstrj -fskip 0.25 0.0 -rcut rcut.save -bo -l 6
printf "Exit $? \n\n"
paste boo.l{4,6}.dat | awk '(NR>1){print $3,$6}' > q4q6.dat
paste boo_ave.l{4,6}.dat | awk '(NR>1){print $3,$6}' > q4q6_ave.dat
python3 ../../python/plot_rdf_average.py
python3 ../../python/plot_coordnum_average_histogram.py
python3 ../../python/plot_ed_q_histogram.py
python3 ../../python/plot_bond_order_ave_histogram.py --l 4
python3 ../../python/plot_bond_order_ave_histogram.py --l 6

echo "RUNME.sh: Test coordination number in transient:"
../../bin/mdtraj -lammpstrj dump.lammpstrj -fskip 0.0 0.5 -cn -rcut rcut.save
printf "Exit $? \n\n"
python3 ../../python/plot_coordnum_average.py

