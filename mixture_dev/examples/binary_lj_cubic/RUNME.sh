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

echo "RUNME.sh: Test liquid g(r) and coordination number at equilibrium:"
../../bin/mdtraj -lammpstrj dump.lammpstrj -fskip 0.25 0.0 -rdf 0.1 -cn -rcut1 3.0
printf "Exit $? \n\n"
python3 ../../python/plot_rdf_average.py
python3 ../../python/plot_coordnum_average_histogram.py

echo "RUNME.sh: Test coordination number in transient:"
../../bin/mdtraj -lammpstrj dump.lammpstrj -fskip 0.0 0.5 -cn -rcut1 3.0
python3 ../../python/plot_coordnum_average.py

