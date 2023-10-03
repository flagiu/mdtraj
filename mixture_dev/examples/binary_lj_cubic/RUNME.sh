#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/mixture_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
echo "Compiling:"
make
echo
cd $back

#mpirun lmp -in in.lmp > out.lmp

echo "Test liquid g(r) and coordination number for binary mixture:"
../../bin/mdtraj -lammpstrj dump.lammpstrj -fskip 0.25 0.0 -rdf 0.1 -cn -rcut1 3.0
printf "Exit $? \n\n"
python3 ../../python/plot_rdf_average.py
python3 ../../python/plot_coordnum_average_histogram.py

../../bin/mdtraj -lammpstrj dump.lammpstrj -cn -rcut1 3.0
python3 ../../python/plot_coordnum_average.py

