#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/notype_dev'
back=$(pwd)

filetype=lammpstrj
file=dump.lammpstrj

cd $PATH_TO_MDTRAJ
echo "Compiling:"
make
echo
cd $back

echo "Test g(r) and coordination number:"
${PATH_TO_MDTRAJ}/bin/mdtraj -$filetype $file -rcut1 1.45 -cn -rdf 0.02
python3 ../../python/plot_coordnum_histogram.py
printf "Exit $? \n\n"
printf " *** Does the liquid g(r) converge to 1 at large r? *** \n\n"

echo "Test g_3(r):"
${PATH_TO_MDTRAJ}/bin/mdtraj -$filetype $file -altbc 0.01 0.85 15 -rcut1 2.40
python3 ../../python/plot_altbc.py
printf "Exit $? \n\n"
