#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/mixture_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
echo "Compiling:"
make
echo
cd $back

echo "Test liquid g(r) for binary mixture:"
${PATH_TO_MDTRAJ}/bin/mdtraj -lammpstrj dump.lammpstrj.save -rdf 0.02
printf "Exit $? \n\n"
python3 ${PATH_TO_MDTRAJ}/python/plot_rdf_average.py
