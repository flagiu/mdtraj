#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/mixture_dev'
back=$(pwd)

RC=5

cd $PATH_TO_MDTRAJ
make
cd $back

LLL=$(python3 ${PATH_TO_MDTRAJ}/test/generate_test.py)

echo "Test g(r):"
${PATH_TO_MDTRAJ}/bin/mdtraj -lammpstrj dump.lammpstrj -rdf 0.04
printf "Exit $? \n\n"
python3 ${PATH_TO_MDTRAJ}/python/plot_rdf_average.py
