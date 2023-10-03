#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/mixture_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
make
cd $back

LLL=$(python3 ./generate_xyz.py)

echo "Test g(r):"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz -box3 $LLL -rdf 0.02
printf "Exit $? \n\n"
python3 ${PATH_TO_MDTRAJ}/python/plot_rdf_average.py

echo "Test coordination number:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz -box3 $LLL -cn
printf "Exit $? \n\n"
python3 ${PATH_TO_MDTRAJ}/python/plot_coordnum_average.py
python3 ${PATH_TO_MDTRAJ}/python/plot_coordnum_average_histogram.py
python3 ${PATH_TO_MDTRAJ}/python/plot_coordnum_histogram.py
