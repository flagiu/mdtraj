#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/mixture_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
echo "Compiling:"
make
echo
cd $back

gunzip trajectory.xyz.gz
echo "Test liquid g(r) for phase-change heterostructure:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz_cp2k trajectory.xyz -box6 30.073 -14.939 0.0 26.1 0.0 24.862 -rdf 0.02
printf "Exit $? \n\n"
# plot only the partial g(r) between Ge,Sb,Te (ignore Ti)
python3 ${PATH_TO_MDTRAJ}/python/plot_rdf_average.py --ignore 3
