#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/mixture_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
echo "RUNME.sh: Compiling:"
bash install_path.sh
make
echo
cd $back

echo "RUNME.sh: Unzipping trajectory:"
[ -e trajectory.xyz.gz ] && gunzip trajectory.xyz.gz

echo "RUNME.sh: Test g(r) & coordination number for phase-change heterostructure:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz_cp2k trajectory.xyz -box6 30.073 -14.939 0.0 26.1 0.0 24.862 -cn -rcut1 3.75 -rdf 0.02
printf "Exit $? \n\n"
python3 ${PATH_TO_MDTRAJ}/python/plot_coordnum_average.py
python3 ${PATH_TO_MDTRAJ}/python/plot_coordnum_average_histogram.py
python3 ${PATH_TO_MDTRAJ}/python/plot_coordnum_histogram.py
# plot only the partial g(r) between Ge,Sb,Te (ignore Ti)
python3 ${PATH_TO_MDTRAJ}/python/plot_rdf_average.py --ignore 3
