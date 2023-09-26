#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/notype_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
echo "Compiling:"
make
echo
cd $back

echo "Test liquid g(r):"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -rdf 0.02
printf "Exit $? \n\n"
printf " *** Does the liquid g(r) converge to 1 at large r? *** \n\n"
bash ../../shell/traj2nxy.sh rdf.traj > rdf.xxx
python3 ../../python/plot_rdf_trajectory.py

echo "Test liquid S(q):"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -sq 2 100 1
printf "Exit $? \n\n"
bash ../../shell/traj2nxy.sh sq.traj > sq.xxx
python3 ../../python/plot_sq_trajectory.py

echo "Test S(q) rom g(r):"
python3 ../../python/calc_sq_from_rdf.py --rho 1.02192 --Nq 50

echo "Test direct correlation:"
python3 ../../python/calc_direct_correlation_from_sq.py --rho 1.02192

echo "Test coordination number:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -cn -rcut1 1.5
printf "Exit $? \n\n"
python3 ../../python/plot_coordnum_histogram.py

#echo "Test liquid ALTBC:"
#${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -altbc 0.01 0.85 15 -rcut1 2.40
#printf "Exit $? \n\n"
#python3 ../../python/plot_altbc.py
