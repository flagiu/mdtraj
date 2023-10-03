#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/notype_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
echo "Compiling:"
bash install_path.sh
make
echo
cd $back

echo "Test Bond parameters:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -bo -l 4 -rcut1 1.5 -rcut2 2.5 -rcut3 4.0
printf "Exit $? \n\n"
exit

echo "Test liquid MSD:"
${PATH_TO_MDTRAJ}/bin/mdtraj -lammpstrj traj_long.lammpstrj.save -msd -period 10000
printf "Exit $? \n\n"
bash ../../shell/traj2nxy.sh msd.traj > msd.xxx
python3 ../../python/plot_msd_trajectory.py

echo "Test liquid S(q,t):"
${PATH_TO_MDTRAJ}/bin/mdtraj -lammpstrj traj_long.lammpstrj.save -sqt 2 50 1 -period 10000
printf "Exit $? \n\n"
python3 ../../python/plot_sqt.py

echo "Test liquid g(r):"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -rdf 0.02
printf "Exit $? \n\n"
bash ../../shell/traj2nxy.sh rdf.traj > rdf.xxx
python3 ../../python/plot_rdf_trajectory.py

echo "Test liquid S(q):"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -sq 2 100 1
printf "Exit $? \n\n"
bash ../../shell/traj2nxy.sh sq.traj > sq.xxx
python3 ../../python/plot_sq_trajectory.py

echo "Test S(q) from g(r):"
python3 ../../python/calc_sq_from_rdf.py --rho 1.02192 --Nq 50

echo "Test direct correlation:"
python3 ../../python/calc_direct_correlation_from_sq.py --rho 1.02192

echo "Test coordination number:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -cn -rcut1 1.5
printf "Exit $? \n\n"
python3 ../../python/plot_coordnum_histogram.py

echo "Test liquid ALTBC:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -altbc 0.01 0.85 15 -rcut1 2.40
printf "Exit $? \n\n"
python3 ../../python/plot_altbc.py

echo "Test liquid ADF (angular distribution):"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -adf 0.01 -rcut1 1.5
printf "Exit $? \n\n"
bash ../../shell/traj2nxy.sh adf.traj > adf.xxx
python3 ../../python/plot_adf_trajectory.py
python3 ../../python/plot_adf_cosine_trajectory.py
