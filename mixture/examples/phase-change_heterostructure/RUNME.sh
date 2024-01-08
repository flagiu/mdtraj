#!/bin/bash
back=$(pwd)

for fol in ../.. ../../../notype_dev
do
	cd $fol
	echo "RUNME.sh: Compiling:"
	bash install_path.sh
	make
	echo
	cd $back
done

echo $0": Unzipping trajectory:"
[ -e trajectory.xyz.gz ] && gunzip trajectory.xyz.gz

#echo $0": Test NOTYPE g(r) & coordination number & Eddington-Debenedetti order parameter for PCH"
#../../../notype_dev/bin/mdtraj -xyz_cp2k trajectory.xyz -box6 30.073 -14.939 0.0 26.1 0.0 24.862 -rcut rcut.save -rdf 0.02 -cn -edq
#printf "Exit $? \n\n"
#bash ../../../notype_dev/shell/traj2nxy.sh rdf.traj > rdf.xxx
#python3 ../../../notype_dev/python/plot_rdf_trajectory.py
#python3 ../../../notype_dev/python/plot_ed_q_histogram.py

echo $0": Test g(r) & coordination number & E-D order parameter for PCH:"
../../bin/mdtraj -xyz_cp2k trajectory.xyz -box6 30.073 -14.939 0.0 26.1 0.0 24.862 -rcut rcut.save -rdf 0.02 -cn -edq
printf "Exit $? \n\n"
# plot only the partial g(r) between Ge,Sb,Te (ignore Ti)
python3 ../../python/plot_rdf_average.py --ignore 3

python3 ../../python/plot_coordnum_average.py
python3 ../../python/plot_coordnum_average_histogram.py
python3 ../../python/plot_coordnum_histogram.py

python3 ../../python/plot_ed_q_histogram.py
