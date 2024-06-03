#!/bin/bash
back=$(pwd)

echo $0": Unzipping trajectory:"
[ -e trajectory.gz ] && gunzip trajectory.gz

echo $0": Test g(r) & coordination number & E-D order parameter for H2O:"
../../bin/mdtraj -xdatcarV trajectory -rcut rcut.save -sq 2 2 1 -sqt 2 2 1 -d -period 500 -fskip .99 0 -d
exit
../../bin/mdtraj -xdatcarV trajectory -rcut rcut.save -rdf 0.02 -1 -cn -edq -bo -l 4 -msd -period 500
../../bin/mdtraj -xdatcarV trajectory -rcut rcut.save -bo -l 6
printf "Exit $? \n\n"

awk '{print $1,$2,$4}' msd.ave > msdO.dat
awk '{print $1,$3,$5}' msd.ave > msdH.dat
awk '{print $1,$6}' msd.ave > msdCM.dat

python3 ../../python/plot_rdf_average.py

python3 ../../python/plot_coordnum_average.py
python3 ../../python/plot_coordnum_average_histogram.py
python3 ../../python/plot_coordnum_histogram.py

python3 ../../python/plot_ed_q_histogram.py
python ../../python/plot_msd_average.py
