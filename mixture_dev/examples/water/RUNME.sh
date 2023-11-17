#!/bin/bash
back=$(pwd)

echo $0": Unzipping trajectory:"
[ -e trajectory.gz ] && gunzip trajectory.gz

echo $0": Test g(r) & coordination number & E-D order parameter for H2O:"
../../bin/mdtraj -xdatcarV trajectory -rcut rcut.save -rdf 0.02 -cn -edq
printf "Exit $? \n\n"
python3 ../../python/plot_rdf_average.py

python3 ../../python/plot_coordnum_average.py
python3 ../../python/plot_coordnum_average_histogram.py
python3 ../../python/plot_coordnum_histogram.py

python3 ../../python/plot_ed_q_histogram.py
