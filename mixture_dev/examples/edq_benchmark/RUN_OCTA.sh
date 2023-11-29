#!/bin/bash

../../bin/mdtraj -poscar octa.vasp -rcut octa.rcut -rdf 0.02 8.0 -adf 0.01 -edq -image_convention -1

python3 ../../python/plot_rdf_average.py --yshift 10
python3 ../../python/plot_adf_average.py --yshift 10
python3 ../../python/plot_ed_q_histogram.py
