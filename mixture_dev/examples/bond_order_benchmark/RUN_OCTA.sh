#!/bin/bash
echo
echo "This program executes a benchmark of methods for computing the bond order"
echo " of an octahedral crystal structure."
echo "Methods are: Errington-Debenedetti, Spherical Harmonics with l=4,6"
echo
../../bin/mdtraj -poscar octa.vasp -rcut octa.rcut -rdf 0.02 8.0 -adf 0.01 -edq -bo -l 4 -image_convention -1
../../bin/mdtraj -poscar octa.vasp -rcut octa.rcut -rdf 0.02 8.0 -adf 0.01 -edq -bo -l 6 -image_convention -1

python3 ../../python/plot_rdf_average.py --yshift 10
python3 ../../python/plot_adf_average.py --yshift 10
python3 ../../python/plot_ed_q_histogram.py
../../python/plot_bond_order_ave_histogram.py --l 4
../../python/plot_bond_order_ave_histogram.py --l 6
