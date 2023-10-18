#!/bin/bash

~/programmi/mdtraj/notype_dev/bin/mdtraj -xyz_cp2k octa.xyz -box1 2 -rcut1 0.6 -adf 0.01 -tag octa
~/programmi/mdtraj/notype_dev/shell/traj2nxy.sh adf.octa.traj > adf.octa.xxx
~/programmi/mdtraj/mixture_dev/python/plot_adf_trajectory.py --intraj adf.octa.xxx --inavg adf.octa.ave

~/programmi/mdtraj/mixture_dev/bin/mdtraj -xyz_cp2k octa.xyz -box1 2 -rcut octa.rcut -edq -tag octa
~/programmi/mdtraj/mixture_dev/python/plot_ed_q_histogram.py --indat ed_q.octa.dat --inlabels labels.octa.dat
