#!/bin/bash

# pxy.dat contains the Pxy pressure component (in GPa)
# of a system of Sb atoms at 600K , 6.12g/cm3, NVE ensemble.
# 50001 steps with resolution 0.02 ps.

# The program "difflux1d" computes:
# - flux autocorrelation
# - flux-integral cross-correlation
# - integral mean-squared-displacement
# given a flux in a single-column file (here Pxy).

# we test different correlation windows "nt"
# and origin-intervals "oi"

delta=0.02 #ps

rm pxy*.out

for nt in 50 500 5000
do

for oi in 1 5 10
do

if(($oi>$nt/2)); then break ; fi

echo "<-- nt=${nt} , oi=${oi} -->"
printf -v out "pxy_nt%04d_oi%02d.out" ${nt} ${oi}
../../bin/difflux1d pxy.dat ${nt} ${oi} ${delta} > $out

done

done
echo "Completed"
