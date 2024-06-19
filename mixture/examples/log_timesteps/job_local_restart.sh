#!/bin/bash
# Running LAMMPS with a scheduled time for dump
# Potential: NN bernasconi-behler with no D2
# Integrator: NVE

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

LMP_TIMES="schedule.times"
LMP_IN="in.lmp"
LMP_OUT="out.lmp"
LMP_ERR="err.lmp"

mpirun -np 2 ~/LAMMPS_2017/src/lmp_mpi -in $LMP_IN >> $LMP_OUT 2>> $LMP_ERR
res=$?
if(($res!=0)); then echo "An error occurred!"; fi
exit $res
