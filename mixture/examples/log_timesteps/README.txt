# python generate_times.py
# bash job_local.sh
# bash ANALYZE.sh

- Logarithmic timesteps must be generated via generate_times.py ;

- An example of how to generate a LAMMPS trajectory according to the
  wanted timesteps is given in lammps.in , job_local.sh , job_local_restart.sh ,
  including an automatic restart procedure.

- MDtraj can only analyze one or more *full* cycles; the first cycle must start
  from the lowest time t=delta ;

- Static observables (g(r),S(q)) will be averaged one every logarithmic cycle;
- Dynamic observables (S(q,t),MSD(t)) will be averaged both within and outside
  each logarithmic cycle;
