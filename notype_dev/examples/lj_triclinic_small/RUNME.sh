#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/notype_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
echo "Compiling:"
make
echo
cd $back

echo "Test liquid RDF:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box6 8.8 -1.7 1.0 7.4 0.6 10.0 -rdf 256
printf "Exit $? \n\n"
printf " *** Does the liquid g(r) converge to 1 at large r? *** \n\n"
