#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/notype_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
echo "Compiling:"
make
echo
cd $back

echo "Test liquid RDF:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -rdf 256
printf "Exit $? \n\n"
printf " *** Does the liquid g(r) converge to 1 at large r? *** \n\n"
