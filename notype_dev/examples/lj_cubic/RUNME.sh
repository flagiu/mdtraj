#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/notype_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
echo "Compiling:"
make
echo
cd $back

#echo "Test liquid RDF:"
#${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -rdf 0.02
#printf "Exit $? \n\n"
#printf " *** Does the liquid g(r) converge to 1 at large r? *** \n\n"

echo "Test liquid S(q):"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -sq 2 100 1
printf "Exit $? \n\n"
exit
echo "Test liquid ALTBC:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box1 8.67 -altbc 0.01 0.85 25 -rcut1 2.25
printf "Exit $? \n\n"
