#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/notype_dev'
back=$(pwd)

cd $PATH_TO_MDTRAJ
echo "Compiling:"
make
echo
cd $back

echo "Test liquid RDF:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.liq.xyz.save -box6 8.8 -3.5 2.1 7.4 1.3 10.0 -rdf 256 -tag "liq"
printf "Exit $? \n\n"
printf " *** Does the liquid g(r) converge to 1 at large r? *** \n\n"

echo "Test crystallized RDF:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xtal.xyz.save -box6 8.8 -3.5 2.1 7.4 1.3 10.0 -rdf 256 -tag "xtal"
printf "Exit $? \n\n"

echo "Test crystallization with l=6 bond order correlation:"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz traj.xyz.save -box6 8.8 -3.5 2.1 7.4 1.3 10.0 -bo -l 6 -rcut1 1.5 -rcut2 2.4
printf "Exit $? \n\n"
