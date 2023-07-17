#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/notype_dev'
back=$(pwd)

RC=5

cd $PATH_TO_MDTRAJ
make
cd $back

LLL=$(python3 ${PATH_TO_MDTRAJ}/test/generate_test.py)

echo "Test Coordination Number"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz test.xyz -box3 $LLL -rcut1 $RC -cn
printf "Exit $? \n\n"

echo "Test RDF"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz test.xyz -box3 $LLL -rdf 200
printf "Exit $? \n\n"

echo "Test ADF"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz test.xyz -box3 $LLL -rcut1 $RC -adf 200
printf "Exit $? \n\n"

for l in 4 5 6
do
echo "Test Bond Orientational parameters with l=$l"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz test.xyz -box3 $LLL -l $l -bo -rcut1 $RC -rcut2 $RC -rcut3 $RC
printf "Exit $? \n\n"
done

echo "Test MSD"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz test.xyz -box3 $LLL -msd
printf "Exit $? \n\n"
