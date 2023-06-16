#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj-dev'
RC=5

LLL=$(python3 generate_test.py)

echo "Test Coordination Number"
${PATH_TO_MDTRAJ}/bin/mdtraj -in test.xyz -L $LLL -rcut1 $RC -cn
printf "Exit $? \n\n"

echo "Test RDF"
${PATH_TO_MDTRAJ}/bin/mdtraj -in test.xyz -L $LLL -rdf 200
printf "Exit $? \n\n"

echo "Test ADF"
${PATH_TO_MDTRAJ}/bin/mdtraj -in test.xyz -L $LLL -rcut1 $RC -adf 200
printf "Exit $? \n\n"

for l in 4 5 6
do
echo "Test Bond Orientational parameters with l=$l"
${PATH_TO_MDTRAJ}/bin/mdtraj -in test.xyz -L $LLL -l $l -bo -rcut1 $RC -rcut2 $RC -rcut3 $RC
printf "Exit $? \n\n"
done

echo "Test MSD"
${PATH_TO_MDTRAJ}/bin/mdtraj -in test.xyz -L $LLL -msd 
printf "Exit $? \n\n"

