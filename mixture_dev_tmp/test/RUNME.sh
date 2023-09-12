#!/bin/bash
PATH_TO_MDTRAJ=$HOME'/programmi/mdtraj/mixture_dev'
back=$(pwd)

RC=5

cd $PATH_TO_MDTRAJ
make
cd $back

LLL=$(python3 ${PATH_TO_MDTRAJ}/test/generate_test.py)

echo "Test RDF"
${PATH_TO_MDTRAJ}/bin/mdtraj -xyz test.xyz -box3 $LLL -rdf 256
printf "Exit $? \n\n"
