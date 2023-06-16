#!/bin/bash

if(($#!=2)); then
	echo "Usage: $0 <LAMMPS data file> <target n. of atoms>"
	echo "Removes atoms from the data file in order to get to the given number."
	echo "Data will be printed on screen, with zero velocity."
	echo
	exit 1
fi

file=$1
Nf=$2

nlines=$(cat $file | wc -l)
# initial number of atoms (15 header lines, 3 lines separating positions and velocities)
Ni=$(echo "($nlines - 15 - 3)/2" | bc)
# atoms to be removed
DN=$(echo "$Ni - $Nf" | bc)
# Keep the first M atoms, and remove one every 2 for the rest
M=$(echo "$Ni - 2*$DN" | bc)

head -n 15 $file > header
sed -i 's/'"$Ni"'/'"$Nf"'/g' header

# Print file

echo "LAMMPS data file. Modified by $0"
tail -n+2 header
tail -n+16 $file | awk -v Nf=$Nf -v M=$M 'BEGIN{i=1}\
{if((NR<=M)||((NR%2==0)&&(i<=Nf))){print i,$2,$3,$4,$5,$6,$7,$8; i++}}\
END{print "\nVelocities\n"; for(i=1;i<=Nf;i++){ print i,0.0,0.0,0.0}}'

rm header
