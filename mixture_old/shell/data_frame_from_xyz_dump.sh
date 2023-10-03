#!/bin/bash

if(($#!=2)); then
	echo "Usage: $0 <xyz dump file> <time fraction>"
	echo "Data will be printed on screen"
	echo
	exit 1
fi

filexyz=$1
frac=$2

# Grep N and the xyz coordinates
N=$(head -n 1 $filexyz)
nlines=$(cat $filexyz | wc -l)
nframes=$( echo "($nlines -1)/($N+2)" | bc )
iframe=$(echo "($frac * $nframes)/1" | bc )
iline=$(echo "2 + $iframe * ($N+2)" | bc )
tail $filexyz -n+$iline | head -n $(($N+1)) > tmp.xyz

# Build the outfile
COMMENT=$(head -n 1 tmp.xyz) # comment line
echo "Frame from "$filexyz". "$COMMENT
cat << EOF
$N atoms
1 atom types

0.0000000000000000e+00 1.5090000000000000e+01 xlo xhi
0.0000000000000000e+00 1.5680000000000000e+01 ylo yhi
0.0000000000000000e+00 1.5800000000000001e+01 zlo zhi

Masses

1 121.75

Atoms # atomic

EOF

tail -n+2 tmp.xyz | awk '{print NR,$0}'		# awk adds the atom number in increasing order
rm tmp.xyz
