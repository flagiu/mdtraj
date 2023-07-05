#!/bin/bash

if(($#!=1)); then
	echo "Usage: $0 <.traj file from mdtraj>"
	echo "Convert block file into column file x,y1,y2,..."
	echo "Assuming the file has a 1-line header."
	echo "Data will be printed on screen."
	echo
	exit 1
fi

file=$1

nfiles=$(tail -n +2 $file | awk '/^$/ {N++;next}{print >"tmp"N}END{print N}')
nlines=$(cat tmp | wc -l)

touch accumulator.txt
ls tmp* | sort -V | xargs -d '\n' -n $(($(ulimit -n) - 10)) sh -c '
       paste accumulator.txt "$@" > accumulator.txt.sav;
       mv accumulator.txt.sav accumulator.txt
' _
cat accumulator.txt
rm tmp* accumulator.txt
