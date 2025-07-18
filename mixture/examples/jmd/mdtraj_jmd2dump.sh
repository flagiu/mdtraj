#!/bin/bash

CONFSpath=$1

run_fraction_of_trajectory() {
tag=$1     # string tag for output files. Use the empty string ('' or "") for no tag
f0=$2      # fraction in [0,1]
f1=$3      # fraction in [0,1]
CONFSpath=$4
typeFile=`dirname $CONFSpath`"/type.dat"

MDTRAJ_PATH="/home/flavio/programmi/mdtraj/mixture"

[[ -e traj.jmd ]] && rm traj.jmd
ls ${CONFSpath}/pos_* | sort -V | while read el; do cat $el >> traj.jmd; done
MDTRAJ="${MDTRAJ_PATH}/bin/mdtraj -jmd traj.jmd $typeFile -fskip $f0 $f1"
[[ $tag ]] && MDTRAJ="${MDTRAJ} -tag $tag" # if tag is not empty, add it as an argument
[[ $tag ]] && tagdot=".${tag}" || tagdot="" # if tag is not empty, prepend a '.' for output files' names

echo $typeFile
$MDTRAJ  -out_lammpsdump -pbc_out
}

run_fraction_of_trajectory "" 0.0 0.0 "$CONFSpath"

rm log labels.dat ndens.dat traj.jmd 
