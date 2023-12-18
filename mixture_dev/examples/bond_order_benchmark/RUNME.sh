#!/bin/bash
LAMMPS_EXE=~/LAMMPS_2017/src/lmp_mpi
echo
echo "This program executes a benchmark of methods for computing the bond order of the given crystal structure."
echo "Methods are: Errington-Debenedetti, Spherical Harmonics with l=4,6"
echo "You need a LAMMPS executable to be specified in this file. Default one is: "$LAMMPS_EXE
echo

if(($#!=1)); then
  printf "\nUSAGE: $0 <sc,fcc,diamond>\n\n"
  printf " - sc with a=1, Lx=3, rcut=1.01: 27 perfect octahedral environments (Nc=6, q=0) for each frame;\n"
  printf " - sc with a=1, Lx=2, rcut=1.01: 8 3-fold defective pyramidal octahedral environments (Nc=3, q=7/8=0.875) for each frame;\n"
  printf " - fcc is not of interest;\n"
  printf " - diamond with a=1, Lx=1, rcut=0.44 (>~sqrt(3)/4): 8 perfect tetrahedral environments (Nc=4, q=1) for each frame.\n\n"
  exit 1
fi

key=$1
[ "$key" == "sc" ] && echo "KEY is sc"
[ "$key" == "fcc" ] && echo "KEY is fcc"
[ "$key" == "diamond" ] && echo "KEY is diamond"

# grep the radial cutoff for 1st shell neighbours
if [ "$key" == "sc" ]; then
  echo '{print $NF}' > awk_rcut
elif [ "$key" == "fcc" ]; then
  echo '{print $NF/sqrt(2)}' > awk_rcut
elif [ "$key" == "diamond" ]; then
  echo '{print $NF * sqrt(3)/4}' > awk_rcut
else
  echo "[ ERROR: input '$key' not accepted. ]\n"
  exit 1
fi
rcut1=$(grep "variable a equal" $key.lmp | awk -f awk_rcut)
echo $rcut1 > rcut.dat   #1st sphere
echo $rcut1 >> rcut.dat  #2nd sphere
echo $rcut1 >> rcut.dat  #3rd sphere
rm awk_rcut
printf "\nRCUT set to $rcut1\n\n"

# generate the .xyz and the .box
mpirun $LAMMPS_EXE < $key.lmp > $key.out
python3 data2xyz.py $key.data $key

# fake a second timestep
#cp $key.xyz tmp
#sed 's@i=\t0@i=\t1@g' $key.xyz >> tmp
#mv tmp $key.xyz

# grep the box
box=$(cat $key.box)

# ANALYSIS
../../bin/mdtraj -xyz_cp2k $key.xyz -box6 $box -rcut rcut.dat -rdf 0.02 -1 -adf 0.01 -cn -edq -bo -l 4 -tag $key
../../bin/mdtraj -xyz_cp2k $key.xyz -box6 $box -rcut rcut.dat -bo -l 6 -tag $key

for el in adf rdf; do
  python3 ../../python/plot_${el}_average.py --inavg ${el}.${key}.ave  --inlabels labels.${key}.dat --yshift 10
  for fmt in png pdf; do
    mv ${el}.${fmt} ${el}.${key}.${fmt}
  done
done

for el in coordnum ed_q; do
  python3 ../../python/plot_${el}_histogram.py --indat ${el}.${key}.dat --inlabels labels.${key}.dat
  for fmt in png pdf; do
    mv ${el}_hist.${fmt} ${el}_hist.${key}.${fmt}
  done
done
