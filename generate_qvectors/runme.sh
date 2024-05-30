#!/bin/bash
gcc -O3 -Wall ./generate_qvectors.c -lm -o ./a.x

mkdir QVECTORS
for i in {2..300..1} ; do
./a.x $i 2> tmp
done

rm ./a.x tmp

exit

# this is to check that max wavenumber is <= q/2 for each q
for i in {2..300..1} ; do
printf -v q "%03d" $i
nmax=$(bash get_max_wavenumber.sh ../QVECTORS/qvector.$q)
if(($nmax>($i/2))); then
  echo $i $nmax
fi
done
