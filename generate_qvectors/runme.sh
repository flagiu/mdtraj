#!/bin/bash
gcc -O3 -Wall ./generate_qvectors.c -lm -o ./a.x

mkdir QVECTORS
for i in {2..10..1} ; do
./a.x $i 2> tmp
done

rm ./a.x tmp
