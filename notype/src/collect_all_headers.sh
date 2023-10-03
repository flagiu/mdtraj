
HEADER=all_headers.hpp
cat > $HEADER << EOF
#ifndef __ALL_HEADERS__
#define __ALL_HEADERS__

#include<cstdlib>
#include<cstring>
#include<complex>
#include<cmath>
#include "lib/utility.hpp"
#include "lib/vecflex.hpp"
#include "lib/particle.hpp"
#include "lib/mymatrix.hpp"
#include "lib/Ycomplex.hpp"
EOF
# First include external objects
for file in $(ls *.hpp | grep -v mdtraj.hpp)
do
    echo "#include \"$file\"" >> $HEADER
done
# Then MDTRAJ, then its input/output methods
for file in mdtraj.hpp io/*.hpp
do
    echo "#include \"$file\"" >> $HEADER
done

echo "#endif" >> $HEADER
