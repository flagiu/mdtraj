
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

for file in rdf.hpp sq.hpp msd_and_ngp.hpp mdtraj.hpp io/*.hpp statics/*.hpp
do
    echo "#include \"$file\"" >> $HEADER
done
echo "#endif" >> $HEADER
