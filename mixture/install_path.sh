#!/bin/sh
actual_path=$(pwd)
sed -i   's@root_path=.*@root_path='"\"$actual_path\""';@' src/mdtraj.hpp
sed -i   's@root_path=.*@root_path='"\"$actual_path\""'@' python/classify_sc_a7/classify_sc_a7.sh

actual_qpath=$(dirname `pwd`)/QVECTORS
sed -i   's@qvectors_path=.*@qvectors_path='"\"$actual_qpath\""';@' src/sq.hpp
sed -i   's@qvectors_path=.*@qvectors_path='"\"$actual_qpath\""';@' src/sqt.hpp
