#!/bin/sh
actual_path=$(pwd)
sed -i   's@root_path=.*@root_path='"\"$actual_path\""';@' src/mdtraj.hpp

actual_qpath=$(dirname `pwd`)/QVECTORS
sed -i   's@qvectors_path=.*@qvectors_path='"\"$actual_qpath\""';@' src/sq.hpp
