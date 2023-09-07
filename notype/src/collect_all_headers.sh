
HEADER=all_headers.hpp
echo "#ifndef __ALL_HEADERS__" > $HEADER
echo "#define __ALL_HEADERS__" >> $HEADER
for file in mdtraj.hpp io/*.hpp statics/*.hpp dynamics/*.hpp
do
    echo "#include \"$file\"" >> $HEADER
done
echo "#endif" >> $HEADER
