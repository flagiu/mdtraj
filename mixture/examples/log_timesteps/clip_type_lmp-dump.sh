#!/bin/bash
TMP_FILE=tmp_awk_197619732
typemax=2
if(($#<1||$#>2)); then
	echo "Input: <lammps dump file> [typemax="$typemax"]"
	echo "I will clip any particle type >=typemax to typemax"
	echo "! Only for dumps having the format: ITEM: ATOMS type x y z"
	echo
	exit 1
fi
if(($#>=2)); then typemax=$2 ; fi

[ -e ${TMP_FILE} ] && echo "ERROR: temporary awk file exists! ${TMP_FILE}" && exit 1

cat > ${TMP_FILE} << EOF
{
if(\$0 ~ /ITEM: NUMBER OF ATOMS/){
  # print, read N from next line, and print next line
  print \$0;
  getline; N=\$1; print;
}
else if(\$0 ~ /ITEM: ATOMS type x y z/){
  # print, read atoms from next N lines, clip and print next N lines
  print \$0;
  for(i=0;i<N;i++){
    getline; type=\$1; x=\$2; y=\$3; z=\$4;
    if(type>typemax){ type=typemax; }; #clip type!
    printf "%d %.7f %.7f %.7f\n", type,x,y,z;
  }
}
else {
  # else just print
  print \$0;
}
}
EOF

awk -v typemax=$typemax -f ${TMP_FILE} $1

rm ${TMP_FILE}
