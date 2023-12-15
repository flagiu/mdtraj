import sys
from data import data
from xyz_cp2k import xyz

"""
#1 First argument must be a LAMMPS 'data' file.
#2 Second argument must be the prefix of the name of the output file (e.g. new --> new.xyz, new.box)
"""

d = data(sys.argv[1])
d.map(1,"id",2,"type",3,"x",4,"y",5,"z",6,"ix",7,"iy",8,"iz")
x = xyz(d)
box = x.one(sys.argv[2])

f = open(sys.argv[2]+".box","w")
for b in box[:-1]:
    f.write(f"{b} ")
f.write(f"{box[-1]}\n")
f.close()

print(f"\nBOX saved into {sys.argv[2]}.box\n")
