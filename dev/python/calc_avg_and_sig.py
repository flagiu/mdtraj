#!/usr/bin/env python3
import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Prints average and st.dev. of the given column of a file.',
                    epilog='End of the summary.'
)
parser.add_argument('--file',  type=argparse.FileType('r'), required=True,
                     help="Input file. Must have at least NCOL+1 columns."
)
parser.add_argument('--ncol', type=int,
                     default=0, required=False,
                     help="Index of the column to be analyzed (0,1,2,...). [default: %(default)s]"
)
parser.add_argument('--fskip0', type=float,
                     default=0.0, required=False,
                     help="Fraction of data to be skipped from the start. [default: %(default)s]"
)
parser.add_argument('--fskip1', type=float,
                     default=0.0, required=False,
                     help="Fraction of data to be skipped from the end. [default: %(default)s]"
)

args = parser.parse_args()

lines=args.file.readlines() #file is automatically opened
n=len(lines)
n0=int(args.fskip0*n)
n1=int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
sublines=lines[n0:n-n1]
s=0.0
s2=0.0
cc=0
for ll in lines:
    cols = ll.strip('\n').split(' ')
    assert len(cols)>args.ncol # the file must contain at least NCOL+1 columns!
    val = float(cols[args.ncol])
    s+=val
    s2+=val*val
    cc += 1
print(s/float(cc),np.sqrt(s2/float(cc)-(s/float(cc))**2))
