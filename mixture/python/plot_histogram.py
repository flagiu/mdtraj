#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Basic histogram of a column of x data.',
                    epilog='End of the summary.\n'
)
parser.add_argument('--file',  type=argparse.FileType('r'), required=True,
                     help="Input file."
)
parser.add_argument('--x', type=int,
                     default=0, required=False,
                     help="Index of the column of x data (0,1,2,...). [default: %(default)s]"
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

x = np.loadtxt(args.file)
n = len(x)
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
x = x[n0:n-n1]

plt.xlabel("x")
plt.ylabel("counts")
if len(x.shape)>1:
  plt.hist(x[:,args.x])
else:
  plt.hist(x)
plt.show()
