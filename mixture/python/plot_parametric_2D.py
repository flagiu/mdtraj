#!/usr/bin/python3
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Basic parametric plot of x(t) vs y(t) data.',
                    epilog='End of the summary.\n'
)
parser.add_argument('fileX',  type=str, help="Input file for x(t) (>=2 columns).")
parser.add_argument('fileY',  type=str, help="Input file for y(t) (>=2 columns).")
parser.add_argument('--x', type=int, default=1, required=False,
                     help="Index of the column of x data (1,2,...). [default: %(default)s]"
)
parser.add_argument('--y', type=int, default=1, required=False,
                     help="Index of the column of y data (1,2,...). [default: %(default)s]"
)

parser.add_argument('--fskip0', type=float, default=0.0, required=False,
                     help="Fraction of data to be skipped from the start. [default: %(default)s]"
)
parser.add_argument('--fskip1', type=float, default=0.0, required=False,
                     help="Fraction of data to be skipped from the end. [default: %(default)s]"
)

args = parser.parse_args()

plt.xlabel('x')
plt.ylabel('y')

t,x  = np.loadtxt(args.fileX, unpack=True, usecols=(0,args.x), comments="#")
ty,y = np.loadtxt(args.fileY, unpack=True, usecols=(0,args.y), comments="#")

assert (t==ty).all()

n = len(x)
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
t = t[n0:n-n1]
x = x[n0:n-n1]
y = y[n0:n-n1]

#plt.plot(x, y, "k-", alpha=0.7)
plt.scatter(x, y, marker="x", c=t, alpha=0.3)
plt.colorbar().set_label("t")

plt.show()
