#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Basic line plot of x-y data, eventually with error bars on x and/or y.',
                    epilog='End of the summary.\n'
)
parser.add_argument('--file',  type=argparse.FileType('r'), required=True,
                     help="Input file."
)
parser.add_argument('--x', type=int,
                     default=0, required=False,
                     help="Index of the column of x data (0,1,2,...). [default: %(default)s]"
)
parser.add_argument('--y', type=int,
                     default=1, required=False,
                     help="Index of the column of y data (0,1,2,...). [default: %(default)s]"
)
parser.add_argument('--dx', type=int,
                     default=-1, required=False,
                     help="Index of the column of x error bars. (0,1,2,...) If negative, don't use x error bars. [default: %(default)s]"
)
parser.add_argument('--dy', type=int,
                     default=-1, required=False,
                     help="Index of the column of y error bars. (0,1,2,...) If negative, don't use y error bars. [default: %(default)s]"
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

plt.xlabel('x')
plt.ylabel('y')
fmt='x-'

if args.dx>=0 and args.dy>=0:
    x, y, dx, dy = np.loadtxt(args.file, unpack=True, usecols=(args.x,args.y,args.dx,args.dy), comments=['#'])
elif args.dx>=0:
    x, y, dx = np.loadtxt(args.file, unpack=True, usecols=(args.x,args.y,args.dx), comments=['#'])
elif args.dy>=0:
    x, y, dy = np.loadtxt(args.file, unpack=True, usecols=(args.x,args.y,args.dy),comments=['#'])
else:
    x, y = np.loadtxt(args.file, unpack=True, usecols=(args.x,args.y))

n = len(x)
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
x = x[n0:n-n1]
y = y[n0:n-n1]

if args.dx>=0 and args.dy>=0:
    dx = dx[n0:n-n1]
    dy = dy[n0:n-n1]
    plt.errorbar(x, y, yerr=dy, xerr=dx, fmt=fmt)
elif args.dx>=0:
    dx = dx[n0:n-n1]
    plt.errorbar(x, y, xerr=dx, fmt=fmt)
elif args.dy>=0:
    dy = dy[n0:n-n1]
    plt.errorbar(x, y, yerr=dy, fmt=fmt)
else:
    plt.plot(x, y, fmt)

plt.show()
