import sys
import argparse
import numpy as np
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Prints the coefficients of a linear fit (with no intercept) of the given column of a file.',
                    epilog='Printing format: slope | slope error.\n'
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
parser.add_argument('--dy', type=int,
                     default=-1, required=False,
                     help="Index of the column of y error bars. (0,1,2,...) If negative, don't use error bars. [default: %(default)s]"
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

if args.dy >=0:
	x, y, dy = np.loadtxt(args.file, unpack=True, usecols=(args.x,args.y,args.dy))
else:
	x, y = np.loadtxt(args.file, unpack=True, usecols=(args.x,args.y))

def f(x,a=1):
	return a*x

n = len(x)
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!

x = x[n0:n-n1]
y = y[n0:n-n1]
if args.dy>=0:
	dy = dy[n0:n-n1]
	popt,pcov = curve_fit(f, x,y, sigma=dy, p0=[1], maxfev=10000)
else:
	popt,pcov = curve_fit(f, x,y, sigma=None, p0=[1], maxfev=10000)

popt_ = np.sqrt(np.diag(pcov))

for a,a_ in zip(popt,popt_):
	print(a,a_)
