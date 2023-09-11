#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np

outpng="coordnum.png"
outpdf="coordnum.pdf"

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Histogram of coordination number of all particles at all times.',
                    epilog='End of the summary.\n'
)
parser.add_argument('--indat',  type=argparse.FileType('r'), required=False, default="coordnum.dat",
                     help="Input file: timestep | particle index | coordination number. [default: %(default)s]"
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

x = np.loadtxt(args.indat)
n = len(x)
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
x = x[n0:n-n1]

cn = x[:,2]

plt.xlabel("Coordination Number")
plt.ylabel("Counts")
plt.hist(cn)

plt.savefig(outpng)
plt.savefig(outpdf)
#plt.show()
