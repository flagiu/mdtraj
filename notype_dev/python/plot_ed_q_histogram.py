#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

outpng="ed_q_hist.png"
outpdf="ed_q_hist.pdf"

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Histogram of Eddington-Debenedettin "q" bond order parameter, of all particles at all times.',
                    epilog='End of the summary.\n'
)
parser.add_argument('--indat',  type=argparse.FileType('r'), required=False, default="ed_q.dat",
                     help="Input file: timestep | particle index | q. [default: %(default)s]"
)
parser.add_argument('--nbins', type=int,
                     default=20, required=False,
                     help="Number of bins for the histogram. [default: %(default)s]"
)
parser.add_argument('--i', type=int,
                     default=-1, required=False,
                     help="Plot only data of particle i; i can be 0,1,...,N-1 . Plot all if i<0. [default: %(default)s]"
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

header=args.indat.readline()
rcut1 = float( header.split("# cutoff =")[1].split(',')[0] )

x = np.loadtxt(args.indat)
indices = x[:,1]
N = int(np.max(indices)+1) # number of particles
n = len(x)//N # number of frames
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
x = x[n0*N:(n-n1)*N]

if args.i>=0:
    x = x[indices==args.i]

data = x[:,2]

fig, ax = plt.subplots(dpi=200)
ax.set_xlabel("q")
ax.set_ylabel("Density")
title=r"$r_{cut}=%.2f$ $\AA$"%rcut1
if args.i>=0:
    title += r", particle $i=%d$"%args.i
ax.set_title(title)
ax.hist(data, bins=args.nbins, ec = "black",  density=True)
#ax.set_yscale("log")
ax.tick_params(axis='y',which='both', direction='in')
#ax.grid(axis='both', which='major')
plt.tight_layout()

plt.savefig(outpng)
plt.savefig(outpdf)
print("%s: Figure saved on %s , %s\n"%(sys.argv[0],outpng, outpdf))
#plt.show()
