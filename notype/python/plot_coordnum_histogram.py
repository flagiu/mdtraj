#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

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
parser.add_argument('--logScale', type=bool,
                     default=False, required=False,
                     help="Use log scale for y axis?. [default: %(default)s]"
)
parser.add_argument('--outname',  type=str,
                     default="coordnum_hist", required=False,
                     help="Prefix for the output file. [default: %(default)s]"
)

args = parser.parse_args()

header=args.indat.readline()
rcut1 = float( header.split("# cutoff =")[1].split(',')[0] )

x = np.loadtxt(args.indat)
n = len(x)
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
x = x[n0:n-n1]

data = x[:,2]
cn_u = np.unique(data)
d = 1 if len(cn_u)==1 else np.diff(cn_u).min()
left_of_first_bin = data.min() - float(d)/2
right_of_last_bin = data.max() + float(d)/2
bins = np.arange(left_of_first_bin, right_of_last_bin + d, d)

fig, ax = plt.subplots(dpi=200)
ax.set_xlabel("Coordination Number")
ax.set_ylabel("Counts")
ax.set_title(r"$r_{cut}=%.2f$ $\AA$"%rcut1)
ax.hist(data, bins=bins, ec = "black", rwidth=0.8, density=True)
ax.set_ylabel("Density")
if args.logScale:
    ax.set_yscale("log")
ax.tick_params(axis='y',which='both', direction='in')
ax.grid(axis='both', which='major')
plt.tight_layout()

fig.savefig(args.outname+".pdf")
fig.savefig(args.outname+".png")
print("Figure saved on %s.png , %s.pdf\n"%(args.outname,args.outname))
#plt.show()
