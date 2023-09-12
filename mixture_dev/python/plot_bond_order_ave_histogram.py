#!/bin/python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Plots the histogram of BOCs q_l^dot among all frames.',
                    epilog='End of the summary.'
)
parser.add_argument('--l', type=int,
                     default=4, required=False,
                     help="Angular momentum value (for labels). [default: %(default)s]"
)
parser.add_argument('--inavg',  type=str,
                     default="NONE", required=False,
                     help="Input file with average q_l^dot (columns:t,q_l^dot,q_l^dot error). If NONE, it will be set to boc.lX.ave with X=l. [default: %(default)s]"
)
parser.add_argument('--scale', type=str,
                     default='lin', required=False,
                     help="Set y scale: 'lin' or 'log'. [default: %(default)s]"
)
parser.add_argument('--nbins', type=int,
                     default=20, required=False,
                     help="Number of bins for the histogram. [default: %(default)s]"
)
parser.add_argument('--qmin', type=float,
                     default=-0.1, required=False,
                     help="Minimum q value for the histogram and for the plot limits. [default: %(default)s]"
)
parser.add_argument('--qmax', type=float,
                     default=1.1, required=False,
                     help="Minimum q value for the histogram and for the plot limits. [default: %(default)s]"
)
x_tolerance=1e-5

#-------------------------------------#
args = parser.parse_args()
outpng="boc.l%d.png"%args.l
outpdf="boc.l%d.pdf"%args.l
if args.scale!='lin' and args.scale!='log':
    print("[ERROR: args.scale must be 'lin' or 'log']\n")
    sys.exit(1)
if args.inavg=="NONE":
	args.inavg = "boc.l%d.ave"%args.l

Xa = np.loadtxt(args.inavg, unpack=False)
t = Xa[:,0]
q = Xa[:,1]
q_ = Xa[:,2]

fig, ax = plt.subplots(dpi=300)
ax.set_xlabel(r"$q_{%d}^{dot}$"%args.l)
ax.set_ylabel("counts")
if args.scale=='log':
	ax.set_yscale('log')
	ax.grid(axis='both', which='major')
	ax.grid(axis='y', which='minor', alpha=0.2)
else:
	ax.grid(axis='both', which='major')
counts,bins = np.histogram(q, bins=args.nbins, range=(args.qmin,args.qmax) )
ax.stairs(counts,bins, baseline=0, fill=True)
ax.set_xlim((args.qmin, args.qmax))
ax.tick_params(which='both', direction='in')
plt.tight_layout()
fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()
