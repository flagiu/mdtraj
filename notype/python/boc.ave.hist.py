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
parser.add_argument('--inavg',  type=argparse.FileType('r'),
                     default="boc.l4.ave", required=False,
                     help="Input file with average q_l^dot (columns:t,q_l^dot,q_l^dot error). [default: %(default)s]"
)
parser.add_argument('-l', type=int,
                     default=4, required=False,
                     help="Angular momentum value. [default: %(default)s]"
)
parser.add_argument('-ylog', type=bool,
                     default=False, required=False,
                     help="Use log y scale?. [default: %(default)s]"
)

x_tolerance=1e-5

#-------------------------------------#
args = parser.parse_args()
outpng="boc.l%d.png"%args.l
outpdf="boc.l%d.pdf"%args.l

Xa = np.loadtxt(args.inavg, unpack=False)
t = Xa[:,0]
q = Xa[:,1]
q_ = Xa[:,2]

fig, ax = plt.subplots(dpi=300)
ax.set_xlabel(r"$q_{%d}^{dot}$"%args.l)
ax.set_ylabel("counts")
if args.ylog:
	ax.set_yscale('log')
	ax.grid(axis='both', which='major')
	ax.grid(axis='y', which='minor', alpha=0.2)
else:
	ax.grid(axis='both', which='major')
ax.hist(q)
ax.set_xlim((-0.1,1.1))
ax.tick_params(which='both', direction='in')
plt.tight_layout()
fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()

