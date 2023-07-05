#!/bin/python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Plots each trajectory r^2(t) and the average MSD=<r^2(t)>.',
                    epilog='End of the summary.'
)
parser.add_argument('--intraj', type=argparse.FileType('r'),
                     default="msd.xxx", required=False,
                     help="Input file with trajectories (columns:t,x1,x2,...). [default: %(default)s]"
)
parser.add_argument('--inavg',  type=argparse.FileType('r'),
                     default="msd.ave", required=False,
                     help="Input file with average trajectory (columns:t,MSD,MSD error). [default: %(default)s]"
)
parser.add_argument('--dt',  type=float, default=0.002, required=False,
                     help="Integration time step (picoseconds)"
)

outpng="msd.png"
outpdf="msd.pdf"

args = parser.parse_args()

Xt = np.loadtxt(args.intraj, unpack=False)
ntraj = Xt.shape[1]-1

Xa = np.loadtxt(args.inavg, unpack=False)
t = args.dt*Xa[:,0]
msd = Xa[:,1]

assert Xt.shape[0]==Xa.shape[0]
assert (Xt[:,0]==Xa[:,0]).all()
assert Xa.shape[1]>=2

fig, ax = plt.subplots()
ax.set_xlabel("t (ps)")
ax.set_ylabel(r"$\langle\, \Delta r^2(t)\, \rangle$ ($\AA^2$)")
ax.tick_params(which='both', direction='in')
for i in range(ntraj):
	ax.plot( t, Xt[:,i+1], 'k', alpha=0.1 )
ax.plot(t, msd, 'r')
ax.grid(axis='both', which='major')
plt.tight_layout()
fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng,outpdf))
#plt.show()

