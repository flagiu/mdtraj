#!/bin/python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid, simpson
#plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Computes c(r) from S(q using Ornstein-Zernicke equation and discrete Fourier transform, and plots it.',
                    epilog='End of the summary.'
)
parser.add_argument('--sq',  type=argparse.FileType('r'),
                     default="sq.ave", required=False,
                     help="Input file with average S(q) (columns:q,S(q),...). [default: %(default)s]"
)
parser.add_argument('--rho',  type=float, required=True,
                     help="Number density."
)
parser.add_argument('--Nr',  type=int,
                     default=100, required=False,
                     help="Number of r points. [default: %(default)s]"
)

outpng="cr_from_sq.png"
outpdf="cr_from_sq.pdf"

#-------------------------------------#
args = parser.parse_args()

Xa = np.loadtxt(args.sq, unpack=False)
q = Xa[:,0]
sq = Xa[:,1]

dq = q[1]-q[0]
qmax = q[-1]
dr = 2*np.pi/(2*qmax) # minimum frequency

r = dr * np.arange(1,1+args.Nr)
cr = np.empty(len(r))
cq = (1-1/sq)/ args.rho
f = cq / (2*np.pi)**3 # function to be integrated

for i in range(len(r)):
    integrand = 4*np.pi * q * f * np.sin(r[i]*q)/r[i]
    integral = simpson(integrand, q)
    cr[i] = integral
np.savetxt("cr_from_sq.ave", np.array([r,cr]).T, header="# r | c(r)")

fig, axes = plt.subplots(2,2, dpi=300)

ax = axes[0][0]
ax.set_xlabel(r"$q$ $[\AA^{-1}]$")
ax.set_ylabel(r"$c(q) = [1-1/S(q)]/\rho$")
ax.plot(q, cq, 'r.-')

ax = axes[1][0]
ax.set_xlabel(r"$q$ $[\AA^{-1}]$")
ax.plot(q, cq, 'r.-')
ax.set_ylim( (-1,1) )

ax = axes[0][1]
ax.set_xlabel(r"$r$ $[\AA]$")
ax.set_ylabel(r"$c(r)$")
ax.plot(r, cr, 'r.-')

ax = axes[1][1]
ax.set_xlabel(r"$r$ $[\AA]$")
ax.plot(r, cr, 'r.-')
ax.set_ylim( (-3,3) )

for axrow in axes:
    for ax in axrow:
        ax.tick_params(which='both', direction='in')
        ax.grid(axis='both', which='major')
plt.tight_layout()
fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()
