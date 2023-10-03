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
                    description = 'Computes g(r) from the Fourier transform of S(q), and plots it.',
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

outpng="rdf_from_sq.png"
outpdf="rdf_from_sq.pdf"

#-------------------------------------#
args = parser.parse_args()

Xa = np.loadtxt(args.sq, unpack=False)
q = Xa[:,0]
sq = Xa[:,1]

dq = q[1]-q[0]
qmax = q[-1]
dr = 2*np.pi/(2*qmax) # minimum frequency

r = dr * np.arange(1,1+args.Nr)
gr = np.empty(len(r))
f = (sq-1)/ args.rho / (2*np.pi)**3 # function to be integrated

for i in range(len(r)):
    integrand = 4*np.pi * q * f * np.sin(r[i]*q)/r[i]
    integral = simpson(integrand, q)
    gr[i] = 1.0 + integral
np.savetxt("rdf_from_sq.ave", np.array([r,gr]).T, header="# r | g(r) from fourier transform")

fig, axes = plt.subplots(2,1, dpi=300)

ax = axes[0]
ax.set_xlabel(r"$q$ $[\AA^{-1}]$")
ax.set_ylabel(r"$S(q)$")
ax.plot(q, sq, 'r')

ax = axes[1]
ax.set_xlabel(r"$r$ $[\AA]$")
ax.set_ylabel(r"$g(r)$")
ax.plot(r, gr, 'r')

for ax in axes:
	ax.tick_params(which='both', direction='in')
	ax.grid(axis='both', which='major')
plt.tight_layout()
fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()
