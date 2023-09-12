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
                    description = 'Computes S(q) from the Fourier transform of g(r), and plots it.',
                    epilog='End of the summary.'
)
parser.add_argument('--rdf',  type=argparse.FileType('r'),
                     default="rdf.ave", required=False,
                     help="Input file with average RDF (columns:r,RDF,...). [default: %(default)s]"
)
parser.add_argument('--rho',  type=float, required=True,
                     help="Number density."
)
parser.add_argument('--Nq',  type=int,
                     default=50, required=False,
                     help="Number of q points. [default: %(default)s]"
)

outpng="sq_from_rdf.png"
outpdf="sq_from_rdf.pdf"

#-------------------------------------#
args = parser.parse_args()

Xa = np.loadtxt(args.rdf, unpack=False)
r = Xa[:,0]
gr = Xa[:,1]

dr = r[1]-r[0]
rmax = r[-1]
dq = 2*np.pi/(2*rmax) # minimum frequency

q = dq * np.arange(1,1+args.Nq)
sq = np.empty(len(q))
f = (gr-1)*args.rho # function to be integrated

for i in range(len(q)):
    integrand = 4*np.pi * r * f * np.sin(q[i]*r)/q[i]
    integral = simpson(integrand, r)
#    integral = np.sum( integrand*dr )
    sq[i] = 1.0 + integral
np.savetxt("sq_from_rdf.ave", np.array([q,sq]).T, header="# q | S(q) from fourier transform")

fig, axes = plt.subplots(2,1, dpi=300)

ax = axes[0]
ax.set_xlabel(r"$r$ $[\AA]$")
ax.set_ylabel(r"$g(r)$")
ax.plot(r, gr, 'r')

ax = axes[1]
ax.set_xlabel(r"$q$ $[\AA^{-1}]$")
ax.set_ylabel(r"$S(q)$")
ax.plot(q, sq, 'r')

for ax in axes:
	ax.tick_params(which='both', direction='in')
	ax.grid(axis='both', which='major')
plt.tight_layout()
fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()
