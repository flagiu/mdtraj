#!/bin/python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Computes S(q) from the Fourier transform of g(r), and plots it.',
                    epilog='End of the summary.'
)
parser.add_argument('--rdf',  type=argparse.FileType('r'),
                     default="rdf.ave", required=False,
                     help="Input file with average RDF (columns:t,RDF,RDF error). [default: %(default)s]"
)
parser.add_argument('--rho',  type=float, required=True,
                     help="Number density."
)

outpng="sq_from_rdf.png"
outpdf="sq_from_rdf.pdf"

#-------------------------------------#
args = parser.parse_args()

Xa = np.loadtxt(args.rdf, unpack=False)
r = Xa[:,0]
gr = Xa[:,1]

r_upper = r+r[0] # = dr, 2*dr, ..., L/2
dr = r_upper[0]
dq = np.pi/r_upper[-1]
q = dq * np.arange(1,len(r)/2)
sq = np.empty(len(q))

for i in range(len(q)):
	integrand = r * gr * np.sin(q[i]*r)/q[i]
	sq[i] = 1.0 + 4*np.pi * args.rho * np.sum( integrand*dr )
print(q)
print(sq)

fig, axes = plt.subplots(2,1, dpi=300)

ax = axes[0]
ax.set_xlabel(r"$r$ / $\AA$")
ax.set_ylabel(r"$g(r)$")
ax.plot(r, gr, 'r')

ax = axes[1]
ax.set_xlabel(r"$q$ / $\AA^{-1}$")
ax.set_ylabel(r"$S(q)$")
ax.plot(q, sq, 'r')

for ax in axes:
	ax.tick_params(which='both', direction='in')
	ax.grid(axis='both', which='major')
plt.tight_layout()
plt.show()
fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()
