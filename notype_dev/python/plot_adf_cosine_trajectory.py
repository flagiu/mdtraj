#!/bin/python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Plots the "trajectory" of angular distribution functions ADF(cosine(angle);t0) and their average ADF(cosine(angle)).',
                    epilog='End of the summary.'
)
parser.add_argument('--intraj', type=argparse.FileType('r'),
                     default="adf.xxx", required=False,
                     help="Input file with trajectories (columns:t,x1,x2,...). [default: %(default)s]"
)
parser.add_argument('--inavg',  type=argparse.FileType('r'),
                     default="adf.ave", required=False,
                     help="Input file with average trajectory (columns:t,RDF,RDF error). [default: %(default)s]"
)

outpng="adf_cosine.png"
outpdf="adf_cosine.pdf"
x_tolerance=5e-5

#-------------------------------------#
args = parser.parse_args()

Xt = np.loadtxt(args.intraj, unpack=False)
ntraj = Xt.shape[1]-1

Xa = np.loadtxt(args.inavg, unpack=False)
x = Xa[:,0]
y = Xa[:,1]

assert Xt.shape[0]==Xa.shape[0] # space binning must have same length
x_isnot_equal = abs( (Xt[:,0]-Xa[:,0])/Xa[:,0] )>=x_tolerance
try:
	assert not x_isnot_equal.any() # space binning must be the same
except AssertionError:
	print("[ ERROR: x values do not match! ]")
	print(" File",args.intraj)
	print(Xt[x_isnot_equal, 0])

	print(" File",args.inavg)
	print(Xa[x_isnot_equal, 0])
	sys.exit(1)

assert Xa.shape[1]>=2 # must have at least x,y

fig, ax = plt.subplots(dpi=300)
ax.set_xlabel(r"$\cos(\theta)$")
ax.set_ylabel(r"Density")
#ax.xaxis.set_ticks(np.arange(0.0, 180.0, 15.0))
ax.tick_params(which='both', direction='in')
for i in range(ntraj):
	ax.plot( x, Xt[:,i+1], 'k', alpha=0.01)
ax.plot(x, y, 'r.-')
ax.set_xlim((-1, 1))
#ax.set_xticks(np.arange(0, 180+15, step=15))
ax.grid(axis='both', which='major')
plt.tight_layout()
fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s , %s\n"%(outpng, outpdf))
#plt.show()
