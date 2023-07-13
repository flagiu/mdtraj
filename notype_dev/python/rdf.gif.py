#!/bin/python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Plots the "trajectory" of radial distribution functions g(r;t0) and their average g(r).',
                    epilog='End of the summary.'
)
parser.add_argument('--intraj', type=argparse.FileType('r'),
                     default="rdf.xxx", required=False,
                     help="Input file with trajectories (columns:t,x1,x2,...). [default: %(default)s]"
)
outpng="rdf.hist.png"
outpdf="rdf.hist.pdf"
x_tolerance=1e-5
gr_threshold=10

#-------------------------------------#
args = parser.parse_args()

Xt = np.loadtxt(args.intraj, unpack=False)
ntraj = Xt.shape[1]-1

def lin(x, x1,x2, y1,y2):
	return y1 + (x-x1) * (y2-y1)/(x2-x1)

fig, axes = plt.subplots(1,2,figsize=(10,6),dpi=300)
for ax in axes:
	ax.set_xlabel(r"$r$ / $\AA$")
	ax.set_ylabel(r"$g(r)$")
	ax.tick_params(which='both', direction='in')
	ax.grid(axis='both', which='major')
n0=0
n1=0
for i in range(ntraj):
	if (Xt[:,i+1]>gr_threshold).any():
		n1+=1
	else:
		n0+=1
axes[0].set_title(r"$g(r)\leq %d$ in %d frames"%(gr_threshold,n0))
axes[1].set_title(r"$g(r)$ exceeded $%d$ in %d frames"%(gr_threshold,n1))

i0=0
i1=0
for i in range(ntraj):
	if (Xt[:,i+1]>gr_threshold).any():
		iax=1
		b2r = lin( i1, 0,n1-1, 0.0,1.0 )
		i1+=1
	else:
		iax=0
		b2r = lin( i0, 0,n0-1, 0.0,1.0 )
		i0+=1
	r = b2r
	g = 0.0
	b = 1.0-b2r
	al = 0.1
	col = (r,g,b,al)
	axes[iax].plot( Xt[:,0], Xt[:,i+1], color=col)

plt.tight_layout()
fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()

