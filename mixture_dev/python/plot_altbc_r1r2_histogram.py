#!/bin/python3

import sys
import argparse
import subprocess
from typing import List
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'
cmap = colormaps['viridis']

def linMap(x,a,b,c,d):
   y=(x-a)/(b-a)*(d-c)+c
   return y

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Histogram of the coordination number for all atoms, for all timesteps, divided according to each pair of types.',
                    epilog='End of the summary.'
)
parser.add_argument('--intraj',  type=argparse.FileType('r'),
                     default="altbc_r1r2.traj", required=False,
                     help="Input file with the list of r_1/r_2 (<=1) for each line. [default: %(default)s]"
)
parser.add_argument('--fskip0', type=float,
                     default=0.0, required=False,
                     help="Fraction of data to be skipped from the start. [default: %(default)s]"
)
parser.add_argument('--fskip1', type=float,
                     default=0.0, required=False,
                     help="Fraction of data to be skipped from the end. [default: %(default)s]"
)

#-----------------------------------------------------------#

args = parser.parse_args()

outname="altbc_r1r2"

lines = args.intraj.readlines()
n = len(lines)
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
lines = lines[n0:n-n1]

#--------------------------------------------------------------------------------------------------------#

fig, axes = plt.subplots(1,2,dpi=300, figsize=(8,4))
ax = axes[0]
ax.set_xlabel(r"$r_-/r_+$ of ALTBC")
ax.set_ylabel(r"Density")
ax.tick_params(which='both', direction='in')

avg=[]
avg_=[]
count=0
bin_edges=np.linspace(0.4, 1.0, 51)
bin_centers=(bin_edges[:-1]+bin_edges[1:])/2
avg_hist = np.zeros(len(bin_centers))
for line in lines:
    if line.split()[0][0]=='#': continue #comment
    y = list(map(float, line.strip('\n').split() ))

    col = cmap( count/len(lines) )
    hist,_bin_edges_ = np.histogram(y, bins=bin_edges, density=True)
    ax.plot(bin_centers,hist, color=col, alpha=0.1)
    avg_hist += hist

    avg.append(np.array(y).mean())
    avg_.append(np.array(y).std()/np.sqrt(len(y)-1))
    count+=1
avg_hist /= count
ax.plot(bin_centers,avg_hist, color='black', alpha=0.5)

ax = axes[1]
ax.set_xlabel("Frame idx.")
ax.set_ylabel(r"$\langle r_-/r_+ \rangle_N$")
ax.errorbar(np.arange(count),avg,avg_,fmt='o', color='black')
ax.tick_params(which='both', direction='in')

plt.tight_layout()

fig.savefig(outname+".png")
print(" Figure saved on %s.png "%(outname))
#plt.show()
#subprocess.call(f"xdg-open {outpng}", shell=True)
