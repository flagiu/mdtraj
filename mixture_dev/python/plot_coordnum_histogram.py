#!/bin/python3

import sys
import argparse
import subprocess
from typing import List
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

def int2types(t, nTypes):
    for x in range(nTypes):
      low = int(x * nTypes - x*(x-1)/2)
      high = int( (x+1) * nTypes - x*(x+1)/2 - 1 )
      if t>=low and t<=high:
        #print(t,x,low,high,x, x + (int(t)-low))
        return x, x + (int(t)-low)
    return None,None

def linMap(x,a,b,c,d):
   y=(x-a)/(b-a)*(d-c)+c
   return y

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Histogram of the coordination number for all atoms, for all timesteps, divided according to each pair of types.',
                    epilog='End of the summary.'
)
parser.add_argument('--indat',  type=argparse.FileType('r'),
                     default="coordnum.dat", required=False,
                     help="Input file with the following columns: Timestep | Particle idx | coordination number for each pair 00, 01, 02, ...  [default: %(default)s]"
)
parser.add_argument('--inlabels',  type=argparse.FileType('r'),
                     default="labels.dat", required=False,
                     help="Input file with one row per type, two columns: atom label, occurrence. [default: %(default)s]"
)
parser.add_argument('--ignore',  type=int, nargs="*",
                     default=[], required=False,
                     help="Ignore the following atom types in the plot. INPUT: one or more integers (0,1,...,ntypes-1). [default: %(default)s]"
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

outpng="coordnum_hist.png"
outpdf="coordnum_hist.pdf"

header=args.indat.readline()
rcut1 = float( header.split("# cutoffs =")[1].split(',')[0] )

X = np.loadtxt(args.indat, dtype=int)
# Apply fskip
n = len(X)
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
X = X[n0:n-n1]

timesteps = X[:,0]
npairs = int( X.shape[1]-2 )
ntypes = int(np.floor( (2*npairs)**0.5 ))

#--------------------------------------------------------------------------------------------------------#
labels=[]
Nt=[]
lines = args.inlabels.readlines()
assert len(lines)==ntypes
for line in lines:
    lab, nt = line.strip('\n').split()
    labels.append(lab)
    Nt.append(int(nt))
Nt=np.array(Nt)
print(" plot_coordnum_histogram.py: Atom types:",labels,". Occurrence:",Nt,". Fraction:",Nt/np.sum(Nt))
ign=np.array(args.ignore)
ign_labels = [ labels[ig] for ig in args.ignore ]
print(" plot_coordnum_histogram.py: List of types to be ignored:",ign_labels)
#---------------------------------------------------------------------------------------------------------#

fig, axes = plt.subplots(npairs,1, figsize=(6,2*npairs), dpi=300, sharex=True)
axes[-1].set_xlabel(r"Coordination Number")
for i in range(npairs):
    ax = axes[i]
    data = X[:,2+i]
    # Bins
    d = np.diff(np.unique(data)).min()
    left_of_first_bin = data.min() - float(d)/2
    right_of_last_bin = data.max() + float(d)/2
    bins = np.arange(left_of_first_bin, right_of_last_bin + d, d)
    # Label
    ti,tj = int2types(i, ntypes)
    if len(labels)>0:
        lab = "%s-%s"%( labels[ti],labels[tj] )
    else:
        lab = "%d-%d"%( ti,tj )
    if len(ign)>0 and (ti==ign).any() or (tj==ign).any():
        continue # ignore this one

    # Color
    red = linMap(ti, 0,ntypes-1, 1,0)
    blue = linMap(tj, 0,ntypes-1, 0,1)
    green = 0.2

    ax.hist(data, bins=bins, label=lab, color=(red,green,blue,0.8), ec = "black", rwidth=0.8, density=True)
    ax.set_ylabel("Density")
    ax.set_yscale("log")
    ax.legend()
    ax.tick_params(axis='y',which='both', direction='in')
    ax.grid(axis='both', which='major')
    if i==0:
       ax.set_title(r"$r_{cut}=%.2f$ $\AA$"%rcut1)
plt.tight_layout()

fig.savefig(outpng)
fig.savefig(outpdf)
print(" plot_coordnum_histogram.py: Figure saved on %s , %s\n"%(outpng, outpdf))
#plt.show()
#subprocess.call(f"xdg-open {outpng}", shell=True)
