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
                    description = 'Plots the coordination number for each pair of types, averaged on all atoms.',
                    epilog='End of the summary.'
)
parser.add_argument('--inavg',  type=argparse.FileType('r'),
                     default="coordnum.ave", required=False,
                     help="Input file with the following columns: Timestep | <coordination number> for each pair 00, 01, 02, ... | Fluctuations for each pair.. [default: %(default)s]"
)
parser.add_argument('--inlabels',  type=argparse.FileType('r'),
                     default="labels.dat", required=False,
                     help="Input file with one row per type, two columns: atom label, occurrence. [default: %(default)s]"
)
parser.add_argument('--ignore',  type=int, nargs="*",
                     default=[], required=False,
                     help="Ignore the following atom types in the plot. INPUT: one or more integers (0,1,...,ntypes-1). [default: %(default)s]"
)

args = parser.parse_args()

outpng="coordnum_ave.png"
outpdf="coordnum_ave.pdf"

header=args.inavg.readline()
rcuts_str = header.split("# cutoffs =")[1].strip('\n').split()
rcuts = [float(rc) for rc in rcuts_str]

X = np.loadtxt(args.inavg)
timesteps = X[:,0]
npairs = int( (X.shape[1]-1)/2 )
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
print(" plot_coordnum_average.py: Atom types:",labels,". Occurrence:",Nt,". Fraction:",Nt/np.sum(Nt))
ign=np.array(args.ignore)
ign_labels = [ labels[ig] for ig in args.ignore ]
print(" plot_coordnum_average.py: List of types to be ignored:",ign_labels)
#---------------------------------------------------------------------------------------------------------#
if len(labels)==0:
   for i in range(ntypes):
      labels.append(str(i))

fig, ax = plt.subplots(dpi=300)
ax.set_xlabel(r"Timestep")
ax.set_ylabel(r"Average Coordination Number")
for i in range(npairs):
    c = X[:,1+i]
    c_ = X[:,1+i+npairs]
    ti,tj = int2types(i, ntypes)
    lab = r"%s-%s, $r_{cut}=%.2f$ $\AA$"%( labels[ti],labels[tj], rcuts[i] )

    if len(ign)>0 and (ti==ign).any() or (tj==ign).any():
        continue # ignore this one

    red = linMap(ti, 0,ntypes-1, 1,0)
    blue = linMap(tj, 0,ntypes-1, 0,1)
    green = 0.2
    ax.errorbar(timesteps,c,c_, fmt='.-', label=lab, color=(red,green,blue,0.7))
ax.legend()
ax.tick_params(which='both', direction='in')
ax.grid(axis='both', which='major')
plt.tight_layout()

fig.savefig(outpng)
fig.savefig(outpdf)
print(" plot_coordnum_average.py: Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()
#subprocess.call(f"xdg-open {outpng}", shell=True)
