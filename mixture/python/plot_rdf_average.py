#!/bin/python3

import sys
import argparse
import subprocess
from typing import List
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'
plt.rcParams['figure.dpi'] = 300

def int2types(t, nTypes):
    for x in range(nTypes):
      low = int(x * nTypes - x*(x-1)/2)
      high = int( (x+1) * nTypes - x*(x+1)/2 - 1 )
      if t>=low and t<=high:
        #print(t,x,low,high,x, x + (int(t)-low))
        return x, x + (int(t)-low)
    return None,None

def linMap(x,a,b,c,d):
   if a==b: y=c
   else:    y=(x-a)/(b-a)*(d-c)+c
   return y
linestyles = ['solid','dashed','dashdot','dotted']

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Plots the RDF for each pair of types.',
                    epilog='End of the summary.'
)
parser.add_argument('--inavg',  type=argparse.FileType('r'),
                     default="rdf.ave", required=False,
                     help="Input file with average RDF (r | g_00 g_01 ... | errors). [default: %(default)s]"
)
parser.add_argument('--inlabels',  type=argparse.FileType('r'),
                     default="labels.dat", required=False,
                     help="Input file with one row per type, two columns: atom label, occurrence. [default: %(default)s]"
)
parser.add_argument('--ignore',  type=int, nargs="*",
                     default=[], required=False,
                     help="Ignore the following atom types in the plot. INPUT: one or more integers (0,1,...,ntypes-1). [default: %(default)s]"
)
parser.add_argument('--yshift',  type=float,
                     default=0.0, required=False,
                     help="Shift vertically the different g(r)'s by this amount. [default: %(default)s]"
)
parser.add_argument('--outname',  type=str,
                     default="rdf", required=False,
                     help="Prefix for the output file. [default: %(default)s]"
)

args = parser.parse_args()

print("Plotting g(r) average ...")

X = np.loadtxt(args.inavg.name)
r = X[:,0]
npairs = int( (X.shape[1]-1)/2 )
ntypes = int(np.floor( (2*npairs)**0.5 ))

labels=[]
Nt=[]
lines = args.inlabels.readlines()
assert len(lines)==ntypes
for line in lines:
    lab, nt = line.strip('\n').split()
    labels.append(lab)
    Nt.append(int(nt))
Nt=np.array(Nt)
xt=Nt/np.sum(Nt) #Fraction

print("  Atom types:",*labels)
print("  Occurrence:",*Nt)
print("  Fraction:",*xt)
ign=np.array(args.ignore)
ign_labels = [ labels[ig] for ig in args.ignore ]
print("  List of types to be ignored:",ign_labels)

fig, ax = plt.subplots()
ax.set_xlabel(r"$r$ [$\AA$]")
ax.set_ylabel(r"$g(r)$")
if args.yshift: ax.set_title(r"Shifted by %.1f"%args.yshift, fontsize=10)
g_tot = np.zeros_like(r)
for i in reversed(range(npairs)): # reverse in order for the legend to match the plots in vertical order
    g = X[:,1+i]
    g_ = X[:,1+i+npairs]
    ti,tj = int2types(i, ntypes)
    g_tot += (xt[ti]*xt[tj]*g if ti==tj else 2*xt[ti]*xt[tj]*g) # total g(r) of all types (also ignored)
    if len(labels)>0:
        lab = "%s-%s"%( labels[ti],labels[tj] )
    else:
        lab = "%d-%d"%( ti,tj )
    if len(ign)>0 and (ti==ign).any() or (tj==ign).any():
        continue # ignore this g(r) (but keep it for g_tot(r))

    ls = linestyles[ti]
    red = linMap(tj, 0,ntypes-1, 0,1)
    blue = linMap(tj, 0,ntypes-1, 1,0)
    green = 0.0 #linMap(ti, 0,npairs-1, 0.2,1)

    #ax.errorbar(r,g,g_, label=lab, color=(red,green,blue,0.7))
    ax.plot(r,g + (i+1)*args.yshift, label=lab, color=(red,green,blue,0.7), linestyle=ls)
ax.plot(r,g_tot, label="total", color=(0,0,0,0.7))
ax.legend()
ax.tick_params(which='both', direction='in')
ax.grid(axis='both', which='major')
plt.tight_layout()

fig.savefig(args.outname+".eps")
fig.savefig(args.outname+".png")
print("Figure saved on %s.png , %s.eps\n"%(args.outname,args.outname))
#plt.show()
#subprocess.call(f"xdg-open {outpng}", shell=True)
