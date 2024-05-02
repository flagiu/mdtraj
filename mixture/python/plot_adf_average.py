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
                    description = 'Plots the ADF for each type of central atom and for each pair of types for the edge atoms.',
                    epilog='End of the summary.'
)
parser.add_argument('--inavg',  type=argparse.FileType('r'),
                     default="adf.ave", required=False,
                     help="Input file with average ADF (cosine | 0-0-0 0-0-1 0-0-2 ... 1-0-1 1-0-2 ... 0-1-0 0-1-1 ... | errors). [default: %(default)s]"
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
                     help="Shift vertically the different ADFs by this amount. [default: %(default)s]"
)
parser.add_argument('--scale',  type=float, nargs="*",
                     default=[1.0], required=False,
                     help="Scale the different ADFs by this amount (same for all, or one for each triplet that is not ignored). [default: %(default)s]"
)
parser.add_argument('--famous_angles',  type=float, nargs="*",
                     default=[109.5], required=False,
                     help="Mark these angles (in degrees) with vertical dotted lines. [default: %(default)s]"
)

args = parser.parse_args()

print("Plotting ADF average ...")
outname="adf"

X = np.loadtxt(args.inavg.name)
x = X[:,0]
ntriplets = int( (X.shape[1]-1)/2 )
ntypes = int(np.floor( (2*ntriplets)**(1/3.) ))
npairs = int(ntypes*(ntypes+1)/2)

if len(args.scale)<1 or len(args.scale)>ntriplets:
    print("ERROR: --yscale must be followed by 1 or n floats, with n = number of triplets that are not ignored")
    print("       nargs = %d ; num_triplets = %d ; num_ignored_types = %d"%(len(args.scale),ntriplets,len(ignored)))
    sys.exit(1)

for ang in args.famous_angles:
    if ang<0. or ang>180.:
        print("ERROR: some argument in --famous_angles is not in degrees")
        sys.exit(1)

#----------- Convert from cos(angle) to angle ----------------#
jacob = np.sqrt(1-x*x)
x = np.arccos( x ) * 180/np.pi # (degrees)
for i in range(ntriplets*2):
    X[:,i+1] *= jacob
#-------------------------------------------------------------#

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
print("  Index | Triplet | Scaling factor")
print("  --------------------------------")

fig, ax = plt.subplots()
ax.set_xlabel(r"Angle [degrees]")
ax.set_ylabel(r"ADF")
if args.yshift: ax.set_title(r"Shifted by %.1f"%args.yshift, fontsize=10)
y_tot = np.zeros_like(x)
num_triplets_notignored=int(0)
for ti in range(ntypes):
    for tp in range(npairs):
        idx = ti*npairs+tp
        y = X[:,1+idx]
        y_ = X[:,1+idx+ntriplets]
        tj,tk = int2types(tp, ntypes)
        y_tot += 0.0 #(xt[ti]*xt[tj]*xt[tk]*y if tj==tk else 2*xt[ti]*xt[tj]*x[tk]*y)
        if len(labels)>0:
            lab = "%s-%s-%s"%( labels[tj],labels[ti],labels[tk] )
        else:
            lab = "%d-%d-%d"%( tj,ti,tk )

        if len(ign)>0 and (ti==ign).any() or (tj==ign).any() or (tk==ign).any():
            continue # ignore this y (but keep it for y_tot)

        try: scale=args.scale[num_triplets_notignored]
        except:
            print("*-*-*-* you probably gave too many arguments to --yscale! *-*-*-*")
            raise

        print("  %6d %s %.2f"%(idx, lab, scale))
        if scale!=1.0:
            lab+= " (%2gx)"%scale


        ls = linestyles[ti]
        red = linMap(tj, 0,ntypes-1, 0,1)
        blue = linMap(tk, 0,ntypes-1, 1,0)
        green = 0.0 #linMap(ti, 0,npairs-1, 0.2,1)

        #ax.errorbar(r,g,g_, label=lab, color=(red,green,blue,0.7))
        y_adjusted = y*scale + (idx+1)*args.yshift
        ax.plot(x,y_adjusted, label=lab, color=(red,green,blue,0.7), linestyle=ls)
        ax.text(180+2, (idx+1)*args.yshift, lab, ha='left', va='center')
        num_triplets_notignored+=1

print("  --------------------------------")
"""
if ntypes>1:
    ax.plot(x,y_tot/ntriplets, label="total (not implemented)", color=(0,0,0,0.7))
    ax.text(180+2, 0.0, "total", horizontalalignment='left', verticalalignment='center')
"""
ax.set_xlim(0,180)
ax.set_xticks([0,30,60,90,120,150,180], minor=False)
ax.set_xticks([15,45,75,105,135,175], minor=True)
for ang in args.famous_angles:
    ax.axvline(ang, linestyle=":", color='gray', zorder=-999)
    ax.text(ang+0.5, 0.99*ax.get_ylim()[1], "%.1f°"%ang,
            horizontalalignment='left', verticalalignment='top', fontsize=8)
#ax.legend()
ax.tick_params(which='both', direction='in')
ax.grid(axis='both', which='major', alpha=0.3)
plt.tight_layout()

fig.savefig(outname+".png")
fig.savefig(outname+".png")
print("Figure saved on %s.png , %s.pdf\n"%(outname,outname))
#plt.show()
#subprocess.call(f"xdg-open {outpng}", shell=True)
