#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

outpng="nnd.png"
outpdf="nnd.pdf"

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Histogram of Nearest Neighbout distances of all particles at all times.',
                    epilog='End of the summary.\n'
)
parser.add_argument('--indat',  type=argparse.FileType('r'), required=False, default="nnd.dat",
                     help="Input file: TimeStep | ParticleIdx | ParticleType | (NeighbourType, DistanceSquared) for the first M neighs, sorted by distance. [default: %(default)s]"
)
parser.add_argument('--inlabels',  type=argparse.FileType('r'), required=False, default="labels.dat",
        help="Two column file: Labels | Number of particles of that type. [default: %(default)s]"
)
parser.add_argument('--bins',  type=int, required=False, default=100,
        help="Number of bins for each neighbour histogram. [default: %(default)s]"
)
parser.add_argument('--yshift',  type=float,
                     default=5.0, required=False,
                     help="Shift vertically the plots for each central atom type by this amount. [default: %(default)s]"
)
parser.add_argument('--xlim',  type=float, nargs=2,
                     default=(2.6,4.9), required=False,
                     help="Set limit for x axis. [default: %(default)s]"
)
parser.add_argument('--ylim',  type=float, nargs=2,
                     default=(0,6), required=False,
                     help="Set limit for y axis. [default: %(default)s]"
)
parser.add_argument('--fskip0', type=float,
                     default=0.0, required=False,
                     help="Fraction of data to be skipped from the start. [default: %(default)s]"
)
parser.add_argument('--fskip1', type=float,
                     default=0.0, required=False,
                     help="Fraction of data to be skipped from the end. [default: %(default)s]"
)

args = parser.parse_args()

labels=[]
Nt=[]
lines = args.inlabels.readlines()
ntypes=len(lines)
for line in lines:
    lab, nt = line.strip('\n').split()
    labels.append(lab)
    Nt.append(int(nt))
Nt=np.array(Nt)
xt=Nt/np.sum(Nt) #Fraction

x = np.loadtxt(args.indat)
n = len(x)
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
x = x[n0:n-n1]

M=int((x.shape[1]-3)/2)
# take the square root of all distances
for nn in range(M): x[:,3+2*nn+1]=x[:,3+2*nn+1]**0.5

fig, ax = plt.subplots(dpi=200)
ax.set_xlabel(r"r [$\AA$]")
ax.set_ylabel("Density")
y_shift=0.0
central_atom_type=x[:,2]

print("#Central atom type, Average distance for each shell:")
for nt in range(ntypes): # for each type of central atom
    y_shift=nt*args.yshift
    print(labels[nt],end='')
    for nn in range(M): # for each closest neighbour
        selection=(central_atom_type==nt)
        dist_forAnyNeighType = x[selection,3+2*nn+1]
        _,edges_forAnyNeighType = np.histogram(dist_forAnyNeighType, bins=args.bins)
        avg=dist_forAnyNeighType.mean()
        print(" %.2f"%avg,end='')
        baseline_values=y_shift+np.zeros(len(edges_forAnyNeighType)-1)
        for nt_other in range(ntypes):
            neigh_type = x[selection,3+2*nn]
            selection2 = (selection)&(neigh_type==nt_other)
            dist = x[selection2,3+2*nn+1]
            values,_ = np.histogram(dist, bins=edges_forAnyNeighType, density=True)
            ax.stairs(values,edges_forAnyNeighType, baseline=baseline_values, fill=False, alpha=0.9)
            baseline_values+=values
    print()

#ax.set_yscale("log")
ax.tick_params(which='both', direction='in')
if args.xlim is not None: ax.set_xlim(args.xlim)
if args.ylim is not None: ax.set_ylim(args.ylim)
if ntypes>1:
    for nt in range(ntypes):
        ax.text(ax.get_xlim()[1], nt*args.yshift, "Around "+labels[nt], ha='left', va='center')
#ax.grid(axis='both', which='major')
plt.tight_layout()

plt.savefig(outpng)
plt.savefig(outpdf)
print(" Figure saved on %s , %s\n"%(outpng, outpdf))
#plt.show()
