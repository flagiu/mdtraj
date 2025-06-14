#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

type_colors=["darkgreen","purple","red","blue","black"]

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
parser.add_argument('--outname',  type=str, required=False, default="nnd",
                     help="Prefix for output name. [default: %(default)s]"
)
parser.add_argument('--yshift',  type=float,
                     default=5.0, required=False,
                     help="Shift vertically the plots for each central atom type by this amount. [default: %(default)s]"
)
parser.add_argument('--xlim',  type=float, nargs=2,
                     default=None, required=False,
                     help="Set limit for x axis. [default: %(default)s]"
)
parser.add_argument('--ylim',  type=float, nargs=2,
                     default=None, required=False,
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
outdat=args.outname+".ave"
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

f=open(outdat,"w")
f.write("#Type of central atom; mean(r) and std(r) for each shell: r1,s1, r2,s2, ... \n")
for nt in range(ntypes): # for each type of central atom
    y_shift=nt*args.yshift
    f.write("%s"%labels[nt])
    for nn in range(M): # for each closest neighbour
        selection=(central_atom_type==nt)
        dist_forAnyNeighType = x[selection,3+2*nn+1]
        _,edges_forAnyNeighType = np.histogram(dist_forAnyNeighType, bins=args.bins)
        centers_forAnyNeighType = (edges_forAnyNeighType[1:]+edges_forAnyNeighType[:-1])/2
        avg=dist_forAnyNeighType.mean()
        std=dist_forAnyNeighType.std()
        f.write(" %.3f %.3f"%(avg,std))
        baseline_values=y_shift+np.zeros(len(centers_forAnyNeighType))
        ax.axhline(y_shift,lw=1,color='k',ls='-',zorder=-1)
        for nt_other in range(ntypes):
            neigh_type = x[:,3+2*nn]
            selection2 = (selection)&(neigh_type==nt_other)
            dist = x[selection2,3+2*nn+1]
            values,_ = np.histogram(dist, bins=edges_forAnyNeighType, density=True)
            ax.fill_between(centers_forAnyNeighType, baseline_values, baseline_values+values, color=type_colors[nt_other])
            #ax.stairs(values,edges_forAnyNeighType, baseline=baseline_values, fill=False, alpha=0.9)
            if nt_other==ntypes-1:
                ax.plot(centers_forAnyNeighType, baseline_values+values, color='k',lw=1)
            baseline_values+=values
    f.write('\n')
f.close()
print("Saved average nnd into %s"%(outdat))

#ax.set_yscale("log")
ax.tick_params(which='both', direction='in')
if args.xlim is not None: ax.set_xlim(args.xlim)
if args.ylim is not None: ax.set_ylim(args.ylim)
if ntypes>1:
    for nt in range(ntypes):
        ax.text(ax.get_xlim()[1], nt*args.yshift, labels[nt], ha='left', va='center', color=type_colors[nt])
#ax.grid(axis='both', which='major')
fig.tight_layout()

fig.savefig(args.outname+".eps")
fig.savefig(args.outname+".png")
print("Figure saved on %s.png , %s.eps\n"%(args.outname,args.outname))
#plt.show()
