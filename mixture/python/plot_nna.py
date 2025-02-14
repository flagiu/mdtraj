#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

outpng="nna.png"
outpdf="nna.pdf"
outdat="nna.ave"
type_colors=["darkgreen","purple","red","blue","black"]
typePair_colors={
    (0,0):type_colors[0],
    (0,1):"cyan",
    (1,0):"cyan",
    (1,1):type_colors[1],
}

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Histogram of Nearest Neighbout Angles of all particles at all times.',
                    epilog='End of the summary.\n'
)
parser.add_argument('--indat',  type=argparse.FileType('r'), required=False, default="nnd.dat",
                     help="Input file: TimeStep | Particle0Idx | Particle0Type | (Neighbour1Type, Neighbour2Type, cos(102)) for the first M neighs, sorted by r1<=r2. [default: %(default)s]"
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

X = np.loadtxt(args.indat)
n = len(X)
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
X = X[n0:n-n1]

Mtriplets=int((X.shape[1]-3)/3) # number of triplets
M=int((np.sqrt(8*Mtriplets+1)+1)/2) # number of pairs

fig, ax = plt.subplots(dpi=200)
ax.set_xlabel(r"$\theta$ $[°]$")
ax.set_ylabel("Probability density")
y_shift=0.0
delta_yshift_triplets=args.yshift/Mtriplets
angle_bins=np.arange(0,180.1,0.1)
angle_centers=(angle_bins[1:]+angle_bins[:-1])/2
central_atom_type=X[:,2]

f=open(outdat,"w")
f.write("#Type of central atom; mean(angle°) and std(angle°) for each shell: m1,s1, m2,s2, ... \n")
for nt in range(ntypes): # for each type of central atom
    selection=(central_atom_type==nt)
    yshift=nt*args.yshift
    f.write("%s"%labels[nt])
    mm=0
    nn1nn2_ordered=[(0,1),(0,2),(0,3),(0,4),(0,5),(1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5)]
    for (nn1,nn2) in [(0,1),(0,2),(1,2), (0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5), (3,4),(3,5),(4,5)]:
        mm=[im for im in range(Mtriplets) if nn1nn2_ordered[im]==(nn1,nn2)][0]
        # write average angle for each ordered triplet, independently of atoms 1,2 type
        cos_forAnyNeighType = X[selection,3+3*mm+2]
        angle_forAnyNeighType = 180/np.pi*np.arccos(cos_forAnyNeighType)
        avg=angle_forAnyNeighType.mean()
        std=angle_forAnyNeighType.std()
        f.write(" %.3f %.3f"%(avg,std))
        #values,_ = np.histogram(angle_forAnyNeighType, bins=angle_bins, density=True)
        baseline_values=yshift+np.zeros(len(angle_centers))
        ax.axhline(baseline_values[0],color='k',lw=0.5)
        neigh1_type = X[:,3+3*mm+0]
        neigh2_type = X[:,3+3*mm+1]
        for nt_other1 in range(ntypes):
            for nt_other2 in range(nt_other1,ntypes):
                selection2 = (selection)&(neigh1_type==nt_other1)&(neigh2_type==nt_other2)
                cos = X[selection2,3+3*mm+2]
                angle = 180/np.pi*np.arccos(cos)
                values,_ = np.histogram(angle, bins=angle_bins, density=True)
                ax.plot(angle_centers,baseline_values+values,'k-')
                #ax.fill_between(angle_centers, baseline_values, baseline_values+values, color=typePair_colors[(nt_other1,nt_other2)])
                if nt_other1==ntypes-1 and nt_other2==ntypes-1:
                    ax.plot(angle_centers, baseline_values+values, color='k',lw=1)
                baseline_values+=values
        yshift+=delta_yshift_triplets

    f.write('\n')
f.close()
print("Saved average nna into %s"%(outdat))

ax.set_xlim(60,180)
ax.set_xticks([60,90,120,150,180], minor=False)
ax.set_xticks([75,105,135,165], minor=True)
ax.grid(axis='x', which='major')

#ax.set_yscale("log")
ax.tick_params(which='both', direction='in')
if args.xlim is not None: ax.set_xlim(args.xlim)
if args.ylim is not None: ax.set_ylim(args.ylim)
if ntypes>1:
    for nt in range(ntypes):
        ax.text(ax.get_xlim()[1], nt*args.yshift, labels[nt], ha='left', va='center', color=type_colors[nt])
#ax.grid(axis='both', which='major')
fig.tight_layout()

fig.savefig(outpng)
fig.savefig(outpdf)
print(" Figure saved on %s , %s\n"%(outpng, outpdf))
#plt.show()
