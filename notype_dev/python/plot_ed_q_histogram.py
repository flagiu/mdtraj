#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import argparse
import numpy as np
from scipy.stats import gaussian_kde
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'
plt.rcParams['figure.dpi'] = 300
cmap = mpl.colormaps['brg']

outpng="ed_q_hist.png"
outpdf="ed_q_hist.pdf"
cnMIN=2
cnMAX=12
class LocalStructure:
    def __init__(self, label:str, explanation:str, expl:str, cn:int, q:float, q_boundaries ):
        self.label = label
        self.explanation = explanation
        self.expl = expl
        self.cn = int(cn)
        self.q = float(q)
        self.q_boundaries = q_boundaries

localStructures = (
    LocalStructure( 'a', '3-fold planar defective octahedral', '3x plan. def. oct.', 3, 3./4., (0.6,0.8) ),
    LocalStructure( 'b', '3-fold pyramidal defective octahedral', '3x pyr. def. oct.', 3, 7./8., (0.8,1.0) ),
    LocalStructure( 'c', '4-fold defective octahedral', '4x def. oct.', 4, 5./8., (.5,.8) ),
    LocalStructure( 'd', '4-fold planar defective octahedral', '4x plan. def. oct.', 4, .5, (0.4,0.5) ),
    LocalStructure( 'e', 'tetrahedral', 'tetr.', 4, 1., (0.8,1.0) ),
    LocalStructure( 'f', '5-fold defective octahedral', '5x def. oct.', 5, 1./3., (.1,.5) ),
    LocalStructure( 'g', 'octahedral', 'oct.', 6, 0., (-.1,.1) )
)

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Histogram of Eddington-Debenedettin "q" bond order parameter, of all particles at all times.',
                    epilog='End of the summary.\n'
)
parser.add_argument('--indat',  type=argparse.FileType('r'), required=False, default="ed_q.dat",
                     help="Input file: timestep | particle index | coordination number | q. [default: %(default)s]"
)
parser.add_argument('--nbins', type=int,
                     default=50, required=False,
                     help="Number of bins for the histogram. [default: %(default)s]"
)
parser.add_argument('--i', type=int,
                     default=-1, required=False,
                     help="Plot only data of particle i; i can be 0,1,...,N-1 . Plot all if i<0. [default: %(default)s]"
)
parser.add_argument('--cn', type=int, nargs=2,
                     default=[cnMIN, cnMAX], required=False,
                     help="Plot only data of particles with coordination number within the interval [min,max]. [default: %(default)s]"
)
parser.add_argument('--fskip0', type=float,
                     default=0.0, required=False,
                     help="Fraction of data to be skipped from the start. [default: %(default)s]"
)
parser.add_argument('--fskip1', type=float,
                     default=0.0, required=False,
                     help="Fraction of data to be skipped from the end. [default: %(default)s]"
)
parser.add_argument('--logScale', type=bool,
                     default=False, required=False,
                     help="Use log scale for y axis?. [default: %(default)s]"
)
parser.add_argument('--density', type=bool,
                     default=False, required=False,
                     help="Plot probability density? Else, plot probability. [default: %(default)s]"
)

args = parser.parse_args()

header=args.indat.readline()
rcut1 = float( header.split("# cutoff =")[1].split(',')[0] )

x = np.loadtxt(args.indat)
#-------- FILTERS -----------------------------#
# filter by particle
indices = x[:,1]
if args.i>=0:
    x = x[indices==args.i]
if len(x)==0:
    print("[ ERROR: no particle with index %d]\n"%args.i)
    exit(1)
# filter by coordination number
coordnum = x[:,2]
if len(args.cn)==2 and args.cn[0]>0 and args.cn[1]>=args.cn[0]:
    selection = (coordnum>=args.cn[0]) & (coordnum<=args.cn[1])
else:
    print("[ ERROR: argument --cn must be followed by 2 positive integers. ]\n")
    exit(1)
x = x[selection]
if len(x)==0:
    print("[ ERROR: no particle with coordination number in the given interval: {args.cn} ]\n")
    exit(1)
# filter by time
t = np.unique(x[:,0])
n = len(t) # number of frames
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
if n-n0-n1 <=0:
    print("[ ERROR: cannot skip all frames! nskip0=%d, nskip1=%d, but there are n=%d frames! ]\n"%(n0,n1,n))
    exit(1)
t0 = t[n0]
t1 = t[n-1-n1]
selection = (x[:,0]>=t0) & (x[:,0]<=t1)
x = x[ selection ]
#----------------------------------------------#
# filtered variables
t = np.unique(x[:,0])
indices = x[:,1]
coordnum = x[:,2]
coordnum_u = np.unique(coordnum)
data = x[:,3]
#---------------------------------------------------------#
# analysis
f = open("ed_q_structures.txt","w")
f.write("[Local structure]\n")
f.write("| Counts | % within same coord.num. |\n")
f.write("|-----------------------------------|\n")
for ls in localStructures:
    f.write(f"[{ls.expl}]\n")
    for typ in types_u:
        selection = (coordnum==ls.cn)
        total = len(data[selection])
        selection = selection & (data>ls.q_boundaries[0]) & (data<=ls.q_boundaries[1])
        count = len(data[selection])
        f.write(f"{count}\t{count/total*100:.2f}\n")
f.write("|-----------------------------------|\n")
f.close()
print("%s: Analysis printed to ed_q_structures.txt \n"%(sys.argv[0]))
#---------------------------------------------------------#
# plot
fig, ax = plt.subplots(figsize=(10,6))
ax.set_xlabel("Order parameter q")
ax.set_ylabel("Probability density") if args.density else: ax.set_ylabel("Probability")
title=r"$r_{cut}=%.2f$ $\AA$"%rcut1
if args.i>=0:
    title += r", particle $i=%d$"%args.i
ax.set_title(title)

cn2col = {}
for i,cn in enumerate(coordnum_u):
    cn2col[cn] = cmap( 0.5 if len(coordnum_u)==1 else i/(len(coordnum_u)-1) )
    selection = coordnum==cn
    counts,edges = np.histogram(data[selection], bins=args.nbins, density=args.density)
    total = len(data[selection])
    if args.density:
        probs = counts
    else:
        probs = counts/total if total>0 else counts
    ax.stairs(probs, edges,
        label=r"$N_c=%d$""\n(%d events)"%(cn, total),
        color=cn2col[cn], alpha=0.7, linewidth=2, zorder=9999 # plot on top of all!
    )
    #bin_centers = 0.5*(x[1:]+x[:-1])
    #density = gaussian_kde(n)
    # normalizzazione sbagliata!
    #ax.plot(bin_centers, density(bin_centers), color=cn2col[cn], alpha=0.7, linewidth=1, zorder=9999 )
Xmin = ax.get_xlim()[0]
Xmax = ax.get_xlim()[1]
for ls in localStructures:
    if ls.cn!=cn:
        continue
    #ax.axvspan(ls.q_boundaries[0], ls.q_boundaries[1], color=cmap(ls.q), alpha=0.1)
    ax.axvline(ls.q, linestyle='--', color=cn2col[ls.cn], alpha=0.7)
    ax.text(
        ls.q-0.01, ax.get_ylim()[-1]-0.01, ls.expl,
        horizontalalignment='right', verticalalignment='bottom', rotation=90, rotation_mode='anchor',
        color=cn2col[ls.cn], alpha=0.7
    )
ax.set_xlim( (min(0.0, Xmin),max(1.0, Xmax)) )
if args.logScale:
    ax.set_yscale("log")
ax.legend() #bbox_to_anchor=(1.01,1.01))
ax.tick_params(axis='y',which='both', direction='in')
#ax.grid(axis='both', which='major')
plt.tight_layout()

plt.savefig(outpng)
plt.savefig(outpdf)
print("%s: Figure saved on %s , %s\n"%(sys.argv[0],outpng, outpdf))
#plt.show()
