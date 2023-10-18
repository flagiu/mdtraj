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
class_labels_ordered = ['g','f','d','c','a','b','e']
class_explanations_ordered =  [
    'octahedral',
    '5-fold defective octahedral',
    '4-fold planar defective octahedral',
    '4-fold defective octahedral',
    '3-fold planar defective octahedral',
    '3-fold pyramidal defective octahedral',
    'tetrahedral',
]
class_expl_ordered =  [
    'oct.',
    '5x def. oct.',
    '4x plan. def. oct.',
    '4x def. oct.',
    '3x plan. def. oct.',
    '3x pyr. def. oct.',
    'tetr.',
]
class_q_values = [ 0., 1./3., 0.5, 5./8., 3./4., 7./8., 1. ]
class_cn_values = [6,5,4,4,3,3,4]

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
parser.add_argument('--cn', type=int,
                     default=-1, required=False,
                     help="Plot only data of particles with given coordination number (between 3 and 6); plot all if cn<0. [default: %(default)s]"
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
selection = (coordnum>=cnMIN) & (coordnum<=cnMAX)
if args.cn>=0:
    selection = selection & (coordnum==args.cn)
x = x[selection]
if len(x)==0:
    print("[ ERROR: no particle with coordination number in [%d,%d] ]\n"%(cnMIN,cnMAX))
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

fig, ax = plt.subplots(figsize=(10,6))
ax.set_xlabel("Order parameter q")
ax.set_ylabel("Probability")
title=r"$r_{cut}=%.2f$ $\AA$"%rcut1
if args.i>=0:
    title += r", particle $i=%d$"%args.i
ax.set_title(title)

cn2col = {}
for i,cn in enumerate(coordnum_u):
    cn2col[cn] = cmap(0.5) # i/(len(coordnum_u)-1) )
    counts,edges = np.histogram(data[selection], bins=args.nbins, density=None)
    ax.stairs(counts/counts.sum(), edges,
        label=r"$N_c=%d$""\n(%d events)"%(cn, counts.sum()),
        color=cn2col[cn], alpha=0.7, linewidth=2, zorder=9999 # plot on top of all!
    )
    #bin_centers = 0.5*(x[1:]+x[:-1])
    #density = gaussian_kde(n)
    # normalizzazione sbagliata!
    #ax.plot(bin_centers, density(bin_centers), color=cn2col[cn], alpha=0.7, linewidth=1, zorder=9999 )
Xmin = ax.get_xlim()[0]
Xmax = ax.get_xlim()[1]
for cn in coordnum_u:
    for i in range(len(class_q_values)):
        if class_cn_values[i]!=cn:
            continue
        xlo = min(0.0, Xmin) if i==0                   else 0.5*(class_q_values[i-1]+class_q_values[i])
        xhi = max(1.0, Xmax) if i==len(class_q_values)-1 else 0.5*(class_q_values[i]+class_q_values[i+1])
        # frac=(i+1)/(len(class_q_values)+1)
        # ax.axvspan(xlo,xhi, color=(frac,frac,frac,0.3))
        ax.axvline(class_q_values[i], linestyle='--', color=cn2col[class_cn_values[i]], alpha=0.7)
        ax.text(
            class_q_values[i]+0.01, ax.get_ylim()[-1]-0.01, class_expl_ordered[i],
            horizontalalignment='right', verticalalignment='top', rotation=90, rotation_mode='anchor',
            color=cn2col[class_cn_values[i]], alpha=0.7
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
