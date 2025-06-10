#!/bin/python3

import sys
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'
plt.rcParams['figure.dpi'] = 200


parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Plots the Intermediate Scattering Function S(q,t).',
                    epilog='Both as a function of t at parametric q, and vice-versa.'
)
parser.add_argument('--inavg',  type=argparse.FileType('r'),
                     default="sqt.ave", required=False,
                     help="Input file with average trajectory (1st block: wave vector bis; 2nd block: Delta timestep bins; 3rd block: 2d matrix of <S(q,t)>; 4th block; 2d matrix of its error. [default: %(default)s]"
)
parser.add_argument('--dt',  type=float, default=0.002, required=False,
                     help="Integration time step (picoseconds). [default: %(default)s]"
)
parser.add_argument('--select_q',  type=int, nargs="+", default=[], required=False,
                     help="Plot vs. t only the q's indexed by the selected integers. [default: %(default)s]"
)
parser.add_argument('--select_t',  type=int, nargs="+", default=[], required=False,
                     help="Plot vs. q only the t's indexed by the selected integers. [default: %(default)s]"
)
parser.add_argument('--yshift',  type=float, default=0.0, required=False,
                     help="Shift vertically the vs-time-plot. [default: %(default)s]"
)
parser.add_argument('--normalize',  type=int, default=1, required=False,
                     help="Normalize the vs-time-plot to S(q,0)? 1: yes ; 0: no. [default: %(default)s]"
)
parser.add_argument('--ylog',  type=int, default=1, required=False,
                     help="Use log y-axis for the vs-q-plot? 1: yes ; 0: no. [default: %(default)s]"
)
parser.add_argument('--xlog',  type=int, default=1, required=False,
                     help="Use log x-axis for the vs-t-plot? 1: yes ; 0: no. [default: %(default)s]"
)
parser.add_argument('--fmt',  type=str, default="-", required=False,
                     help="Use the given format for errorbar plotting. [default: %(default)s]"
)
parser.add_argument('--outname',  type=str,
                     default="sqt", required=False,
                     help="Prefix for the output file. [default: %(default)s]"
)

#-------------------------------------#
args = parser.parse_args()
assert args.normalize==0 or args.normalize==1
assert args.ylog==0 or args.ylog==1

print("Plotting monospecies S(q,t) ...")
lines = args.inavg.readlines()
# Parse wavevectors from 1st block:
i = 0
block=0
q = []
while block==0:
    line = lines[i]
    words = line.strip('\n').split(' ')
    if line[0]=='#':
        words=None
    else:
        assert len(words)==1
        if words[0]=='':
            block=1 # separation empty line
        else:
            q.append(float(words[0]))
    i+=1
q = np.array(q)
Nq=len(q)

# Parse delta timesteps from 2nd block:
t = []
while block==1:
    line = lines[i]
    words = line.strip('\n').split(' ')
    if line[0]=='#':
        words=None
    else:
        assert len(words)==1
        if words[0]=='':
            block=2 # separation empty line
        else:
            t.append(float(words[0]))
    i+=1
t = args.dt * np.array(t)
Nt=len(t)

# Load 3rd block
sqt = np.genfromtxt(args.inavg.name, skip_header=1+Nq+1+Nt+1, skip_footer=Nq)
# Load 4th block
sqt_ = np.genfromtxt(args.inavg.name, skip_header=1+Nq+1+Nt+1+Nq+1, skip_footer=0)
if Nq==1 or Nt==1:
    assert len(sqt)==Nt or len(sqt)==Nq
    assert len(sqt_)==Nt or len(sqt_)==Nq
    sqt = sqt.reshape(Nq,Nt)
    sqt_ = sqt_.reshape(Nq,Nt)
else:
    assert sqt.shape==(Nq,Nt)
    assert sqt_.shape==(Nq,Nt)
print("  Found %d t's and %d q's."%(Nt,Nq))

# Extract S(q,0)=S(q)
sq = sqt[:,0]
sq_ = sqt_[:,0]

# Apply selection for t and q
if args.select_t!=[]:
    t_selected_idx = args.select_t
    t_selected = [t[el] for el in t_selected_idx ]
    print("  Selected t's:",*args.select_t,".")
    print("          i.e.:",*t_selected,".")
else:
    t_selected_idx = np.arange(len(t)) # all indexes
    t_selected = t

if args.select_q!=[]:
    q_selected_idx = args.select_q
    q_selected = [q[el] for el in q_selected_idx ]
    print("  Selected q's:",*args.select_q,".")
    print("          i.e.:",*q_selected,".")
else:
    q_selected_idx = np.arange(len(q)) # all indexes
    q_selected = q

def add_upperxticks_to_cbar(cbar,ticks):
    ax_upper=cbar.ax.twiny()
    ax_upper.set_xlim(cbar.ax.get_xlim())
    ax_upper.set_xticks(ticks)
    ax_upper.set_xticklabels(["%.2f"%el for el in ticks])
#-------------------------------------------------------------------------------------------------------------------#

fig, axes = plt.subplots(2,2, figsize=(10,5),gridspec_kw={'height_ratios': [1, 9]} ) # figsize is experimental
ax = axes[1][0]
ax.axhline(0, color="k", ls="--", lw=0.5, zorder=-999)
ax.set_xlabel(r"$q$ $[\AA^{-1}]$")
ax.set_ylabel(r"$S(q,t)$")
cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=t.min(), vmax=t.max())
colors = cmap(norm(t))
for ii in range(len(t_selected_idx)):
    i = t_selected_idx[ii]
    lab = r"$t=%.1f$ $ps$"%t_selected[ii]
    col = colors[i]
    y = sqt[:,i]
    y_ = sqt_[:,i]
    ax.errorbar(q,y,y_, fmt=args.fmt, color=col, label=lab, alpha=0.5)
cbar_t=fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=axes[0][0],
                    orientation="horizontal", label="t [ps]")
if args.select_t!=[]: add_upperxticks_to_cbar(cbar_t,t_selected)
ax.tick_params(which='both', direction='in')
if args.ylog:
    ax.set_yscale("log")
    ax.set_ylim(max(ax.get_ylim()[0],1e-3),ax.get_ylim()[1])


ax = axes[1][1]
ax.set_xlabel(r"$t$ [$ps$]")
if args.normalize==1:
	ax.set_ylabel(r"$S(q,t)$ / $S(q,0)$")
else:
	ax.set_ylabel(r"$S(q,t)$")
cmap = mpl.cm.magma
norm = mpl.colors.Normalize(vmin=q.min(), vmax=q.max())
colors = cmap(norm(q))
for ii in range(len(q_selected_idx)):
    i = q_selected_idx[ii]
    lab = r"$q=%.2f$ $\AA^{-1}$"%q_selected[ii]
    axes[1][0].axvline(q_selected[ii],ls='--',lw=0.5,color='k',zorder=-1)
    col = colors[i]
    if args.normalize==1:
        y = sqt[i] / sq[i]
        y_ = y * np.sqrt( (sqt_[i]/sqt[i])**2 + (sq_[i]/sq[i])**2 ) # gaussian error propagation
        #print(sq[i],sq_[i])
        #print(sqt[i],sqt_[i])
        #print(y,y_)
    else:
        y = sqt[i]
        y_ = sqt_[i]
    ax.errorbar(t,y+ii*args.yshift,y_, fmt=args.fmt, color=col, label=lab, alpha=0.5)
if args.normalize==1:
    ax.set_ylim(-0.1,1.1)
    ax.axhline(np.exp(-1), color="k", ls="--", lw=0.5, zorder=-999)
ax.tick_params(which='both', direction='in')
ax.axhline(0, color="k", ls="--", lw=0.5, zorder=-999)
#ax.grid(axis='both', which='major')
cbar_q=fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=axes[0][1],
                    orientation="horizontal", label=r"q [$\AA^{-1}$]")
if args.select_q!=[]: add_upperxticks_to_cbar(cbar_q,q_selected)
if args.xlog:
    ax.set_xscale("log")

plt.tight_layout()

fig.savefig(args.outname+".eps")
fig.savefig(args.outname+".png")
print("Figure saved on %s.png , %s.eps\n"%(args.outname,args.outname))
#plt.show()
