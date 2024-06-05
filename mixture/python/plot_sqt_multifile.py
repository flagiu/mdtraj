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
my_path="/".join(sys.argv[0].split("/")[:-1])
sys.path.append(my_path)
from read_sqt import read_sqt


parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Plots the Intermediate Scattering Function S(q,t).',
                    epilog='Both as a function of t at parametric q, and vice-versa.'
)
parser.add_argument('--files',  type=argparse.FileType('r'), nargs="+",
                     default=["sqt.ave"], required=False,
                     help="Input file(s) with average trajectory (1st block: wave vector bis; 2nd block: Delta timestep bins; 3rd block: 2d matrix of <S(q,t)>; 4th block; 2d matrix of its error. [default: %(default)s]"
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

outname="sqt"

#-------------------------------------#
args = parser.parse_args()
assert args.normalize==0 or args.normalize==1
assert args.ylog==0 or args.ylog==1
assert args.xlog==0 or args.xlog==1

fig, axes = plt.subplots(2,2, figsize=(10,5),gridspec_kw={'height_ratios': [1, 9]} ) # figsize is experimental
#----------------------------------------------------------#
cmap_t = mpl.cm.viridis
cmap_q = mpl.cm.magma
min_q=9999
max_q=-9999
min_t=9999
max_t=-9999
for file in args.files:
    q,t,sqt,sqt_=read_sqt(file.name,dt=args.dt)
    min_t=min(min_t,t.min())
    max_t=max(max_t,t.max())
    min_q=min(min_q,q.min())
    max_q=max(max_q,q.max())
norm_q = mpl.colors.Normalize(vmin=min_q, vmax=max_q)
norm_t = mpl.colors.Normalize(vmin=min_t, vmax=max_t)
#----------------------------------------------------------#

for file in args.files:
    q,t,sqt,sqt_=read_sqt(file.name,dt=args.dt)

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

    ax = axes[1][0]
    ax.set_xlabel(r"$q$ $[\AA^{-1}]$")
    ax.set_ylabel(r"$S(q,t)$")
    colors = cmap_t(norm_t(t))
    for ii in range(len(t_selected_idx)):
        i = t_selected_idx[ii]
        lab = r"$t=%.1f$ $ps$"%t_selected[ii]
        col = colors[i]
        x=q
        y = sqt[:,i]
        y_ = sqt_[:,i]
        ax.errorbar(x,y,y_, fmt=args.fmt, color=col, label=lab, alpha=0.5)

    ax = axes[1][1]
    ax.set_xlabel(r"$t$ [$ps$]")
    if args.normalize==1: ax.set_ylabel(r"$S(q,t)$ / $S(q,0)$")
    else: ax.set_ylabel(r"$S(q,t)$")
    colors = cmap_q(norm_q(q))
    for ii in range(len(q_selected_idx)):
        i = q_selected_idx[ii]
        lab = r"$q=%.2f$ $\AA^{-1}$"%q_selected[ii]
        col = colors[i]
        x=t
        if args.normalize==1:
            y = sqt[i] / sq[i]
            y_ = y * np.sqrt( (sqt_[i]/sqt[i])**2 + (sq_[i]/sq[i])**2 ) # gaussian error propagation
        else:
            y = sqt[i]
            y_ = sqt_[i]
        if args.xlog==1: # don't plot t=0
            x=x[1:]
            y=y[1:]
            y_=y_[1:]
        ax.errorbar(x,y+i*args.yshift,y_, fmt=args.fmt, color=col, label=lab, alpha=0.5)

ax.set_ylim(-0.1,1.1)
ax.axhline(np.exp(-1), color="k", linestyle="--", linewidth=0.5, zorder=-999)
ax.tick_params(which='both', direction='in')
#ax.grid(axis='both', which='major')

ax=axes[1][0]
cbar_t=fig.colorbar(mpl.cm.ScalarMappable(norm=norm_t, cmap=cmap_t), cax=axes[0][0],
                    orientation="horizontal", label="t [ps]")
if args.select_t!=[]: add_upperxticks_to_cbar(cbar_t,t_selected)
ax.tick_params(which='both', direction='in')
if args.ylog==1:
    ax.set_yscale("log")
    ax.set_ylim(max(ax.get_ylim()[0],1e-3),ax.get_ylim()[1])

ax=axes[1][1]
cbar_q=fig.colorbar(mpl.cm.ScalarMappable(norm=norm_q, cmap=cmap_q), cax=axes[0][1],
                    orientation="horizontal", label=r"q [$\AA^{-1}$]")
if args.select_q!=[]: add_upperxticks_to_cbar(cbar_q,q_selected)
if args.xlog==1:
    ax.set_xscale("log")
plt.tight_layout()

fig.savefig(outname+".png")
fig.savefig(outname+".pdf")
print("Figure saved on %s.png , %s.pdf\n"%(outname, outname))
#plt.show()
