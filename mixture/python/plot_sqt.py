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
                    description = 'Plots the Angular-Limited Three-Body Correlation as a 2D heatmap.',
                    epilog='End of the summary.'
)
parser.add_argument('--inavg',  type=argparse.FileType('r'),
                     default="sqt.ave", required=False,
                     help="Input file with average trajectory (1st block: wave vector bis; 2nd block: Delta timestep bins; 3rd block: 2d matrix of <S(q,t)>; 4th block; 2d matrix of its error. [default: %(default)s]"
)
parser.add_argument('--dt',  type=float, default=0.002, required=False,
                     help="Integration time step (picoseconds). [default: %(default)s]"
)

parser.add_argument('--yshift',  type=float, default=0.0, required=False,
                     help="Shift vertically the vs-time-plot. [default: %(default)s]"
)
parser.add_argument('--normalize',  type=int, default=1, required=False,
                     help="Normalize the vs-time-plot to S(q,0)? 1: yes ; 0: no. [default: %(default)s]"
)


outname="sqt"

#-------------------------------------#
args = parser.parse_args()
assert args.normalize==0 or args.normalize==1

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

#-------------------------------------------------------------------------------------------------------------------#

fig, axes = plt.subplots(2,2, figsize=(10,5),gridspec_kw={'height_ratios': [1, 9]} ) # figsize is experimental
ax = axes[1][0]
ax.set_xlabel(r"$q$ $[\AA^{-1}]$")
ax.set_ylabel(r"$S(q,t)$")
cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=t.min(), vmax=t.max())
colors = cmap(norm(t))
for i in range(Nt):
    lab = r"$t=%.1f$ $ps$"%t[i]
    col = colors[i]
    y = sqt[:,i] + i*args.yshift
    y_ = sqt_[:,i]
    ax.errorbar(q,y,y_, fmt='.-', color=col, label=lab, alpha=0.5)
ax.tick_params(which='both', direction='in')
ax.set_yscale("log")
#ax.grid(axis='both', which='major')
fig.colorbar( mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=axes[0][0], orientation="horizontal", label="t [ps]")

ax = axes[1][1]
ax.set_xlabel(r"$t$ [$ps$]")
if args.normalize==1:
	ax.set_ylabel(r"$S(q,t)$ / $S(q,0)$")
else:
	ax.set_ylabel(r"$S(q,t)$")
cmap = mpl.cm.magma
norm = mpl.colors.Normalize(vmin=q.min(), vmax=q.max())
colors = cmap(norm(q))
for i in range(Nq):
    lab = r"$q=%.2f$ $\AA^{-1}$"%q[i]
    col = colors[i]
    if args.normalize==1:
        y = sqt[i] / sqt[i,0] + i*args.yshift
        y_ = y * np.sqrt( (sqt_[i]/sqt[i])**2 + (sqt_[i,0]/sqt[i,0])**2 ) # gaussian error propagation
    else:
        y = sqt[i] + i*args.yshift
        y_ = sqt_[i]
    ax.errorbar(t,y,y_, fmt='.-', color=col, label=lab, alpha=0.5)
ax.tick_params(which='both', direction='in')
#ax.grid(axis='both', which='major')
fig.colorbar( mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=axes[0][1], orientation="horizontal", label=r"q [$\AA^{-1}$]")

plt.tight_layout()

fig.savefig(outname+".png")
fig.savefig(outname+".pdf")
print("Figure saved on %s.png , %s.pdf\n"%(outname, outname))
#plt.show()
