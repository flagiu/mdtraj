#!/bin/python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'
plt.rcParams['figure.dpi'] = 200
plt.rcParams['figure.autolayout'] = True

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Plots each trajectory r^2(t) and the average MSD=<r^2(t)>.',
                    epilog='End of the summary.'
)
parser.add_argument('--inavg',  type=argparse.FileType('r'),
                     default="msd.ave", required=False,
                     help="Input file with average trajectory (columns:t, MSD for each type, MSD error for each type, MSD of Center-of-Mass). [default: %(default)s]"
)
parser.add_argument('--dt',  type=float, default=0.002, required=False,
                     help="Integration time step (picoseconds)"
)
parser.add_argument('--inlabels',  type=argparse.FileType('r'),
                     default="labels.dat", required=False,
                     help="Input file with one row per type, two columns: atom label, occurrence. [default: %(default)s]"
)
parser.add_argument('--ignore',  type=str, nargs="*",
                     default=[], required=False,
                     help="Ignore the following atom types in the plot. INPUT: one or more strings (0,1,...,ntypes-1, or labels). [default: %(default)s]"
)
parser.add_argument('--fitD', type=bool, required=False, default=False, action=argparse.BooleanOptionalAction,
                    help="Do a linear fit to deduce the diffusion coefficient D. Remember to check dt! [default: %(default)s]"
)

parser.add_argument('--outname',  type=str,
                     default="msd", required=False,
                     help="Prefix for the output file. [default: %(default)s]"
)

nminD=int(3)

args = parser.parse_args()

print("Plotting MSD average ...")
#--------------------------------------------------------------------------------------------------------#
labels=[]
Nt=[]
lab2idx = {}
lines = args.inlabels.readlines()
for i,line in enumerate(lines):
    lab, nt = line.strip('\n').split()
    labels.append(lab)
    Nt.append(int(nt))
    lab2idx[lab] = i
Nt=np.array(Nt)
xt=Nt/np.sum(Nt) #Fraction
print("  Atom types:",*labels)
print("  Occurrence:",*Nt)
print("  Fraction:",*xt)
ign_labels=np.array(args.ignore, dtype=str)
ign_types=np.zeros(len(ign_labels), dtype=int)
print("  List of types to be ignored:",ign_labels)
#---------------------------------------------------------------------------------------------------------#

X = np.loadtxt(args.inavg, unpack=False).T
t = args.dt*X[0]
assert (X.shape[0]-2)//2==len(Nt) # number of types

fig, ax = plt.subplots()
ax.set_xlabel("t [ps]")
ax.set_ylabel(r"$\langle\, \Delta r^2(t)\, \rangle$ [$\AA^2$]")

if args.fitD:
    figD, axD = plt.subplots()
    axD.set_xlabel(r"Fitted on last $n$ points")
    axD.set_ylabel(r"$D$ [$\AA^2/ps$]")
    Ddata=np.empty((len(t)-nminD+1, 1+2*len(Nt)))
    Ddata[:,0]=np.arange(len(t),nminD-1,-1)

for i in range(len(Nt)):
    msd=X[1+i]
    msd_=X[1+len(Nt)+i]
    if args.fitD:
        print("# Progressive linear fit (no intercept) on the last n points:")
        def f(x,a=1): return a*x
        p0 = (msd[-1]-msd[0])/(t[1]-t[0])
        print("# Type",labels[i])
        print("# n, D, D_error")
        for ni,n in enumerate(Ddata[:,0]):
            n=int(n)
            x=t[-n:]
            y=msd[-n:]
            y_=msd_[-n:]
            popt,pcov = curve_fit(f, x,y, sigma=y_, p0=[p0], maxfev=10000)
            popt_ = np.sqrt(np.diag(pcov))
            Ddata[ni,1+i] = popt[0]
            Ddata[ni,1+len(Nt)+i] = popt_[0]
            print("%d %f %f"%(n,popt[0],popt_[0]))
            xl=np.linspace(0.9*x[0],1.1*x[-1],3)
            #ax.plot(xl,f(xl,*popt),'k--',alpha=0.7)
        print()
        axD.errorbar(Ddata[:,0],Ddata[:,1+i],Ddata[:,1+len(Nt)+i], fmt="s-", alpha=0.7, label=labels[i])
    ax.errorbar(t, msd, msd_, fmt="o-", alpha=0.5, label=labels[i] )

ax.legend()
ax.tick_params(which='both', direction='in')
#ax.grid(axis='both', which='major')
fig.savefig(args.outname+".eps")
fig.savefig(args.outname+".png")
print("Figure saved on %s.png , %s.eps\n"%(args.outname,args.outname))

if args.fitD:
    fmt="%.0f"
    for i in range(len(Nt)): fmt+=" %f %f"
    np.savetxt(args.outname+"_D.dat", Ddata, header="n (last points for fit), D for each type, D error for each type", fmt=fmt)
    print("Diffusion Coefficient data saved on %s_D.dat\n"%(args.outname))
    axD.legend()
    axD.tick_params(which='both', direction='in')
    figD.savefig(args.outname+"_D.eps")
    figD.savefig(args.outname+"_D.png")
    print("Figure saved on %s_D.png , %s_D.eps\n"%(args.outname,args.outname))

#plt.show()
