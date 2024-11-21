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
sys.path.append("my_path")
from read_sqt import read_sqt
from read_lammps_output import read_from_lammps_output
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp1d


parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Plots the self-overlap parameter Q(t) and its susceptibility.',
                    epilog='as a function of t in log scale.'
)
parser.add_argument('--files',  type=str, nargs="+",
                     default=["Qoverlap.ave"], required=False,
                     help="Input file(s) with average trajectory (Delta timestep bins; Qs Qss Qd Qdd Qsd . [default: %(default)s]"
)
parser.add_argument('--norm',  type=bool,
                     default=False, required=False,
                     help="Normalize to V/kT ?. [default: %(default)s]"
)

cmap_files = mpl.cm.coolwarm
outname="Qoverlap"
#taufile=".dat"

# Physical constants in conventional units ---------------#
Navog=6.02214076e23 #atoms/mol
kB=1.380649e-23 #J/K     #8.61733326e-5 #eV/K
e=1.60217663e-19 # J/eV

# Derived constants --------------------------------------#
kB_eVK=kB/e         # kB is useful also in eV/K

def smooth(x, n=3):
    w = np.ones(n)/n
    return np.convolve(x,w, mode="same")

#-------------------------------------#
args = parser.parse_args()

fig, axes = plt.subplots(4,1,sharex=True, figsize=(8,10)) # figsize is experimental
ax = axes[0]
if args.norm:
	ax.set_ylabel(r"$\chi_4^{ss}(t)$ [nm$^3$/eV]") # = \beta V \langle \delta q(t)^2 \rangle
else:
	ax.set_ylabel(r"$\langle \delta q^{ss}(t)^2 \rangle$")
#ax.set_ylabel(r"$\langle q(t) \rangle$")
#ax.axhline(np.exp(-1),xmin=0, color="k", linestyle="--", label=r"$1/e$", linewidth=0.5, zorder=-999)
ax = axes[1]
if args.norm:
	ax.set_ylabel(r"$\chi_4^{dd}(t)$ [nm$^3$/eV]") # = \beta V \langle \delta q(t)^2 \rangle
else:
	ax.set_ylabel(r"$\langle \delta q^{dd}(t)^2 \rangle$")
#ax.set_ylabel(r"$\langle q^2(t) \rangle$")
#ax.axhline(np.exp(-1),xmin=0, color="k", linestyle="--", label=r"$1/e$", linewidth=0.5, zorder=-999)
ax = axes[2]
if args.norm:
	ax.set_ylabel(r"$\chi_4^{sd}(t)$ [nm$^3$/eV]") # = \beta V \langle \delta q(t)^2 \rangle
else:
	ax.set_ylabel(r"$\langle \delta q^{sd}(t)^2 \rangle$")
#ax.set_ylabel(r"$\langle q(t) \rangle\cdot \langle q(t) \rangle$")
#ax.axhline(np.exp(-1),xmin=0, color="k", linestyle="--", label=r"$1/e$", linewidth=0.5, zorder=-999)
ax = axes[3]
if args.norm:
	ax.set_ylabel(r"$\chi_4(t)$ [nm$^3$/eV]") # = \beta V \langle \delta q(t)^2 \rangle
else:
	ax.set_ylabel(r"$\langle \delta q(t)^2 \rangle$")

axes[-1].set_xlabel(r"$t$ $[ps]$")
#----------------------------------------------------------#

for i,file in enumerate(args.files):
    try:
        timesteps,Qs,Qss,Qd,Qdd,Qsd=np.loadtxt(file, unpack=True)
    except FileNotFoundError:
        print("  Warning: file %s not found!"%file)
        continue
    except ValueError:
        print("  Warning: file %s empty or corrupted!"%file)
        continue

    col = cmap_files(i/(len(args.files)-1) if len(args.files)>1 else 0.0)

    T_nvt = float(file.split("K")[0])
    filepath="/".join(file.split("/")[:-1])
    try:
        thermo_file=filepath+"/../out.lmp"
        box, N, dt, thermo_names, thermo = read_from_lammps_output(thermo_file)
        V = box[0]*box[3]*box[5] if len(box)==6 else box[0]*box[1]*box[2]
        temp = np.array(thermo["Temp"])
        T_nve=temp.mean()
        print("  T_nvt = %.2f K ; T_nve = %.2f K averaged on all thermo data"%(T_nvt,T_nve))
        lab = r"$T_{NVT}=%.0f$K , $T_{NVE}=%.0f$K"%(T_nvt,T_nve)
    except FileNotFoundError:
        print("  Warning: file %s not found!"%thermo_file)
        dt=args.dt
        T_nve=T_nvt
        lab = r"$T_{NVT}=%.0f$K , $T_{NVE}=T_{NVT}$"%(T_nvt)
        N = 4096
        V = 51.3388**3 # angstrom^3
    
    if args.norm:
        chifactor = (V * 1e-3)/(kB_eVK*T_nve)
    else:
        chifactor = 1.0

    #-------------------------------------------------------------------------------------------------------------------#
    
    t=timesteps*dt
    Q=Qs+Qd
    Q2=Qss+Qdd+Qsd
    chi4 = (Q2-Q*Q) * chifactor
    
    chi4ss = (Qss-Qs*Qs) * chifactor
    chi4dd = (Qdd-Qd*Qd) * chifactor
    chi4sd = (Qsd-2*Qs*Qd) * chifactor
    
    """
    axes[0].plot(t,Qs, "-", label="s", color=col, zorder=-i)
    axes[0].plot(t,Qd, "--", label="d", color=col, zorder=-i)
    
    axes[1].plot(t,Qss, "-", label="ss", color=col, zorder=-i)
    axes[1].plot(t,Qdd, "--", label="dd", color=col, zorder=-i)
    axes[1].plot(t,Qsd, ":", label="sd", color=col, zorder=-i)
    
    axes[2].plot(t,Qs*Qs, "-", label="s*s", color=col, zorder=-i)
    axes[2].plot(t,Qd*Qd, "--", label="d*d", color=col, zorder=-i)
    axes[2].plot(t,Qs*Qd+Qd*Qs, ":", label="s*d + d*s", color=col, zorder=-i)
    #axes[2].plot(t,Qs*Qs+Qd*Qd+2*Qs*Qd, "x-", label="total (sum)", color=col, zorder=-i)
    
    axes[3].plot(t,chi4ss, "-", color=col, label="ss", zorder=-i)
    axes[3].plot(t,chi4dd, "--", color=col, label="dd", zorder=-i)
    axes[3].plot(t,chi4sd, ":", color=col, label="sd+ds", zorder=-i)
    #axes[3].plot(t,chi4ss+chi4dd+chi4sd, "x-", color=col, label="total (sum)", zorder=-i)
    """
    
    #axes[0].plot(t,Q, ".-", label=lab, color=col, zorder=-i)
    #axes[1].plot(t,Q2, ".-", label=lab, color=col, zorder=-i)
    #axes[2].plot(t,Q*Q, ".-", label=lab, color=col, zorder=-i)
    axes[0].plot(t,chi4ss, ".-", label=lab, color=col, zorder=-i)
    axes[1].plot(t,chi4dd, ".-", label=lab, color=col, zorder=-i)
    axes[2].plot(t,chi4sd, ".-", label=lab, color=col, zorder=-i)
    axes[3].plot(t,chi4, ".-", label=lab, color=col, zorder=-i)
    
    #axes[3].plot(t,2+smooth(chi4self,n=3), args.fmt, color=col, label=None, zorder=-i)


#ax.grid(axis='both', which='major')
for ax in axes:
    ax.tick_params(which='both', direction='in')
    ax.set_xscale("log")
    ax.set_xlim((1e-1,1e3)) # picoseconds

axes[-1].legend(frameon=False, fontsize=6)

#for ax in axes[:3]:
#    ax.set_ylim(0,1.2)

fig.tight_layout()

fig.savefig(outname+".png")
fig.savefig(outname+".pdf")
print("Figure saved on %s.png , %s.pdf\n"%(outname, outname))
#plt.show()
