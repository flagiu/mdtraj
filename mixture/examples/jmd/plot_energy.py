#!/bin/python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use("~/.matplotlib/my.mplstyle")
plt.rcParams['lines.linewidth'] = 0.1
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=["k", "r", "g", "b", "m", "c", "y"]) 
import jmd_info as ji

m=121.75 #g/mol #MONOSPECIE
bar_to_GPa = 1e-4
eVang3_to_GPa = 160.218

N,V_initial = ji.get_NV()
T = ji.get_T()
P = ji.get_P()
#rho=10/6.023*m*N/V  # g/mol , ang^3 --> g/cm^3
T/=8.6167e-05       # eV --> K
P*=eVang3_to_GPa    # eV/ang^3 --> GPa
title = r"N=%d , T=%.1f K , $P$=%.1f GPa"%(N,T,P)

ref_dt  =2e-6 #ns
inf_dt  =0.2*10.179e-3*1e-3 #ns
outname="energy"

#-------------------------------------#

fig, axes = plt.subplots(4,1,dpi=300, sharex=True, figsize=(6,6))
axes[0].set_ylabel(r"$E_{pot}$ [eV/at.]")
axes[1].set_ylabel(r"$P$ [GPa]")
axes[2].set_ylabel(r"$\rho$ [g/cm$^3$]")
axes[3].set_ylabel(r"$T$ [K]")
axes[-1].set_xlabel(r"$t$ [ns]")
fig.suptitle(title)

## Reference
x,y = np.loadtxt("REFERENCE/energy_potential.dat", unpack=True) # just to get data count
skiprows_ref=int(0.01*len(x))
x,y = np.loadtxt("REFERENCE/energy_potential.dat", unpack=True, skiprows=skiprows_ref) # timestep, eV
t = (x-x[0])*ref_dt # ns
y = y/N # eV --> eV/at.
axes[0].plot(t, y, '-', alpha=0.7, linewidth=1.0, label="Reference")

x,y = np.loadtxt("REFERENCE/pressure.dat", unpack=True, skiprows=skiprows_ref) # timestep, kbar
t = (x-x[0])*ref_dt # ns
y = y*bar_to_GPa # bar --> GPa
axes[1].plot(t, y, '-', alpha=0.7, linewidth=1.0, label="Reference")

x,y = np.loadtxt("REFERENCE/density.dat", unpack=True, skiprows=skiprows_ref) # timestep, g/cm3
t = (x-x[0])*ref_dt # ns
axes[2].plot(t, y, '-', alpha=0.7, linewidth=1.0, label="Reference")

x,y = np.loadtxt("REFERENCE/temperature.dat", unpack=True, skiprows=skiprows_ref) # timestep, Kelvin
t = (x-x[0])*ref_dt # ns
axes[3].plot(t, y, '-', alpha=0.7, linewidth=1.0, label="Reference")
###

### Inference
x,y = np.loadtxt("energy_potential.dat", unpack=True) # timestep, eV/at.
t = (x-x[0])*inf_dt # ns
axes[0].plot(t, y, '-', alpha=0.7, linewidth=1.0, label="Inference")

x,y = np.loadtxt("pressure.dat", unpack=True) # timestep, eV/ang^3
t = (x-x[0])*inf_dt # ns
y = y*eVang3_to_GPa # eV/ang^3 --> GPa
axes[1].plot(t, y, '-', alpha=0.7, linewidth=1.0, label="Inference")

x,y = np.loadtxt("density.dat", unpack=True) # timestep, atoms/ang^3
t = (x-x[0])*inf_dt # ns
y = y*10/6.023*m # atoms/ang^3 --> g/cm3
axes[2].plot(t, y, '-', alpha=0.7, linewidth=1.0, label="Inference")

x,y = np.loadtxt("temperature.dat", unpack=True) # timestep, atoms/ang^3
t = (x-x[0])*inf_dt # ns
y = y/8.6167e-5 # eV --> Kelvin
axes[3].plot(t, y, '-', alpha=0.7, linewidth=1.0, label="Inference")
###

axes[0].legend()
for ax in axes: ax.tick_params(which='both', direction='in')
#ax.grid(axis='both', which='major')
plt.tight_layout()
fig.savefig(outname+".png")
print("Figure saved on %s.png \n"%(outname))
#plt.show()
