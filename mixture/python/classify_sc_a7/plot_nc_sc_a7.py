#!/bin/python3
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.autolayout'] = True
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14

q4dot_th=0.65  # threshold for crystallinity
class_dict={"SC":0,"A7":1} # classes assigned by the classifier
type_lmp_dict={"SC":2,"A7":3} # type assigned to LAMMPS output
dt=0.002 #ps
nc_q4dot_file="nc.l4.dat"
sc_a7_file="sc_a7.dat"
sc_a7_ave_file="sc_a7_ave.dat"

steps,nc=np.loadtxt(nc_q4dot_file, unpack=True)
#_,sc,a7,outlier=np.loadtxt(sc_a7_file, unpack=True)
#assert (steps==_).all()
_,sc_ave,a7_ave,outlier_ave=np.loadtxt(sc_a7_ave_file, unpack=True)
assert (steps==_).all()
t=steps*dt

plot_args=dict(markerfacecolor="none",linewidth=1,alpha=0.9)
fig,ax = plt.subplots(figsize=(6,6))
ax.set_xlabel("t [ps]")
ax.set_ylabel("Labeled particles")
ax.tick_params(which="both",direction="in")
ax.plot(t,nc,'k-',label="XTAL",**plot_args)
#ax.plot(t,sc,'b--',label="SC 1",**plot_args)
#ax.plot(t,a7,'r--',label="A7 1",**plot_args)
#ax.plot(t,outlier,'g--',label="Outlier 1",**plot_args)
ax.plot(t,sc_ave,'b-',label="SC",**plot_args)
ax.plot(t,a7_ave,'r-',label="A7",**plot_args)
ax.plot(t,outlier_ave,'g-',label="Outlier",**plot_args)
ax.set_ylim((0,ax.get_ylim()[1])
ax.legend(frameon=False)
fig.savefig("sc_a7.png")
