#!/usr/bin/env python
import sys # for argv
from subprocess import call #for shell usage
import numpy as np
import fileinput # for file lines overwriting
import matplotlib.pyplot as plt
from math import floor,pi

import os, fnmatch, re

def find_files(directory, pattern): #Generator: finds all coordinate files and their respective step
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                stepval = int(''.join(filter(str.isdigit, basename)))
                yield stepval,filename
        break #prevents going into subfolders

#######################################################
def get_r_molgl(line): #get r from molgl line
    words = line.split()
    return np.array(words[:3], dtype=float)
######################################################
def get_rs(source): #get all r's from molgl file
    with open(source,"r") as f:
        lines = f.readlines()
        N = len(lines)
        rs = np.zeros((N,3))
        for i in range(N):
            rs[i,:] = get_r_molgl(lines[i])
    return rs

##########################################################################################
# Calcola il MSD di tutte le configurazioni che trova nella cartella assegnata

debug=0
start_step=-1
end_step=sys.maxsize * 2 + 1

if len(sys.argv)<2:
    sys.exit("\nUSAGE: python "+sys.argv[0]+" <folder path> [start_step] [end_step]\n");
else:
    directory=sys.argv[1]
    if len(sys.argv)>2:
        start_step=int(sys.argv[2])
        if len(sys.argv)>3:
            end_step=int(sys.argv[3])

N=int(''.join(filter(str.isdigit, directory.split('/')[1].split('N')[-1])))

steps=[]
fnames=[]
Ls=[]
Nsteps=0
for step,filename in find_files(directory, 'mcsim_coords_s*'):
    if step>=start_step and step <= end_step:
        Nsteps=Nsteps+1
        fnames.append(filename)
        steps.append(step)
        with open(directory+"mcsim_state_s"+filename.split("_s")[1],"r") as f:
            lines=f.readlines()
            for line in lines:
                words=line.split(' ')
                if words[0]=='L':
                    Ls.append(float(words[1]))
zipped = list(zip(steps,fnames,Ls)) #steps go first!
zipped.sort()
steps.sort()
print("Total n. of files:",Nsteps)
dts=[]
for i in range(Nsteps):
    for j in range(i):
        dts.append(steps[i]-steps[j])
dts=np.sort(np.array(list(set(dts)))) #remove duplicates and sort
msd=np.zeros_like(dts)
msd2=np.zeros_like(dts)
Navg=np.zeros_like(dts)
r=np.zeros((Nsteps,N,3)) #all the trajectories

print("\nReading and averaging (r(t)-r(t0))^2 over particles,coordinates and t0...")

for i in range(Nsteps):#step,fname,L in zipped:
    step,fname,L = zipped[i]
    print("step=%10d L=%f"%(step,L))
    r[i,:,:] = get_rs(fname)
    if i>0:
        dr = r[i,:,:] - r[i-1,:,:]
        dr = (dr-L*(dr/L).round()) #PBC
        r[i,:,:] = r[i-1,:,:] + dr
        for j in range(i):
            dt=steps[i]-steps[j]
            idx=np.where(dts==dt)[0][0]
            dr = r[i,:,:] - r[j,:,:]
            dr2=dr**2
            dr4=dr2**2
            msd[idx] += dr2.mean()
            msd2[idx] += dr4.mean()
            Navg[idx] += 1

msd=msd/Navg
msd2=msd2/Navg
msd2=np.sqrt(msd2-msd**2)/N/3
msd2[Navg>1]=msd2[Navg>1]/np.sqrt(Navg[Navg>1]-1)
#np.savetxt(directory+"msd.dat", np.transpose([dts,msd,msd2]))
print("done")
D=(msd[-1]-msd[0])/(dts[-1]-dts[0])
print("Stima grezza (pendenza globale): D=%f\n"%D)
fig=plt.figure()
ax=plt.axes()
ax.errorbar(dts,msd,msd2,fmt='x-')
#plt.axhline(4.0,label="$\sigma^2$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("t")
ax.set_ylabel("MSD")
ax.set_aspect("equal")
ax.grid(True, which='both')
ax.tick_params(which='both', direction='in')
plt.tight_layout()
plt.savefig("msd.eps")

if debug:
    plt.figure()
    plt.hist(dr_raw[:,0], label="raw $\Delta x$", alpha=0.5)
    plt.hist(dr[:,0], label="PBC $\Delta x$",alpha=0.5)
    plt.title("Time %d"%t[-1])
    plt.legend()

plt.show()

#################################################
def plot_MSD(t,msd,phi,L,s):
    #plot data
    plt.figure()
    plt.grid(True)
    plt.plot(t-t[0],msd,"x-")
    plt.xlabel("t")
    plt.ylabel("MSD(t)")
    tit ="$\phi$="+str(phi)+" L="+str(L)
    plt.title(tit)
    plt.savefig(fname("MSD/msd",L,N,s,".pdf"))
    plt.show()
    plt.close()

