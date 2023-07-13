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
#######################################################
def get_r_molgl(line): #get r from molgl line
    words = line.split()
    return np.array(words[:3], dtype=float)
######################################################
def get_min_dist(source,L):
    minR = 2*L
    with open(source,"r") as f:
        lines = f.readlines()
        N = len(lines)
        for i in range(N):
            ri = get_r_molgl(lines[i])
            for j in range(i+1,N):
                rj = get_r_molgl(lines[j])
                rij = ri - rj
                rij -= L*np.array([floor(rij[0]/L),
                                   floor(rij[1]/L),
                                   floor(rij[2]/L)])
                rij = np.linalg.norm(rij)
                if rij<minR:
                    minR = rij
    return minR
######################################################            
def calc_g_r(source,L):
    err = False
    b = 1/10 #(0.5*L)/100
    Nb = floor(0.5*L/b)
    rs = np.zeros(Nb)
    bins = np.zeros(Nb)
    with open(source,"r") as f:
        lines = f.readlines()
        N = len(lines)
        for i in range(N):
            ri = get_r_molgl(lines[i])
            for j in range(i+1,N):
                rj = get_r_molgl(lines[j])
                rij = ri-rj
                rij -= L*np.array([floor(rij[0]/L),
                                   floor(rij[1]/L),
                                   floor(rij[2]/L)])
                rij = np.linalg.norm(rij)
                k = floor(rij/b)
                if (k<Nb):
                    bins[k] += 1;
                    if (rij<2.0):
                        print("Overlapping: ",i,j,ri,rj,rij)
                        err=True
                        break
    c0 = (N-1)/(L*L*L) * 4.0*pi/3.0 * N/2.0 / 8.0 # why do I need 8?
    v1 = 0.0
    for k in range(Nb):
        rs[k] = (k+0.5)*b
        v2 = ((k+1)*b)*((k+1)*b)*((k+1)*b)
        bins[k] /= c0*(v2-v1) # normalize to ideal gas g(r)
        v1 = v2
    return rs,bins,err
#################################################
def plot_g(r,g,gerr,L,s):
    #plot data
    plt.figure()
    plt.grid(True)
    plt.errorbar(r,g,gerr,marker=".")
    plt.xticks(np.arange(0,L/2+1, 1.0))
    plt.ticks_params(which='both', direction='in')
    plt.xlabel("r")
    plt.ylabel("g(r)")
    if s<0:
        figname=directory+"g/g(r)_average.png"
    else:
        figname=directory+"g/g(r)_s"+str(s)+".png"
    plt.savefig(figname)
    plt.show()
    plt.close()
#################################################
def get_g(path,L,step, read=0):
    try:
        r,g,zeros = np.loadtxt(directory+"g/g(r)_s"+str(step)+".dat", unpack=True)
        return (r,g,zeros),False
    except:
        r,g,overlap = calc_g_r(path,L)
        zeros = np.zeros_like(g)
        if overlap:
            print("Overlapping occurred at s="+str(step)+". g(r) will not be saved.\n")
        else:
            np.savetxt(directory+"g/g(r)_s"+str(step)+".dat", np.array([r,g,zeros]).T)
        return (r,g,zeros),overlap
    
###################################################
##########    MAIN    #############################
###################################################

only_avg=0
if len(sys.argv)<2:
    sys.exit("USAGE: python "+sys.argv[0]+" <folder path> [only_average=0]\n\nQuesto programma calcola la g(r) di ogni configurazione salvata nella cartella e ne fa la media.\nSe only_average=1, fa solo la media delle g(r) già salvate.\n");
else:
    directory=sys.argv[1]
    if len(sys.argv)>2:
        only_avg=int(sys.argv[2])
try:
    os.mkdir(directory+"g/")
except FileExistsError:
    print("g-folder already exists.")
L=int(''.join(filter(str.isdigit, directory.split('/')[1].split('N')[0])))
N=int(''.join(filter(str.isdigit, directory.split('/')[1].split('N')[-1])))

fnames=[]
steps=[]
for step,filename in find_files(directory, 'mcsim_coords_s*'):
    if step>0:
        fnames.append(filename)
        steps.append(step)
Nfiles=len(fnames)
Navg=Nfiles
zipped = list(zip(steps,fnames)) #steps go first!
zipped.sort()
first=True
i=0
for step,fname in zipped:
    if (i%10 == 0):
        print("> done "+str(i)+" / "+str(Nfiles))
    i=i+1
    (r,g,zeros),overlap = get_g(fname, L, step, read=only_avg)
    if first:
        first=False
        G=np.zeros(len(g))
        G2=np.zeros(len(g))
    if overlap:
        Navg=Navg-1
    else:
        G += g
        G2 += g**2
G /= Navg # <g>
G2 /= Navg # <g^2>
G2 = np.sqrt((G2 - G**2)/Navg) # st.dev. of <g>
np.savetxt(directory+"g/g(r)_average_n"+str(Navg)+".dat", np.array([r,G,G2]).T)
plot_g(r,G,G2,L,-1)
print("g(r)'s and their average were saved. Check the g-folder.")

#print("Problematic rmin: ",rmin[rmin<2.0])
