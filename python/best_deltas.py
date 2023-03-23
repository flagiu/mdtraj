#!/usr/bin/env python

from subprocess import call
import numpy as np
import fileinput
import matplotlib.pyplot as plt
from math import floor,pi

# find the best deltas for translational, rotational and box MC moves at fixed density
# for uniaxial hard ellipsoids with X0=2.
# in units of the shorter semiaxis
#######################################################################
def fname(string,L,N,s,tag=".dat"):
    return string+"_L"+str(L)+"_N"+str(N)+"_s"+str(s)+tag
#######################################################################à
def set_conf_file(N,L,delta,lt=-1,fname="mcsim_config.dat"):
    for line in fileinput.input(fname, inplace=True):
        words=line.split()
        if words[0]=="N":
            print("N",N)
        elif words[0]=="L":
            print("L",L)
        elif words[0]=="latt_type":
            print("latt_type",lt)
        else:
            print(line,end="")
########################################################################
def get_rates(source):
    with open(source,"r") as f:
        for line in f:
            words=line.split()
            if words[0]=="tra_acc_freq":
                x = float(words[1])
            elif words[0]=="rot_acc_freq":
                y = float(words[1])
            elif words[0]=="box_acc_freq":
                z = float(words[1])
    return x,y,z
######################################################

phis = [0.1,0.2,0.3]
Ls = [20,40]#,80]#,160,320]
nsteps=10000
savesteps=nsteps
deltas = np.linspace(2.5,5, 12)
print("deltas: ",deltas)

for phi in phis:
    for L in Ls:
        N = floor(phi*L*L*L/(4*pi/3 *2))
        p=str(phi)
        l=str(L)
        tra=[]
        rot=[]
        box=[]
        set_conf_file(N,L,0,1,savesteps,0) # initialize the particles once for all
        call("./main.exe", shell=True)
        for i in range(len(deltas)):
            call("echo 'Running "+str(i+1)+" out of "+str(len(deltas))+" phi="+p+" L="+l+" delta="+str(deltas[i])+"'", shell=True)
            set_conf_file(N,L,1,nsteps-1,savesteps,deltas[i])
            call("./main.exe", shell=True)
            
            x,y,z = get_rates(fname("NVT/data/mcsim_state",L,N,nsteps))
            tra.append(x)
            rot.append(y)
            box.append(z)

        tra = np.array(tra)
        rot = np.array(rot)
        box = np.array(box)
        np.savetxt("NVT/deltas/deltas_phi"+p+"_L"+l+"_nsteps"+str(nsteps)+".dat",
                   np.array([deltas,tra,rot,box]).T.round(3),
                   header="# 1:delta 2:tra_acc_rate 3:rot_acc_rate 4: box_acc_rate\n")
        plt.figure()
        plt.plot(deltas,tra,"-o",label="tra")
        plt.plot(deltas,rot,"-o",label="rot")
        #plt.plot(deltas,box,"-o",label="box")
        plt.xlabel("delta")
        plt.ylabel("acceptance rate")
        plt.title("phi="+p+" L="+l)
        plt.legend()
        plt.savefig("NVT/deltas/deltas_phi"+p+"_L"+l+" nsteps"+str(nsteps)+".pdf")
        plt.close()
