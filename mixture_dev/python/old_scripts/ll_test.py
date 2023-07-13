#!/usr/bin/env python

from subprocess import call
import fileinput
import numpy as np
import matplotlib.pyplot as plt
from math import floor,pi

#######################################################################
def fname(string,L,N,s,tag=".dat"):
    return string+"_L"+str(L)+"_N"+str(N)+"_s"+str(s)+tag
####################################################################
def set_conf_file(N,L,delta,lt=-1,fname="mcsim_config.dat"):
    for line in fileinput.input(fname, inplace=True):
        words=line.split()
        if words[0]=="N":
            print("N",N)
        elif words[0]=="L":
            print("L",L)
        elif words[0]=="del_tra":
            print("del_tra",delta)
        elif words[0]=="del_rot":
            print("del_rot",delta)
        elif words[0]=="del_box":
            print("del_box",delta)
        elif words[0]=="latt_type":
            print("latt_type",lt)
        else:
            print(line,end="")

phi = 0.1
N=int(8)
Nmax=1e4
p=str(phi)
nmax=str(Nmax)
while N<Nmax:
    n=str(N)
    L = (N*(4*pi/3 *2)/phi)**(1.0/3.0)
    l=str(L)
    set_conf_file(N,L,1,0)
    print("Running "+n+" out of "+nmax+" phi="+p+" L="+l)
    call("./main.exe -v -tll 0 "+n+" "+n, shell=True)
    N = N*2
