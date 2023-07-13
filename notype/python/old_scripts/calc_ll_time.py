import numpy as np
from matplotlib import pyplot as plt
import sys, os, fnmatch

def find_files(directory, pattern): #Generator: finds all coordinate files and their respective N
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                Nval = int(''.join(filter(str.isdigit, basename)))
                yield Nval,filename
#######################################################################################
def get_statistics(source):
    with open(source,"r") as f:
        lines = f.readlines()
        n=len(lines)
        t0=np.zeros(n-1)#skip header
        t1=np.zeros(n-1)
        for i in range(1,n): #skip header
            words = lines[i].split()
            t0[i-1] = float(words[1])
            t1[i-1] = float(words[2])
    return t0.mean(),t0.std(), t1.mean(),t1.std()

Ns=[]
mT0s=[] # mT0: average time with no linked list
sT0s=[] # sT0: std.dev. of ''
mT1s=[] # mT1: average time with linked list
sT1s=[] # ...
for N,filename in find_files('./', 'll_test_N*.dat'):
    Ns.append(N)
    mT0,sT0, mT1,sT1 = get_statistics(filename)
    mT0s.append(mT0)
    sT0s.append(sT0)
    mT1s.append(mT1)
    sT1s.append(sT1)
NN=len(Ns)
zipped = list(zip(Ns,mT0s,sT0s,mT1s,sT1s)) #Ns go first!
zipped.sort()

totlist = list(zip(*zipped))
np.savetxt("ll_test_results.dat", np.transpose(totlist), header="# 1:N 2:avg_t_noLL 3:std_t_noLL 4:avg_t_LL 5:std_t_LL")
N,mT0,sT0,mT1,sT1 = totlist

fig,ax = plt.subplots(1,1)
ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')
plt.errorbar(N,mT0,sT0,fmt="x--", label="no LL")
plt.errorbar(N,mT1,sT1,fmt="x--", label="LL")
ax.set_xlabel("N")
ax.set_ylabel("t (ms)")
#ax.set_title("Errore stimato da 10 ripetizioni")
plt.gca().set_aspect('equal')
ax.set_ylim(bottom=0.1)
fig.tight_layout()
plt.legend()
plt.savefig("ll_test_results.png")
plt.show()
