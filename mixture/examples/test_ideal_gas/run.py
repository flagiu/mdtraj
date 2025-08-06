#!/bin/python3
import os
import numpy as np
from ase import Atoms
import matplotlib.pyplot as plt
np.random.seed(42)

nframes=10
N=1000
dim=int(3)
rho=1 #atoms/volume
outname="pos.xyz"

V=N/rho
L=V**(1.0/dim) # cubic box
print("Generating ideal gas trajectory")
os.system(f"[ -f {outname} ] && rm {outname}")
for i in range(nframes):
    pos=np.random.uniform(0,L, size=(N,dim))
    atoms = Atoms(('H',)*N, positions=pos, cell=(L,)*dim, pbc=(True,True,True))
    atoms.write(outname, format='xyz', append=True, comment=f"Timestep = {i}")

# Analyze
mdtraj="/home/flavio/programmi/mdtraj/mixture/bin/mdtraj"
#print("Computing g(r)")
#os.system(f"{mdtraj} -xyz {outname} -box1 {L:.10f} -rdf 0.02 -1")
print("Computing q_tetr")
with open("rcut.dat","w") as f: f.write("1.8\n")
os.system(f"{mdtraj} -xyz {outname} -box1 {L:.10f} -rcut rcut.dat -edq -cn")
print("Computing q_oct")
with open("rcut.dat","w") as f: f.write("1.8\n")
os.system(f"{mdtraj} -xyz {outname} -box1 {L:.10f} -rcut rcut.dat -oct -cn")

step,pidx,ptype,cn,q_tetr = np.loadtxt("ed_q.dat", unpack=True)
sel=(cn==4)
print(f"Average coordination number during calculation of q_tetr: {cn.mean():.3f}")
np.savetxt("Q_TETR.dat",np.array([q_tetr[sel]]).T)

step,pidx,ptype,q_oct,rsrl = np.loadtxt("q_oct.dat", unpack=True)
np.savetxt("Q_OCT.dat",np.array([q_oct]).T)

os.system("rm ed_q.dat q_oct.dat")

plt.hist(q_tetr, alpha=0.8, label="tetr. %.5f"%q_tetr.mean())
plt.hist(q_oct, alpha=0.8, label="oct. %.5f"%q_oct.mean())
plt.xlabel("q")
plt.ylabel("counts")
plt.legend()
plt.savefig("q.png")
plt.show()

#os.system("rm coordnum.ave labels.dat  ndens.dat  pos.xyz ed_q.dat log q_oct.dat q_oct_time.ave rcut.dat")
