#!/bin/env/python3

import numpy as np

# reticolo cubico (ortorombico) + rumore di scala epsilon
Nx=Ny=Nz=10
a=b=c = 3 # lattice spacing
eps = 0.1*a
nframes = 50
def choose_type(i,j,k):
	#return np.random.randint(3) # 3 types 0,1,2, uniformly chosen
	return np.random.randint(4)%3 # 3 types 0,1,2, but 0 twice more probable than 1,2
#-----------------------------------------------------#

N = Nx*Ny*Nz
L = np.array( [(Nx+1)*a, (Ny+1)*b, (Nz+1)*c] )
print(*L)

# xyz format
with open("traj.xyz", "w") as f:
	for t in range(nframes):
		f.write(f"{N}\n")
		f.write(f"Atoms. Timestep: {t:d}\n")
		for i in range(Nx):
			for j in range(Nx):
				for k in range(Nx):
					type = choose_type(i,j,k)
					x = i*a + np.random.randint(100)/100 *eps
					y = j*b + np.random.randint(100)/100 *eps
					z = k*c + np.random.randint(100)/100 *eps
					#z = np.random.randint(100)/100 *a
					f.write(f"{type} {x} {y} {z}\n")
