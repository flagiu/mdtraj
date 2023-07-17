#!/bin/env/python3

import numpy as np

# reticolo cubico + rumore di scala epsilon
Nx=10
N = Nx*Nx*Nx
a = 3 # lattice spacing
eps = 0.1*a
nframes = 50
L = (Nx+1)*a * np.ones(3)
print(*L)

# 2 tipi
nType0 = int(0.8*N)

# xyz format
n=0
with open("test.xyz", "w") as f:
	for t in range(nframes):
		f.write(f"{N}\n")
		f.write(f"Atoms. Timestep: {t:d}\n")
		for i in range(Nx):
			for j in range(Nx):
				for k in range(Nx):
					type = 0 if n < nType0 else 1 # type must be ordered
					x = i*a + np.random.randint(100)/100 *eps
					y = j*a + np.random.randint(100)/100 *eps
					z = k*a + np.random.randint(100)/100 *eps
					#z = np.random.randint(100)/100 *a
					f.write(f"{type} {x} {y} {z}\n")
					n+=1
