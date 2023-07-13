#!/bin/env/python3

import numpy as np

Nx=5
N = Nx*Nx
a = 3 # angstrom
eps = 0.1*a
nframes = 6
type = 1
L = (Nx+1)*a * np.ones(3)
print(*L)

# reticolo quadrato + rumore di scala epsilon
# xyz format
with open("test.xyz", "w") as f:
	for t in range(nframes):
		f.write(f"{N}\n")
		f.write(f"Atoms. Timestep: {t:d}\n")
		for i in range(Nx):
			for j in range(Nx):
				x = i*a + np.random.randint(100)/100 *eps
				y = j*a + np.random.randint(100)/100 *eps
				z = np.random.randint(100)/100 *a
				f.write(f"{type} {x} {y} {z}\n")
