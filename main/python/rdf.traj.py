import sys
import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'


intraj="rdf.traj"
inavg="rdf.ave"
dt = 0.002 #ps
outpng="rdf.png"

if len(sys.argv)>4:
	print("Usage: %s [intraj] [inavg] [dt]\n - intraj: \t input file for trajectories (columns: t,x1,x2,...). Default: %s\n - inavg: \t input file for average (first 2 columns: t,RDF). Default: %s\n - dt: \t timestep (ps) (USELESS). Default: %1e\n"%(sys.argv[0],intraj,inavg,dt))
	sys.exit(1)

if len(sys.argv) >= 2:
	intraj = sys.argv[1]
	if len(sys.argv) >= 3:
		inavg = sys.argv[2]
		if len(sys.argv) >= 4:
			dt=float(sys.argv[3])

Xt = np.loadtxt(intraj, unpack=False)
ntraj = Xt.shape[1]-1

Xa = np.loadtxt(inavg, unpack=False)
r = Xa[:,0]
y = Xa[:,1]

assert Xt.shape[0]==Xa.shape[0]
assert (Xt[:,0]==Xa[:,0]).all()
assert Xa.shape[1]>=2

fig, ax = plt.subplots()
ax.set_xlabel(r"$r$ / $\AA$")
ax.set_ylabel(r"$g(r)$")
ax.tick_params(which='both', direction='in')
for i in range(ntraj):
	ax.plot( r, Xt[:,i+1], 'k', alpha=0.1 )
ax.plot(r, y, 'r')
ax.grid(axis='both', which='major')
plt.tight_layout()
fig.savefig(outpng)
print("Figure saved on %s\n"%(outpng))
#plt.show()

