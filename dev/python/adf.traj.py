import sys
import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'


intraj="adf.traj"
inavg="adf.ave"
dt = 0.002 #ps
outpng="adf.png"
outpdf="adf.pdf"
x_tolerance=0.01 # error tolerance on the bins (fraction)

if len(sys.argv)>4:
	print("Usage: %s [intraj] [inavg] [dt]\n - intraj: \t input file for trajectories (columns: t,x1,x2,...). Default: %s\n - inavg: \t input file for average (first 2 columns: t,ADF). Default: %s\n - dt: \t timestep (ps) (USELESS). Default: %1e\n"%(sys.argv[0],intraj,inavg,dt))
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

assert Xt.shape[0]==Xa.shape[0] # space binning must have same length
x_isnot_equal = abs( (Xt[:,0]-Xa[:,0])/Xa[:,0] )>=x_tolerance
try:
	assert not x_isnot_equal.any() # space binning must be the same
except AssertionError:
	print("[ ERROR: x values do not match! ]")
	print(" File",intraj)
	print(Xt[x_isnot_equal, 0])
	
	print(" File",inavg)
	print(Xa[x_isnot_equal, 0])
	sys.exit(1)

assert Xa.shape[1]>=2 # must have at least x,y

fig, ax = plt.subplots(dpi=300)
ax.set_xlabel(r"$\cos\theta$")
ax.set_ylabel(r"ADF")
ax.tick_params(which='both', direction='in')
for i in range(ntraj):
	ax.plot( r, Xt[:,i+1], 'k', alpha=0.01)
ax.plot(r, y, 'r')
ax.grid(axis='both', which='major')
plt.tight_layout()
fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()

