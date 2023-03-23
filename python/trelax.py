import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

outpng="trelax.png"
m=121.750
amu2kg=1.66054
parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Evaluates the relaxation time of the system at the time at which MSD(t)=R^2, being R the position of the first peak of the radial distribution function.',
                    epilog='End of the summary.'
)
parser.add_argument('pathsFile', type=str,
                     help="File containing a column of paths ot folders, each containing RDF and MSD."
)
parser.add_argument('--rdf', type=str,
                     default="RDF.ave", required=False,
                     help="Pattern for the RDF file name."
)
parser.add_argument('--msd',  type=str,
                     default="msd.ave", required=False,
                     help="Pattern for the MSD file name"
)
parser.add_argument('--dt',  type=float, default=0.002, required=False,
                     help="Integration time step (picoseconds)"
)
parser.add_argument('--density',  type=float, default=6.49, required=False,
                     help="Numerical density (g/cm3)."
)
parser.add_argument('-p','--plot',  action='store_true',
                     help="Plot g(r), S(q) and MSD(t)."
)

args = parser.parse_args()
rho = args.density /(m*amu2kg) # 1/angstrom^3
if args.plot:
	fig, axes = plt.subplots(1,3, figsize=(16,6))
	axes[0].set_xlabel(r"$r$ ($\AA$)")
	axes[0].set_ylabel(r"$g(r)$")
	axes[1].set_xlabel(r"$q$ ($\AA^{-1}$)")
	axes[1].set_ylabel(r"$S(q)$")
	axes[2].set_xlabel(r"$t$ (ps)")
	axes[2].set_ylabel(r"$\langle\, \Delta r^2(t)\, \rangle$ ($\AA^2$)")
print(f"# path | g(r) peak (Ang) | t_relax (ps) | upper bound? 1:0 ")
with open(args.pathsFile, "r") as f:
	paths = f.readlines()
	for line in paths:
		path=line[:-1] #remove '\n'
		infile=path+'/'+args.rdf
		r,g,g_ = np.loadtxt(infile, unpack=True, usecols=(0,1,2))
		r0 = r[ g.argmax() ]
		#----------- structure factor ---------------#
		dr = r[1]-r[0]
		dq = 2*np.pi / (2*r[-1])
		qmax = 15 #2*np.pi / r[sum(g==0)-1]
		q = np.linspace(dq, qmax, int(qmax/dq) )
		Sq = np.empty(len(q))
		for i,Q in enumerate(q):
			Sq[i] = 1.0 + 4*np.pi*rho/Q * np.sum( dr*r*np.sin(Q*r)*(g-1.0) )
		#--------------------------------------------#
		infile=path+'/'+args.msd
		t,msd,msd_ = np.loadtxt(infile, unpack=True, usecols=(0,1,2))
		t *= args.dt
		t0 = t[ ((msd-r0*r0)**2).argmin() ]
		tag= ( 1 if t0==t[0] else 0 )
		print(f"{path} {r0:.3g} {t0:.3g} {tag}")
		idx = (t<=10*t0)
		if args.plot:
			axes[0].plot(r,g)
			axes[0].axvline(r0, color='k', alpha=0.3)
			axes[1].plot(q,Sq)
			axes[2].plot(t[idx],msd[idx])
			axes[2].axhline(r0*r0, color='k', alpha=0.3)
			axes[2].axvline(t0,    color='k', alpha=0.3)
if args.plot:	
	for ax in axes:
		ax.tick_params(which='both', direction='in')
		ax.grid(axis='both', which='major')
	plt.tight_layout()
	fig.savefig(outpng)
	print("Figure saved on %s\n"%(outpng))
	plt.show()

