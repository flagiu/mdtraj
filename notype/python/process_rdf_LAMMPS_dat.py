import numpy as np
import matplotlib.pyplot as plt
import sys

fname="rdf.dat"
fskip1 = 0.0 #skip the first 0.1 fraction of frames
fskip2 = 0.0 #skip the last 0.1 fraction of frames

if len(sys.argv) > 4:
	print("Usage: %s [RDF file = %s] [fskip_from_beginning = %.1f] [fskip_from_end = %.1f]"%(sys.argv[0],fname,fskip1,fskip2))
	print("Separate the RDF file produced by LAMMPS according to the timestep.")
	print("Also separates each file into a g(r) file and a coordination number file.")
	print("Note: works only if nbins and time intervals are constant in time.")
	print("")
	sys.exit(1)

if len(sys.argv) >= 2:
	fname=sys.argv[1]
	if len(sys.argv) >= 3:
		fskip1 = float(sys.argv[2])
		if len(sys.argv) == 4:
			fskip2 = float(sys.argv[3])

def myparse(input, types):
	# input: string; types: tuple of python types
	return [ t(s) for t,s in zip(types, input.split()) ]

print(f"File: {fname}")
print(f"fskip1: {fskip1:.1f}")
print(f"fskip2: {fskip2:.1f}")
with open(fname, "r") as f:
	lines = f.readlines()
	line=lines[3]
	(t, nbins) = myparse( line, (int,int) )
	ntimes_total = (len(lines)-3)//(nbins+1)
	nskip1 = int(fskip1*ntimes_total)
	nskip2 = int(fskip2*ntimes_total)
	ntimes = ntimes_total-nskip1-nskip2
	print(f"Choosen frames {nskip1+1} to {nskip1+ntimes} ({ntimes} out of {ntimes_total})")
	R = np.empty(nbins, dtype=np.float32 )		# radial distance
	RDF = np.empty(nbins, dtype=np.float32 )	# g(r)
	COORD = np.empty(nbins, dtype=np.float32 )	# integral from 0 to r of g(r)
	rdfs = []
	coords = []
	for ti in range(nskip1,nskip1+ntimes):
		irow = 3 + ti*(nbins+1)
		(t, nbins) = myparse( lines[irow], (int,int) )
		for i in range(nbins):
			irow += 1
			line = lines[irow]
			(ii, R[i], RDF[i], COORD[i]) = myparse( line, (int, float, float, float) )
#		np.savetxt(f"rdf_{t}.dat", np.array([R,RDF,COORD]).T)
		rdfs.append( RDF.copy() )
		coords.append( COORD.copy() )
	np.savetxt(f"rdf.traj", np.array([R,*rdfs]).T, header="# r, g(r,t1), g(r,t2),... at different timesteps")
	np.savetxt(f"rdf.integral.traj", np.array([R,*coords]).T, header="# r, integral of g(r,t1), g(r,t2),... at different timesteps")
