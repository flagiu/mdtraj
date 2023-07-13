import numpy as np
import matplotlib.pyplot as plt
import sys

fname="rdf.dat"
fskip1 = 0.0 #skip the first 0.1 fraction of frames
fskip2 = 0.0 #skip the last 0.1 fraction of frames
outfile_rdf = "rdf.ave"
outfile_rdf_integral = "rdf.integral.ave"

if len(sys.argv) > 4:
	print("Usage: %s [RDF file = %s] [fskip_from_beginning = %.1f] [fskip_from_end = %.1f]"%(sys.argv[0],fname,fskip1,fskip2))
	print("Calculates the average of the RDF and of its integral, both produced by LAMMPS.")
	print("Results are saved into %s and %s"%(outfile_rdf,outfile_rdf_integral))
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
	RDF = np.zeros(nbins, dtype=np.float32 )	# <g(r)>
	RDF2 = np.zeros(nbins, dtype=np.float32 )	# <g(r)^2>
	COORD = np.zeros(nbins, dtype=np.float32 )	# <integral from 0 to r of g(r)>
	COORD2 = np.zeros(nbins, dtype=np.float32 )	# <integral from 0 to r of g(r) )^2>
	for ti in range(nskip1,nskip1+ntimes):
		irow = 3 + ti*(nbins+1)
		(t, nbins) = myparse( lines[irow], (int,int) )
		for i in range(nbins):
			irow += 1
			line = lines[irow]
			(ii, R[i], g, c) = myparse( line, (int, float, float, float) )
			RDF[i] += g / ntimes
			RDF2[i] += g*g / ntimes
			COORD[i] += c / ntimes
			COORD2[i] += c*c / ntimes

	RDF2 = np.nan_to_num( np.sqrt((RDF2-RDF*RDF)/ntimes), nan=0.0 )
	data = np.array([R,RDF,RDF2]).T
	np.savetxt(outfile_rdf, data, header="r,<g(r)>,fluctuation")

	COORD2 = np.nan_to_num( np.sqrt((COORD2-COORD*COORD)/ntimes), nan=0.0 )
	data = np.array([R,COORD,COORD2]).T
	np.savetxt(outfile_rdf_integral, data, header="r,<integral of g(r)>,fluctuation")

