import numpy as np
import matplotlib.pyplot as plt
import sys

dt = 0.002 #picoseconds
n_plot = 1

if len(sys.argv)<2 or len(sys.argv)>4:
    print("\nUsage: %s <input file> [number of plots = 1] [time step = 1e-2ps]\n"%sys.argv[0])
    print(" The file should have 3 columns: time index, particle index, value.")
    print(" Comments must be marked by '#'.\n")
    sys.exit(1)

if len(sys.argv)>2:
	n_plot = int(sys.argv[2])
	if len(sys.argv)>3:
		dt = float(sys.argv[3])

ts = [] # list of time steps

with open(sys.argv[1], "r") as f:
	lines = f.readlines()
	for line in lines:
		if line[0]=='#':
			continue
		ts.append( int( line.split()[0] ) )

ndata = len(ts)
ts = np.unique( np.array(ts, dtype=np.int32) )
ntimes = len(ts)

print(f"Number of data lines: {ndata}")
print(f"Number of distinct timesteps: {ntimes}")

ti = 0
xs = [ [], ] # list of ( list of values per particle ) per time steps
with open(sys.argv[1], "r") as f:
	lines = f.readlines()
	for line in lines:
		if line[0]=='#':
			continue
		a,b,c = line.split()
		current_timestep = int(a)
		if current_timestep != ts[ti]:
			ti+=1
			xs.append([])
		xs[ti].append( float(c) )


for i in range(0,ntimes, ntimes//n_plot):
	print(i)
	t = ts[i]*dt
	x = xs[i]
	plt.xlabel("x")
	plt.ylabel("density")
	plt.hist(x, label=f"{t:.1f} ps", density=True, alpha=0.5)

plt.legend()
plt.show()
