import numpy as np
import matplotlib.pyplot as plt
import sys

dt = 0.002 #picoseconds
timestep_to_plot = 0

if len(sys.argv)<2 or len(sys.argv)>5:
    print("\nUsage: %s <input file 1 > <input file 2> [index of timestep to be plotted = %d] [time step = %fps]\n"%(sys.argv[0],timestep_to_plot,dt))
    print(" The files should have 3 columns: time index, particle index, value.")
    print(" Comments must be marked by '#'")
    print(" The two files must derive from the same trajectory.\n")
    sys.exit(1)

file1 = sys.argv[1]
file2 = sys.argv[2]
if len(sys.argv)>3:
	timestep_to_plot = int(sys.argv[3])
	if len(sys.argv)>4:
		dt = float(sys.argv[4])

#--------- READ TIME STEPS FROM FILE 1 ----------#
ts = [] # list of time steps

with open(file1, "r") as f:
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
#--------------------------------------------------#

#--------- READ FILE 1 ----------#
ti = 0
xs = [ [], ] # list of ( list of values per particle ) per time steps
with open(file1, "r") as f:
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
#--------------------------------#

#--------- READ FILE 2 ----------#
ti = 0
ys = [ [], ] # list of ( list of values per particle ) per time steps
with open(file2, "r") as f:
	lines = f.readlines()
	for line in lines:
		if line[0]=='#':
			continue
		a,b,c = line.split()
		current_timestep = int(a)
		if current_timestep != ts[ti]:
			ti+=1
			ys.append([])
		ys[ti].append( float(c) )
#--------------------------------#

i = timestep_to_plot
t = ts[i]*dt
x = xs[i]
y = ys[i]
plt.xlabel("x")
plt.ylabel("y")
plt.scatter(x, y, label=f"{t:.1f} ps", alpha=0.7)

plt.legend()
plt.show()
