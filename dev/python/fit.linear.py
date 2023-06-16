import sys
import numpy as np
from scipy.optimize import curve_fit

xcol=int(0)
ycol=int(1)
dycol=int(2)
fskip1 = 0.0
fskip2 = 0.0

if len(sys.argv)<2 or len(sys.argv)>7:
	print("Usage: %s <input file (2 columns)> [x-column index = %d] [y-column index = %d] [dy-column index = %d] [fskip from start = %.2f] [fskip from end = %.2f]\n"%(sys.argv[0],xcol,ycol,fskip1,fskip2))
	sys.exit(1)

fname = sys.argv[1]
if len(sys.argv) >= 3:
	xcol = int(sys.argv[2])
	if len(sys.argv) >= 4:
		ycol = int(sys.argv[3])
		if len(sys.argv) >= 5:
			dycol = int(sys.argv[4])
			if len(sys.argv) == 6:
				fskip1 = float(sys.argv[5])
				if len(sys.argv) == 7:
					fskip2 = float(sys.argv[6])
"""
print(f"File: {fname}")
print(f"x-column index: {xcol}")
print(f"y-column index: {ycol}")
print(f"dy-column index: {dycol}")
print(f"fskip1: {fskip1:.1f}")
print(f"fskip2: {fskip2:.1f}")
"""

x, y, dy = np.loadtxt(fname, unpack=True, usecols=(xcol,ycol,dycol))

n = len(x)
n1 = int(fskip1*n)
n2 = int(fskip2*n)

x = x[n1: n-n2]
y = y[n1: n-n2]
dy = dy[n1: n-n2]

def f(x,a=1):
	return a*x

popt,pcov = curve_fit(f, x,y, sigma=dy, p0=[1], maxfev=10000)

popt_ = np.sqrt(np.diag(pcov))

for a,a_ in zip(popt,popt_):
	print(a,a_)
