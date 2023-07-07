#!/bin/python3

import sys
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'


parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Plots the Angular-Limited Three-Body Correlation as a 2D heatmap.',
                    epilog='End of the summary.'
)
parser.add_argument('--inavg',  type=argparse.FileType('r'),
                     default="altbc.ave", required=False,
                     help="Input file with average trajectory (1st block: radial bins; 2nd block: 2d matrix of ALTBC average; 3rd block: 2d matrix of ALTBC error). [default: %(default)s]"
)

outpng="altbc.png"
outpdf="altbc.pdf"

#-------------------------------------#
args = parser.parse_args()

lines = args.inavg.readlines()
i = 0
block=0
r = []
while block==0:
	line = lines[i]
	words = line.strip('\n').split(' ')
	if line[0]=='#':
		if i==0:
			angle_th = float(words[-2])
	else:
		assert len(words)==1
		if words[0]=='':
			block=1 # separation empty line
		else:
			r.append(float(words[0]))			
	i+=1
r = np.array(r)
nbins=len(r)

altbc = np.loadtxt(args.inavg.name, skiprows=nbins+2)
assert altbc.shape[0]==nbins
assert altbc.shape[1]==nbins

x,y = np.meshgrid(r,r)
z_min = altbc.min()
z_max = altbc.max()

fig, ax = plt.subplots(dpi=300)
ax.set_xlabel(r"$r_1$ [$\AA$]")
ax.set_ylabel(r"$r_2$ [$\AA$]")
ax.set_title(r"ALTBC @ %.1f°"%angle_th)
c = ax.pcolormesh(x, y, altbc, cmap='jet', vmin=z_min, vmax=z_max)
ax.axis([x.min(), x.max(), y.min(), y.max()])
fig.colorbar(c, ax=ax)

plt.tight_layout()
fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()
subprocess.call(f"xdg-open {outpng}", shell=True)
