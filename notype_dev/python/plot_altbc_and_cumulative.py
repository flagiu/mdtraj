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
                     help="Input file with average trajectory (1st block: radial bins; 2nd block: 2d matrix of ALTBC average). [default: %(default)s]"
)

outpng="altbc_and_cumulative.png"
outpdf="altbc_and_cumulative.pdf"

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
		block=1 # separation empty line
		else:
			r.append(float(words[0]))
            if words[0]=='':
	i+=1
r = np.array(r)
nbins=len(r)

altbc = np.loadtxt(args.inavg.name, skiprows=nbins+2)
assert altbc.shape[0]==nbins
assert altbc.shape[1]==nbins

dr = r[1]-r[0]
x,y = np.meshgrid(r,r)
z = altbc
my_cmap='jet' #'rainbow' #'magma'

fig = plt.figure(dpi=300, figsize=(8,8))
gs = fig.add_gridspec(
	2, 2,
	width_ratios=(4, 1), height_ratios=(1, 4), # ratio btw marginal and main axes
	left=0.1, right=0.9, bottom=0.1, top=0.9,
    wspace=0.05, hspace=0.05
)
ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
cax = fig.add_subplot(gs[1, 1])

ax_histx.tick_params(axis="x", labelbottom=False)

ax.set_xlabel(r"$r_1$ [$\AA$]")
ax.set_ylabel(r"$r_2$ [$\AA$]")
zz = z / (4*np.pi*x**2*dr * 4*np.pi*y**2*dr)
zz/=np.max(zz)

c = ax.pcolormesh(x, y, zz, cmap=my_cmap, vmin=0.0, vmax=1.0)
ax.axis([x.min(), x.max(), y.min(), y.max()])
ax.set_aspect('equal')#, adjustable='box')
cb = fig.colorbar(c, ax=cax, location='left')
cb.ax.yaxis.set_ticks_position('right')
cb.ax.yaxis.set_label_position('right')
fig.suptitle(r"ALTBC($r_1,r_2$) @ %.1f°"%angle_th)
cax.axis('off')

cumulative_x = np.sum(z, axis=1) * dr / (4*np.pi* r**2 * dr)
cumulative_x /= max(cumulative_x)
ax_histx.plot(r, cumulative_x, '.-', color='k', label = "Cumulative (a.u.)")
ax_histx.legend()

fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()
#subprocess.call(f"xdg-open {outpng}", shell=True)
