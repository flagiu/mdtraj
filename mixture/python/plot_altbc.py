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
parser.add_argument('--normalize',  type=bool,
                     default=False, required=False,
                     help="Normalize to the maximum bin value? [default: %(default)s]"
)
parser.add_argument('--outname',  type=str,
                     default="altbc", required=False,
                     help="Prefix for the output file. [default: %(default)s]"
)

#-------------------------------------#
args = parser.parse_args()

print("Plotting ALTBC ...")

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

dr = r[1]-r[0]
x,y = np.meshgrid(r,r)
z = altbc
my_cmap='jet' #'rainbow' #'magma'

fig, ax = plt.subplots(dpi=300, figsize=(6,5)) # (6,5) is experimental
ax.set_xlabel(r"$r_1$ [$\AA$]")
ax.set_ylabel(r"$r_2$ [$\AA$]")
ax.set_title(r"ALTBC($r_1,r_2,%.1f$°)"%angle_th)
if args.normalize:
    zz = z / np.max(z)
else:
    zz = z # / (4*np.pi*x**2*dr * 4*np.pi*y**2*dr) # this is done in the altbc.hpp file
c = ax.pcolormesh(x, y, zz, cmap=my_cmap, vmin=zz.min(), vmax=zz.max())
ax.axis([x.min(), x.max(), y.min(), y.max()])
ax.set_aspect('equal', adjustable='box')
cb = fig.colorbar(c, ax=ax, shrink=0.8, aspect=20*0.8) # 0.8 is experimental
plt.tight_layout()

fig.savefig(args.outname+".pdf")
fig.savefig(args.outname+".png")
print("Figure saved on %s.png , %s.pdf\n"%(args.outname,args.outname))
#plt.show()
#subprocess.call(f"xdg-open {outpng}", shell=True)
