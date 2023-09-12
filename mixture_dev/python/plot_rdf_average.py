#!/bin/python3

import sys
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'

def int2types(t, nTypes):
    for x in range(nTypes):
      low = int(x * nTypes - x*(x-1)/2)
      high = int( (x+1) * nTypes - x*(x+1)/2 - 1 )
      if t>=low and t<=high:
        #print(t,x,low,high,x, x + (int(t)-low))
        return x, x + (int(t)-low)
    return None,None

def linMap(x,a,b,c,d):
   y=(x-a)/(b-a)*(d-c)+c
   return y

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Plots the RDF for each pair of types.',
                    epilog='End of the summary.'
)
parser.add_argument('--inavg',  type=argparse.FileType('r'),
                     default="rdf.ave", required=False,
                     help="Input file with average RDF (r | g_00 g_01 ... | errors). [default: %(default)s]"
)

args = parser.parse_args()

outpng="rdf.png"
outpdf="rdf.pdf"

X = np.loadtxt(args.inavg.name)
r = X[:,0]
npairs = int( (X.shape[1]-1)/2 )
ntypes = int(np.floor( (2*npairs)**0.5 ))

fig, ax = plt.subplots(dpi=300)
ax.set_xlabel(r"$r$ [$\AA$]")
ax.set_ylabel(r"$g(r)$")
for i in range(npairs):
    g = X[:,1+i]
    g_ = X[:,1+i+npairs]
    ti,tj = int2types(i, ntypes)
    lab = "%d-%d"%( ti,tj )

    red = linMap(ti, 0,ntypes-1, 1,0)
    blue = linMap(tj, 0,ntypes-1, 0,1)
    green = 0.2
    ax.errorbar(r,g,g_, label=lab, color=(red,green,blue,0.7))
ax.legend()
ax.tick_params(which='both', direction='in')
ax.grid(axis='both', which='major')
plt.tight_layout()

fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()
#subprocess.call(f"xdg-open {outpng}", shell=True)
