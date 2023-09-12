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

fig, ax = plt.subplots(dpi=300)
ax.set_xlabel(r"$r$ [$\AA$]")
ax.set_ylabel(r"Partial $g(r)$")
for i in range(npairs):
    g = X[:,1+i]
    g_ = X[:,1+i+npairs]
    ax.errorbar(r,g,g_, label=i)
ax.legend()
plt.tight_layout()

fig.savefig(outpng)
fig.savefig(outpdf)
print("Figure saved on %s, %s\n"%(outpng, outpdf))
#plt.show()
#subprocess.call(f"xdg-open {outpng}", shell=True)
