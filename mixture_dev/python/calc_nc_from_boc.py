#!/usr/bin/env python3
import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Computes the number of crystalline particles per frame, given a BOC file.',
                    epilog='End of the summary.'
)
parser.add_argument('--file',  type=argparse.FileType('r'), required=True,
                     help="Input file (columns: time,particle_index,ql_dot)."
)
parser.add_argument('--threshold', type=float,
                     default=0.65, required=False,
                     help="If ql_dot > threshold, the particle is defined as crystalline. [default: %(default)s]"
)
parser.add_argument('--n', type=int,
                     default=0, required=False,
                     help="Output the index of crystalline particles in the n.th frame (non negative) in a separate file. [default: %(default)s]"
)

args = parser.parse_args()
if args.n<0:
    print("n can not be negative.")
    sys.exit(1)

t,idx,qldot = np.loadtxt(args.file, unpack=True)

t_u = np.unique(t)
nc = np.empty(len(t_u), dtype=int)
for i in range(len(nc)):
    sel = (t==t_u[i])
    xtalline = (qldot[sel]>args.threshold)
    nc[i] = np.sum( xtalline )

np.savetxt("nc.dat", np.array([t_u,nc]).T, header="#Timestep,Nc", fmt="%f %d")
print(" N. of xtalline particles saved into nc.dat")
#---------------------------------------------------------------------------------------------#
sel = (t==t_u[args.n])
xtalline = (qldot[sel]>args.threshold)
xtalline_idx = idx[sel][xtalline]
np.savetxt("xtalline.%d.dat"%args.n, xtalline_idx.reshape(1,-1), fmt="%d")
print(" Indices of xtalline particles for %d-th frame saved into xtalline.%d.dat"%(args.n,args.n))
