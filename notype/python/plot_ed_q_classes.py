#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'
plt.rcParams['figure.dpi'] = 200

outpng="ed_q_classes.png"
outpdf="ed_q_classes.pdf"
class2label = {-1:'a', -2:'b', -3:'c', -4:'d', -5:'e', -6:'f', -7:'g'}
class_labels = ['a','b','c','d','e','f','g']
class_explanations =  [
    '3-fold planar defective octahedral',
    '3-fold pyramidal defective octahedral',
    '4-fold defective octahedral',
    '4-fold planar defective octahedral',
    'tetrahedral',
    '5-fold defective octahedral',
    'octahedral'
]
def class2idx(cla: int): # map -1,-2,... to 0,1,...
    return int(abs(cla)-1)

parser = argparse.ArgumentParser(
                    prog = sys.argv[0],
                    description = 'Time series of number of each local structure motifs according to Eddington-Debenedettin "q" bond order parameter.',
                    epilog='End of the summary.\n'
)
parser.add_argument('--indat',  type=argparse.FileType('r'), required=False, default="ed_q_classes.dat",
                     help="Input file: timestep | particle index | E-D class (from -1 to -7). [default: %(default)s]"
)

parser.add_argument('--fskip0', type=float,
                     default=0.0, required=False,
                     help="Fraction of data to be skipped from the start. [default: %(default)s]"
)
parser.add_argument('--fskip1', type=float,
                     default=0.0, required=False,
                     help="Fraction of data to be skipped from the end. [default: %(default)s]"
)
parser.add_argument('--dt',  type=float, default=0.002, required=False,
                     help="Integration time step (picoseconds)"
)

args = parser.parse_args()

header=args.indat.readline()
rcut1 = float( header.split("# cutoff =")[1].split(',')[0] )

x = np.loadtxt(args.indat)
# filter by time
t = np.unique(x[:,0])
indices = x[:,1]
n = len(t) # number of frames
n0 = int(args.fskip0*n)
n1 = int(args.fskip1*n)
assert n-n0-n1 > 0 # cannot skip all lines!
t0 = t[n0]
t1 = t[n-1-n1]
selection = (x[:,0]>=t0) & (x[:,0]<=t1)
x = x[ selection ]

t = np.unique(x[:,0])
data = x[:,2]

classes_vs_time=np.zeros( (len(class2label), len(t)) )
for i in range(len(t)):
    classes, counts = np.unique( data[ x[:,0]==t[i] ], return_counts=True)
    for cla,count in zip(classes,counts):
        classes_vs_time[ class2idx(cla), i ] = count

fig, ax = plt.subplots(figsize=(10,6))
ax.set_xlabel("t [ps]")
ax.set_ylabel("Counts")
title=r"$r_{cut}=%.2f$ $\AA$"%rcut1
ax.set_title(title)
for i in range(len(class2label)):
    ax.plot(t*args.dt, classes_vs_time[i], '-', label=class_labels[i], alpha=0.7)
ax.legend(bbox_to_anchor=(1.01,1.01))
#ax.set_yscale("log")
ax.tick_params(axis='y',which='both', direction='in')
#ax.grid(axis='both', which='major')
plt.tight_layout()

plt.savefig(outpng)
plt.savefig(outpdf)
print("%s: Figure saved on %s , %s\n"%(sys.argv[0],outpng, outpdf))
#plt.show()
