#!/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import joblib

ls=np.arange(4,22+1) #np.array([4,5,6,7,8],dtype=int)
class_dict={"SC":0,"A7":1}
classifier_str="classifier.pkl"
q4dot_th=0.65
typeSClammps=2
typeA7lammps=3

if len(sys.argv)!=6:
    print("Input: <classifier.pickle> <file containing q4,...,q22 columns> <file containing q4dot> <lammps trajectory> <output data file>")
    sys.exit(1)

classifier_str=sys.argv[1]
q_file=sys.argv[2]
q4dot_file=sys.argv[3]
lmpdump=sys.argv[4]
outfile=sys.argv[5]

clf = joblib.load(classifier_str)

Xdata = np.loadtxt(q_file)
# shape: (N*len_unique_timesteps,len(ls))
print(" data shape:",Xdata.shape, file=sys.stderr)

y_pred = clf.predict(Xdata)

for label in class_dict.keys():
    cl=class_dict[label]
    n=sum(y_pred==cl)
    print("Predicted %d %s"%(n,label), file=sys.stderr)

q4dot = np.loadtxt(q4dot_file)

isXtal=(q4dot>=q4dot_th)

types=np.ones(len(y_pred)) # liquid type: 1
types[isXtal&(y_pred==class_dict["SC"])] = typeSClammps
types[isXtal&(y_pred==class_dict["A7"])] = typeA7lammps

ofile=open(outfile,"w")
print("#timestep, N_SC, N_A7", file=ofile)

with open(lmpdump,"r") as f:
	lines=f.readlines()
	iline=0
	N=-1
	step=0
	while iline<len(lines):
		line=lines[iline].strip('\n')
		words=line.split()
		if line=="ITEM: TIMESTEP":
			timestep=int(lines[iline+1])
		if line=="ITEM: NUMBER OF ATOMS":
			N=int(lines[iline+1])
		
		if line=="ITEM: ATOMS type x y z":
			print(line)
			Nsc=0; Na7=0
			iline+=1
			for i in range(N):
				line=lines[iline].strip('\n')
				words=line.split()
				typ=types[step*N+i]
				if typ==typeSClammps: Nsc+=1
				elif typ==typeA7lammps: Na7+=1
				print("%d %s %s %s"%(typ,words[1],words[2],words[3]))
				iline+=1
			step+=1
			print("%d %d %d"%(timestep,Nsc,Na7), file=ofile)
		else:
			print(line)
			iline+=1

ofile.close()
