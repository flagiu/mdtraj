#!/bin/python3
import sys
import numpy as np
import joblib

q4dot_th=0.65  # threshold for crystallinity
class_dict={"SC":0,"A7":1} # classes assigned by the classifier
type_lmp_dict={"liq":1,"SC":2,"A7":3,"outlier":4} # type assigned to LAMMPS output
prob_th=0.8 # don't trust the classification if the max probability is < prob_th

if len(sys.argv)!=6:
    print("Input: <classifier.pickle> <file containing q4,...,q8 in columns> <file containing q4dot> <lammps trajectory> <output Nsc Na7 Noutlier file>")
    sys.exit(1)

classifier_str=sys.argv[1]
q_file=sys.argv[2]
q4dot_file=sys.argv[3]
lmpdump=sys.argv[4]
outfile=sys.argv[5]

clf = joblib.load(classifier_str)
print(" classes:",clf.classes_, file=sys.stderr)
print(" mapping:",class_dict, file=sys.stderr)
Xdata = np.loadtxt(q_file)
print(" data shape (N*timesteps,n_features):",Xdata.shape, file=sys.stderr)

y_pred = clf.predict(Xdata)
yprob_pred = clf.predict_proba(Xdata) # nsamples x nclasses
yprob_pred_max = np.max(yprob_pred, axis=1)

q4dot = np.loadtxt(q4dot_file)

isXtal=(q4dot>=q4dot_th)
confident=(yprob_pred_max>=prob_th)
isOutlier=isXtal&(~confident)
isSC=isXtal&confident&(y_pred==class_dict["SC"])
isA7=isXtal&confident&(y_pred==class_dict["A7"])

print("Total classification counts:", file=sys.stderr)
print(" %d liq"%sum(~isXtal), file=sys.stderr)
print(" %d xtal"%sum(isXtal), file=sys.stderr)
print("   %d SC"%sum(isSC), file=sys.stderr)
print("   %d A7"%sum(isA7), file=sys.stderr)
print("   %d outliers"%sum(isOutlier), file=sys.stderr)

types=np.ones(len(y_pred))*type_lmp_dict["liq"]
types[isOutlier] = type_lmp_dict["outlier"]
types[isSC] = type_lmp_dict["SC"]
types[isA7] = type_lmp_dict["A7"]

ofile=open(outfile,"w")
print("#timestep, N_SC, N_A7, N_outlier", file=ofile)

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
            Nsc=0; Na7=0; Noutlier=0
            iline+=1
            for i in range(N):
                line=lines[iline].strip('\n')
                words=line.split()
                typ=types[step*N+i]
                if typ==type_lmp_dict["SC"]:
                    Nsc+=1
                    typ_str="SC"
                elif typ==type_lmp_dict["A7"]:
                    Na7+=1
                    typ_str="A7"
                elif typ==type_lmp_dict["outlier"]:
                    Noutlier+=1
                    typ_str="Outlier"
                else:
                    typ_str="Liquid"
                print("%s %s %s %s"%(typ_str,words[1],words[2],words[3]))
                iline+=1
            step+=1
            print("%d %d %d %d"%(timestep,Nsc,Na7,Noutlier), file=ofile)
        else:
            print(line)
            iline+=1

ofile.close()
