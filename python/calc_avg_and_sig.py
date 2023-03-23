#!/usr/bin/env python3
import sys,os
import numpy as np
if (len(sys.argv)>1):
    ob=sys.argv[1]
else:
    print("Usage:\n calc_avg_and_sig.py <data file (two columns)> <fraction of data to discard (default:0.9)>")
    quit()
if len(sys.argv) > 2:
    fract=float(sys.argv[2].strip('\n'))
else:
    fract=0.9
with open(ob) as f:
    lob=f.readlines()
#print('lphi=', lphi)
nob=len(lob)
snob=int(fract*nob)
#print ('press=', press, 'nphi=', nphi, ' snphi=', snphi)
#print('snphi=', snphi, ' fine=', nphi)
if snob==nob:
    snob=nob-1
subl=lob[snob:nob]
#print('subl=', subl)
sumob=0.0
sumob2=0.0
cc=0.0
for ll in subl:
    obval=float(ll.strip('\n').split(' ')[1])
    #print('phi=',phi)
    sumob+=obval
    sumob2+=obval**2
    cc += 1.0
print(sumob/float(cc),np.sqrt(sumob2/float(cc)-(sumob/float(cc))**2))
