#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

argv =sys.argv
if len(argv) < 2:
    print("Please supply a file with data to plot!")
    quit()
fn=argv[1]
plt.title('Dati')
plt.xlabel('x')
plt.ylabel('y')
x, y = np.loadtxt(fn, comments=['#'], usecols=(0,1), unpack=True)
plt.plot(x, y, 'x-',label='DAti')
#plt.savefig('traiettoria.png')
plt.show()
