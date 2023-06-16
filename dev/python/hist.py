import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv)!=2:
    print("\nUsage: %s <1-column file>\n"%sys.argv[0])

x = np.loadtxt(sys.argv[1])

plt.xlabel("x");
plt.ylabel("counts")
plt.hist(x)
plt.show()
