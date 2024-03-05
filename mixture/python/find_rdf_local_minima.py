#!/bin/python3

import sys
import numpy as np
from scipy.signal import argrelextrema, find_peaks

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


def main(rdf_file: str, output_file: str, w: str, n: int, delta: float, debug: bool):
    X = np.loadtxt(rdf_file)
    p = (X.shape[1]-1)//2
    x = X[:,0]

    rcut=np.empty((3,p))
    if debug: import matplotlib.pyplot as plt

    for ip in range(p):
        y = X[:,1+ip]
        
        if debug: plt.plot(x,y,label='data',color='black')
        
        ok=False
        while not ok:
            print("Finding g(r) minima: "
                  "delta=%.2f resolution: n=%d smoothing points"%(delta,n)
                  )
            xs = smooth(x,n,w)
            ys = smooth(y,n,w)
            
            if debug: plt.plot(xs,ys,label="pair %d , n=%d"%(ip+1,n), alpha=0.5)
            # find local minima and maxima
            imax,_ = find_peaks( ys, prominence=0.01*max(ys))
            imin,_ = find_peaks(-ys, prominence=0.01*max(ys))
            xmin = xs[imin]
            xmax = xs[imax]
            ymin = ys[imin]
            ymax = ys[imax]
            # criterion: must resolve at least 2 minima
            ok=len(xmin)>=2
            if ok: ok=ok&((xmin[1]-xmin[0])>delta)
            if len(xmin)>2: ok=ok&((xmin[2]-xmin[1])>delta)
            # this may be False if g(r) is noisy --> increase smoothing
            n+=2

        print("Finding g(r) minima: DONE type pair [ %d / %d ]\n"%(ip+1,p))
        rcut[0,ip] = xmin[0]
        rcut[1,ip] = xmin[1]
        rcut[2,ip] = (xmin[2] if len(xmin)>2 else xmin[1])

    if debug:
        plt.scatter(xmin,ymin, marker='x',color='red',alpha=0.7, label='minima')
        plt.scatter(xmax,ymax, marker='o',color='green',alpha=0.7, label='maxima')
        plt.legend()
        plt.show()

    # print the first 3 local minima in column
    np.savetxt(output_file, rcut, fmt="%.3f")

    return


if __name__=='__main__':
    w = 'hanning'
    n = 10
    delta = 0.3 # angstrom
    debug = False
    if len(sys.argv)<3 or len(sys.argv)>7:
        print("Input: <g(r) file> <r_cut output file> [window=%s] [n=%d] [delta=%f] [debug=%d]"%(w,n,delta,debug))
        print("  The g(r) file must have the following 1+2*p columns,")
        print("  where p is the number of type pairs:")
        print("    r | g_00(r), g_01(r), ... | error for each g(r)")
        print()
        print("  I will smooth each g(r) using a moving the given window")
        print("  of n points, then find the position of the first 3")
        print("  local minima using `scipy.signal.argrelextrema` and")
        print("  save them as a (3 x n) matrix to the output file.")
        print()
        print("  While any two minima are found to be closer than delta")
        print("  (distance units), I will assume that this is due to noise,")
        print("  so I will repeat with a +2 larger smoothing window.")
        print()
        sys.exit(1)
    if len(sys.argv)>3: w = sys.argv[3]
    if len(sys.argv)>4: n = int(sys.argv[4])
    if len(sys.argv)>5: delta = float(sys.argv[5])
    if len(sys.argv)>6: debug = bool(sys.argv[6])
    main(sys.argv[1], sys.argv[2], w, n, delta, debug)
