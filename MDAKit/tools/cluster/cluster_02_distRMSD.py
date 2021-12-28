#!/usr/bin/env python
'''
@Author ymh
@Email  maohuay@hotmail.com
@Date   2019-09-30 11:17:24
@Web    https://github.com/Aunity
'''
import os
import sys
import numpy as np

def hist_rmsd(rmsdf, rmin, rmax, n):
    rmsd = np.load(rmsdf)['a']
    txx = np.concatenate(rmsd)
    h,e = np.histogram(txx, bins=n, range=(rmin,rmax))
    e = [(e[i]+e[i+1])/2 for i in range(len(e)-1)]
    return h,e

def save_hist(h,e,outf):
    h0 = h*1.0/sum(h)
    data = np.array([e,h,h0]).T
    np.savetxt(outf, data, fmt="%.4e", header="RMSD count density")

def main():
    if len(sys.argv[1:]) != 4:
        print('Usage:python %s <npzfp> <rmin> <rmax> <nsize>'%sys.argv[0])
        sys.exit(0)
    npzfp, rmin, rmax, n = sys.argv[1], np.float(sys.argv[2]), np.float(sys.argv[3]), int(sys.argv[4])
    h0, e0 = [], []
    if os.path.isfile(npzfp):
        outf = os.path.split(npzfp)[1]+".prob"
        h,e = hist_rmsd(npzfp, rmin, rmax, n)
        save_hist(h,e,outf)
    else:
        outf = "rmsdmatrix.prob"
        npzfs = [os.path.join(npzfp,_) for _ in os.listdir(npzfp) if _.endswith("npz")]
        for npzf in npzfs:
            h,e = hist_rmsd(npzf, rmin, rmax, n)
            if len(e0) == 0:
                h0,e0 = h,e
            else:
                h0,e0 = h0+h,e
        save_hist(h,e,outf)

if __name__ == '__main__':
    main()
