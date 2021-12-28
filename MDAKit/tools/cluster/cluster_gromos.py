#!/usr/bin/env python
'''
@Author ymh
@Email  maohuay@hotmail.com
@Date   2019-09-29 23:53:33
@Web    https://github.com/Aunity
'''
import os
import sys
import argparse
import numpy as np
import mdtraj as md
from msmbuilder import libdistance

def maxmatrix(mat):
    ndx0 = mat.argmax(axis=1)
    maxv = mat.max(axis=1)
    y = maxv.argmax()
    v = maxv.max()
    x = ndx0[y]
    return x,y,v

def gromos():
    trajs
    rmsdmax = -9999
    centers = {}
    for traj0 in trajs:
        for traj1 in trajs:
            rmsd = libdistance.cdist(traj0, traj1, metric="rmsd")
            x,y,v = maxmatrix(rmsd)
            if v>=rmsdmax:
                pass

def main():
    if len(sys.argv[1:]) != {}:
        print('Usage:python %s <> <>'%sys.argv[0])
        sys.exit(0)
    xtcfp, top, cutoff, stride
if __name__ == '__main__':
    main()
