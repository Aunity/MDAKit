#!/usr/bin/env python
'''
@Author ymh
@Email  maohuay@hotmail.com
@Date   2020-05-21 14:19:20
@Web    https://github.com/Aunity
'''

import os
import sys
import warnings
warnings.filterwarnings("ignore")
import MDAnalysis as mda

def cat_traj(trajfile:str, topfile:str) -> None:
    u = mda.Universe(os.path.abspath(topfile), os.path.abspath(trajfile))
    ts = u.trajectory
    print("Number of frames: %d"%ts.n_frames)
    print("Time step: %.2f (%s)"%(ts.dt, ts.units['time']))
    print("Total time: %.2f (%s)"%(ts.totaltime,ts.units['time']))
    print("Number of atoms: %d"%ts.n_atoms)
    print("Cell information:")
    print("                ", ts.ts.dimensions)

def main():
    if len(sys.argv[1:]) != 2:
        print('Usage:python %s <xtcf> <topf>' % sys.argv[0])
        sys.exit(0)
    xtcf, topf = sys.argv[1:]
    cat_traj(xtcf, topf)

if __name__ == '__main__':
    main()
