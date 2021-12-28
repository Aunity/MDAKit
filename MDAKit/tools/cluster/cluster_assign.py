#!/usr/bin/env python
'''
@Author ymh
@Email  maohuay@hotmail.com
@Date   2019-09-30 09:30:01
@Web    https://github.com/Aunity
'''
import os
import sys
import argparse
import numpy as np
import mdtraj as md
import pandas as pd
from msmbuilder import libdistance
from multiprocess import Pool

def parallel(args):
    xtcf, cf, top, sel, dt = args
    seq = assign_(xtcf, cf, top, sel, dt=1, outfname=None)
    return seq

def assign_(xtcf, cf, top, sel, dt=1, outfname=None):
    if not outfname:
        outfname = os.path.split(xtcf)[1].split(".")[0]
    top = md.load_pdb(top)
    ndx = top.top.select(sel)
    centerpdbs = md.load(cf, top=top, atom_indices=ndx)
    xtc = md.load_xtc(xtcf, top=top, atom_indices=ndx, stride=dt)
    seq,iter0 = libdistance.assign_nearest(xtc, centerpdbs, metric="rmsd")
    np.savez_compressed(outfname+'.npz', a=seq)
    return seq

def parse_args():
    parser = argparse.ArgumentParser(description="Cal RMSD matrix for trajectory(s).")
    parser.add_argument("-i", dest="inp", help="Trajectory file or floder contains trajectories, or meta pkl3 file contains infor.",required=True)
    parser.add_argument("-t", dest="top", help="Topology file for trajectory.", default=".")
    parser.add_argument("-c", help="Cluster center pdbs or xtcfile.", required=True)
    parser.add_argument("-s", dest="dt", help="every dt time",  type=int, default=1)
    parser.add_argument("-sel", help="mdtraj select str to extract atoms for cal RMSD. defatult: name CA", default="name CA")
    parser.add_argument("-o", dest="oup", help="A floder or a name to save the RMSD matrix npy file.", default="rmsd_mat")
    parser.add_argument("-T", help="Number of thread to cal RMSD. default=1", type=int, default=1)

    args = parser.parse_args()
    return os.path.abspath(args.inp), os.path.abspath(args.top), args.sel, args.oup, args.T, args.dt, os.path.abspath(args.c)

def main():
    xtcfp,top,sel,outp,T,dt,cf = parse_args()

    if os.path.isfile(xtcfp) and xtcfp.endswith(".xtc"):
        #seq = cal_rmsdmatrix(xtcfp, cf, top, sel, dt)
        assign_(xtcfp, cf, top, sel, dt=dt)
    else:
        if not os.path.exists(outp):
            os.mkdir(outp)
        cwd0 = os.getcwd()
        os.chdir(outp)
        if xtcfp.endswith("pkl3") or xtcfp.endswith("pkl") or xtcfp.endswith("pkl2"):
            meta = pd.read_pickle(xtcfp)
            xtcfs = meta["traj_fn"].to_numpy()
            top = meta["top_fn"].iloc[0]
        else:
            xtcfs = [os.path.join(xtcfp,_) for _ in os.listdir(xtcfp) if _.endswith(".xtc")]
        n_trajs = len(xtcfs)
        dtrajs = []
        if T >1:
            pool = Pool(T)
            args = [(xtcfs[i], cf, top, sel, dt) for i in range(n_trajs)]
            dtrajs = pool.map(parallel, args)
        else:
            for i in range(n_trajs):
                dtraj = assign_(xtcfs[i], cf, top, sel, dt=dt, outfname=None)
                dtrajs.append(dtraj)
        os.chdir(cwd0)
        pd.to_pickle(dtrajs, "dtrajs.pkl3")
if __name__ == '__main__':
    main()
