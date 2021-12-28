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

def maxmatrix(mat):
    ndx0 = mat.argmax(axis=1)
    maxv = mat.max(axis=1)
    y = maxv.argmax()
    v = maxv.max()
    x = ndx0[y]
    return x,y,v

def parallel_job(args):
    xtcf0, xtcf1, top, sel = args
    rmsd,v,x,y = cal_rmsdmatrix(xtcf0, xtcf1, top, sel)
    return rmsd,v,x,y

def cal_rmsdmatrix(xtcf0, xtcf1, top, sel, dt=1, outfname=None):
    if not outfname:
        fname0,fname1 = os.path.split(xtcf0)[1].split(".")[0],os.path.split(xtcf1)[1].split(".")[0]
        outfname = "%s_%s"%(fname0,fname1)
    top = md.load_pdb(top)
    ndx = top.top.select(sel)
    xtc0 = md.load_xtc(xtcf0, top=top, atom_indices=ndx, stride=dt)
    if xtcf0==xtcf1:
        xtc1 = xtc0
    else:
        xtc1 = md.load_xtc(xtcf1, top=top, atom_indices=ndx, stride=dt)
    rmsd = libdistance.cdist(xtc0, xtc1, metric="rmsd")
    rmsd = np.triu(rmsd)
    #np.save(outfname+".npy", rmsd)
    #from scipy import sparse
    #b = sparse.csr_matrix(rmsd)
    #sparse.save_npz('b_compressed.npz', b, True)
    x,y,v = maxmatrix(rmsd)
    np.savez_compressed(outfname+'.npz', a=rmsd, b=np.array([x,y,v]))
    return rmsd,v,x,y

def parse_args():
    parser = argparse.ArgumentParser(description="Cal RMSD matrix for trajectory(s).")
    parser.add_argument("-i", dest="inp", help="Trajectory file or floder contains trajectories, or meta pkl3 file contains infor.",required=True)
    parser.add_argument("-t", dest="top", help="Topology file for trajectory.", default=".")
    parser.add_argument("-s", dest="dt", help="every dt time",  type=int, default=1)
    parser.add_argument("-sel", help="mdtraj select str to extract atoms for cal RMSD. defatult: name CA", default="name CA")
    parser.add_argument("-o", dest="oup", help="A floder or a name to save the RMSD matrix npy file.", default="rmsd_mat")
    parser.add_argument("-T", help="Number of thread to cal RMSD. default=1", type=int, default=1)

    args = parser.parse_args()
    return os.path.abspath(args.inp), os.path.abspath(args.top), args.sel, args.oup, args.T, args.dt

def main():
    xtcfp,top,sel,outp,T,dt = parse_args()
    vmax = -999
    if os.path.isfile(xtcfp) and xtcfp.endswith(".xtc"):
        rmsd,v,fi,fj = cal_rmsdmatrix(xtcfp, xtcfp, top, sel, dt)
        if v>=vmax:
            trajmax = [xtcfp,fi,fj,v]
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
        for i in range(n_trajs):
            for j in range(i, n_trajs):
                rmsd,v,x,y = cal_rmsdmatrix(xtcfs[i], xtcfs[j], top, sel, dt)
                if v>=vmax:
                    trajmax = [xtcfs[i], xtcfs[j], v]
                    vmax = v
        os.chdir(cwd0)
    print(trajmax)
if __name__ == '__main__':
    main()
