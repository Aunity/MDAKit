#!/usr/bin/env python
import os
import sys

import argparse
import numpy as np
import pandas as pd
import mdtraj as md

from multiprocessing import Pool

def func(args):
    traj, top, reffs, sel = args
    if not isinstance(traj, md.core.trajectory.Trajectory):
        traj = md.load(traj, top=top)
    atm_ndx = traj.top.select(sel)
    traj = traj.atom_slice(atm_ndx)
    rmsds = []
    for reff in reffs:
        ref = md.load_pdb(reff)
        atm_ndx = ref.top.select(sel)
        ref = ref.atom_slice(atm_ndx)
        traj = traj.superpose(ref)
        rmsd = md.rmsd(traj, ref)
        rmsds.append(rmsd)

    return np.array(rmsds)

def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", dest="inp", help="Input file. single pdb file, trajectory file or meta trajectories", required=True)
    parser.add_argument("-t", dest="top", help="Topology file for the trajectory file, needed if the inp is single trajectory")
    parser.add_argument("-r", dest="ref", help="Reference structure to cal RMSD.", nargs="+")
    parser.add_argument("-sel", help="Atom select to align and cal RMSD. default: protein and name CA", default="protein and name CA")
    parser.add_argument("-o", dest="oup", help="Output file or floder to save result")

    args = parser.parse_args()

    return args.inp, args.top, args.ref, args.sel, args.oup

def main():
    print("start cal.")
    inp, topf, reff, sel, oup = parse_args()
    T = 4
    if inp.endswith(".pdb"):
        pdb = md.load_pdb(inp)
        dCNH = func([pdb,inp,reff,sel])
        print(dCNH)
    elif inp.endswith(".xtc"):
        xtc = md.load_xtc(inp, top=topf)
        dCNH = func([xtc, topf,reff, sel])
        print(dCNH)
        np.savetxt(oup, dCNH.T, fmt="%.4f",delimiter=',')
    elif inp.startswith("meta"):
        meta = pd.read_pickle(inp)
        topfn = meta["top_fn"].iloc[0]
        files = list(meta['traj_fn'])
        args = [(fi, topfn, reff, sel) for fi in files]
        pool = Pool(T)
        dCNHs = pool.map(func, args)
        if not os.path.exists(oup):
            os.mkdir(oup)
        names = []
        for ref in reff:
            name = os.path.split(ref)[-1]
            names.append(name)
        title = ",".join(names)
        for dCNH, arg in zip(dCNHs, args):
            fi = arg[0]
            fname = os.path.split(fi)[-1] + ".txt"
            fname = os.path.join(oup, fname)
            np.savetxt(fname, dCNH.T, fmt="%.4f",delimiter=',',header=title)
    elif os.path.isdir(inp):
        files = [ os.path.join(inp, _) for _ in os.listdir(inp) if _.endswith(".xtc") ]
        args = [(fi, topf,reff,sel) for fi in files]
        pool = Pool(T)
        #func(args[0])
        #return
        dCNHs = pool.map(func, args)
        if not os.path.exists(oup):
            os.mkdir(oup)
        names = []
        for ref in reff:
            name = os.path.split(ref)[-1]
            names.append(name)
        title = ",".join(names)
        for dCNH, arg in zip(dCNHs, args):
            fi = arg[0]
            fname = os.path.split(fi)[-1] + ".txt"
            fname = os.path.join(oup, fname)
            np.savetxt(fname, dCNH.T, fmt="%.4f",delimiter=',',header=title)
    else:
        print("Error: check the input file.")
    #meta = pd.read_pickle(metaf)
    #top =

if __name__ == "__main__":
    main()
