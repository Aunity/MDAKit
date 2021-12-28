#!/usr/bin/env python
import os
import sys
import argparse

import numpy as np
import mdtraj as md
import pandas as pd

from multiprocessing import Pool

def func(arg):
    xtcf, topf, reffs, p, sel = arg
    toppdb = md.load_pdb(topf)
    atm0 = toppdb.top.select(sel)
    xtc = md.load_xtc(xtcf, topf)
    xtc = xtc.atom_slice(atm0)

    rmsd = []
    for reff in reffs:
        refpdb = md.load_pdb(reff)
        atm1 = refpdb.top.select(sel)
        ref = refpdb.atom_slice(atm1)
        rmsd.append(md.rmsd(xtc, ref))
    rmsd.append([p]*len(xtc))
    return np.array(rmsd).T

def parse_args():
    parser = argparse.ArgumentParser(description="Cal a MSM ensemble's RMSD distribution.")
    parser.add_argument('-m', dest='msmf', help="msmbuilder 3*'s MSM pickl model. type:pickl", required=True)
    parser.add_argument('-s', dest='samplesp', help="MSM samples xtc files' floder.")
    parser.add_argument('-t', dest='topf', help="The top file for the xtc file.")
    parser.add_argument('-r', dest='reff', help="The refrence to cal rmsd. default is the topf.", nargs="+", default=None)
    parser.add_argument('-sel', dest='prosel', help="The selection to cal RMSD. default=name CA and protein", default="protein and name CA")
    parser.add_argument('-o', dest='outf', help="Output file name")
    parser.add_argument('-T', dest='tlp', help="The cpu core number default:1", default=1, type=int)

    args = parser.parse_args()

    return args.msmf, args.samplesp, args.topf, args.prosel, args.tlp, args.outf, args.reff

def main():
    msmf, samplesp, topf, sel, T, outf, reff = parse_args()
    if reff is None:
        reff = [topf]
    M = pd.read_pickle(msmf)
    if hasattr(M, "populations_"):
        populations = M.populations_/M.populations_.sum()
        fmt = 1000
        args = [ (os.path.join(samplesp,"%d.xtc"%(i+fmt)), topf, reff, populations[i],sel) for i in range(M.n_states_) ]
    else:
        populations = M.stationary_distribution/M.stationary_distribution.sum()
        trajfs = sorted([ os.path.join(samplesp, _) for _ in os.listdir(samplesp) if _.endswith(".xtc") ])
        args = [ (trajfs[i], topf, reff, populations[i],sel) for i in range(M.nstates) ]
    #top = md.load_pdb(reff)
    #atom_index = top.top.select(sel)
    pool = Pool(T)
    result = pool.map(func, args)
    names = [ os.path.split(ref)[1].replace(".pdb","") for ref in reff ]
    header = " ".join(names) + " weight"
    np.savetxt(outf, np.concatenate(result), header=header, fmt="%.5e")
    #pd.to_pickle(result, "result.pkl3")
    #txx = np.concatenate(np.array(result), axis=1).T
    #pd.to_pickle(txx, outf)

if __name__ == "__main__":
    main()
