#!/usr/bin/env python

import os
import sys
import prody
import argparse

import numpy as np
import pandas as pd

from multiprocessing import Pool

def assign_pcs(args):
    fn, topf, eda, pcs, sel, outf = args

    if fn.endswith("pdb"):
        pdb = prody.parsePDB(fn)
        pdb = pdb.select(sel).copy()

        ensemble = prody.Ensemble('A single pdb file ensemble')
        ensemble.setCoords( pdb.getCoords() )
        ensemble.addCoordset( pdb.getCoordsets() )
        ensemble.iterpose()

        PCs = prody.calcProjection(ensemble, eda[pcs])
        print(PCs)
        return
    elif fn.endswith(".dcd"):

        structure = prody.parsePDB(topf)
        str_sel   = structure.select(sel)

        #dcd = prody.DCDFile(fn)
        dcd = prody.Trajectory(fn)
        dcd.link(structure)
        dcd.setCoords(structure)
        dcd.setAtoms(str_sel)

        PCs = prody.calcProjection(dcd, eda[pcs])
        if outf is not None:
            header = " ".join([ "PC%d"%(i+1) for i in pcs ])
            np.savetxt(outf, PCs, fmt="%.4f", header=header, comments="")
    else:
        print("Unsupport file type: %s"%fn)
        return None
    return PCs

def parse_args():
    parser = argparse.ArgumentParser(description="Assign the MD trajectory or pdb file to PCs.")

    parser.add_argument('-i', dest="inp", help="A single dcd or pdb file to assign or a floder contains dcd files.", required=True)
    parser.add_argument('-t', dest="top", help="Top file for the dcd trajectory file.")
    parser.add_argument('-e', dest="eda", help="The eda pickl file.", required=True)
    parser.add_argument('-s', dest="sel", help="The select string for protein.", required=True)
    parser.add_argument('-p', dest="pcs", type=int, nargs='+', help="The pcs to assign start with 0, default: all", default=None)
    parser.add_argument('-o', dest="out", help="The output file path. default=./temp", default="temp")
    parser.add_argument("-T", help="Number of thread to fun this job. default:1", type=int, default=1)
    args = parser.parse_args()

    return args.inp, args.top, args.eda, args.sel, args.pcs, args.out, args.T

if __name__ == "__main__":
    inpf, topf, edaf, sel, pcs, outf, T = parse_args()

    eda = pd.read_pickle(edaf)
    if pcs is None:
        pcs = np.arange(eda.numModes())

    if os.path.isfile(inpf):
        assign_pcs((inpf, topf, eda, pcs, sel, os.path.split(inpf)[-1]+".PC"))
    else:
        if not os.path.exists(outf):
            os.mkdir(outf)
        args = [ ( os.path.join(inpf, dcdf), topf, eda, pcs, sel, os.path.join(outf, dcdf+".PC")) for dcdf in os.listdir(inpf) ]
        pool = Pool(T)
        PCs  = pool.map(assign_pcs, args)
        pd.to_pickle(PCs, os.path.join(outf, "PCs-all.pickl"))
