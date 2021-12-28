#!/usr/bin/env python
'''
@Author ymh
@Email  maohuay@hotmail.com
@Date   2019-10-21 21:02:27
@Web    https://github.com/Aunity
'''
import os
import sys
import argparse
import pandas as pd
import mdtraj as md

def parse_args():
    parser = argparse.ArgumentParser(description="Combin xtc files to dcd file.")
    parser.add_argument("-i", dest="inp", help="Trajs floder of meta pickle file.", required=True)
    parser.add_argument("-t", dest="top", help="Top file for trajs.", default=None)
    parser.add_argument("-dt", help="every dt step loaded. default: 1", type=int, default=1)
    parser.add_argument("-s", dest="sel", help="mdtraj protein selection.")
    parser.add_argument("-o", dest="outn", help="outfile name",required=True)

    args = parser.parse_args()
    return args.inp, args.top, args.sel, args.outn, args.dt

def main():
    inp, topf, sel, outn, dt = parse_args()
    if os.path.isdir(inp):
        xtcfs = sorted([os.path.join(inp,_) for _ in os.listdir(inp)])
        if topf is None:
            print("Eorror: top file is required!")
            exit(0)
    else:
        meta = pd.read_pickle(inp)
        topf = meta["top_fn"].values[0]
        xtcfs = meta["traj_fn"].values
    top = md.load_pdb(topf)
    ndx = top.top.select(sel)
    top = top.atom_slice(ndx)

    xtcall = None
    for xtcf in xtcfs:
        xtc = md.load_xtc(xtcf, top=topf, atom_indices=ndx, stride=dt)
        if xtcall is None:
            xtcall = xtc
        else:
            xtc.topology = xtcall.topology
            xtcall = xtcall.join(xtc)

    outtop = "%s.pdb"%outn
    outdcd = "%s.dcd"%outn
    top.save_pdb(outtop)
    xtcall.superpose(top)
    xtcall.save_dcd(outdcd)

if __name__ == '__main__':
    main()
