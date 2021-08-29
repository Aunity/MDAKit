#!/usr/bin/env python
'''
@Author ymh
@Email  maohuay@hotmail.com
@Date   2020-10-20 15:33:29
@Web    https://github.com/Aunity
'''
import os
import sys
import argparse

import numpy as np
import mdtraj as md
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Cal hydrogen bond for xtc or pdb file.")
    parser.add_argument("-i", dest="inp", help="Trajectory file for aligen to reference structure. xtc, dcd etc.", required=True)
    parser.add_argument("-t", dest="topf", help="Top file for the trajectory file.",required=True)
    #parser.add_argument("-r", dest="reff", help="Reference for the trajectories align to.",required=True)
    parser.add_argument("-s", dest="sel", help="Selection for topology", default='all')
    parser.add_argument("-o", dest="oup", help="The result floder. default:.", default="aligned")

    args = parser.parse_args()

    return args.inp, args.topf, args.sel, args.oup

def superpose(trajf, topf, outf, sel):
    xtc = md.load(trajf, topf)
    ndx = xtc.top.select(sel)
    xtc = xtc.atom_slice(ndx)
    ref = md.load_pdb(topf)
    ref = ref.atom_slice(ndx)
    xtc = xtc.superpose(ref)
    xtc.save_xtc(outf)

def main():
    inp, topf, sel, outf = parse_args()
    if os.path.isdir(inp):
        if not os.path.eixsts(oup):
            os.mkdir(oup)
        fnames = [_ for _ in os.listdir(inp)]
        for fname in fnames:
            trajf = os.path.join(inp, fname)
            outf = os.path.join(oup, fname)
            superpose(trajf, topf, outf, sel)
    else:
        superpose(inp, topf, oup, sel)

if __name__ == "__main__":
    main()
