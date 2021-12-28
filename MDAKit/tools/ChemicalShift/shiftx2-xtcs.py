#!/software/anaconda2/envs/msmb/bin/python

import os
import sys
import argparse

import mdtraj as md
import pandas as pd
import numpy as np

from multiprocessing import Pool


def job(cmd):
    os.system(cmd)

def shift_xtc(xtc_p, top,ph,shiftx1, T=1):
    xtcs = [os.path.join(xtc_p,xtc) for xtc in os.listdir(xtc_p) if xtc.endswith('.xtc')]
    cmds = []
    if shiftx1:
        shiftx1 = "-shiftx1"
    else:
        shiftx1 = ""
    for xtc in xtcs:
        cmd = "shiftx2-xtc.py -x %s -t %s -p %f %s"%(xtc,top,ph,shiftx1)
        cmds.append(cmd)
    p = Pool(T)
    p.map(job,cmds)

def parse_arg():
    parser = argparse.ArgumentParser(description='Calculate chemical shift for xtc')
    parser.add_argument('-x', dest='xtc', help="xtcis file for cal chemical shift", required=True)
    parser.add_argument('-t', dest='top', help="top file for the xtc file, default=/home/ymh/GGBP/2fvy.pdb",default='/home/ymh/GGBP/2fvy.pdb')
    parser.add_argument('-T', dest='tlp', help="number of threads", default=1, type=int)
    parser.add_argument('-p', dest="ph", help="PH value, default PH=7.0", default=7.0, type=float)
    parser.add_argument('-shiftx1', help="Version of shiftx1.", action="store_true")
    args = parser.parse_args()
    return args.xtc, args.tlp, args.top, args.ph,args.shiftx1

def main():
    xtc_p,T,top,ph,shiftx1 = parse_arg()
    shift_xtc(xtc_p,top,ph,shiftx1,T)
if __name__ == "__main__":
    main()
