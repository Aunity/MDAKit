#!/usr/bin/env python
import os
import sys
import argparse

import numpy as np
import mdtraj as md
import pandas as pd

from cal_hb import cal_hb_xtc

def parse_args():
    parser = argparse.ArgumentParser(description="Cal hydrogen bond for msm ensemble.")
    parser.add_argument("-m", dest="msmf", help="The msm pickl file.", required=True)
    parser.add_argument("-s", dest="samf", help="The samples floder of each cluster.", required=True)
    parser.add_argument("-t", dest="topf", help="The top file for the trajectory file.", required=True)
    parser.add_argument("-o", dest="outf", help="The result floder.", default="hb-msm")

    args = parser.parse_args()

    return args.msmf, args.samf, args.topf, args.outf

def main():
    msmf, samf, topf, outf = parse_args()

    traj  = md.load(topf)
    M     = pd.read_pickle(msmf)
    if hasattr(M, "n_states_"):
        nstates = M.n_states_
        populations = M.populations_
    else:
        nstates = M.nstates
        populations = M.stationary_distribution

    if not os.path.exists(outf):
        os.mkdir(outf)
    angles_dict = {}
    distan_dict = {}
    hb_dict = {}
    hb_keys = []

    print("Cal hydrogen bond")
    n = 1

    for trajf in sorted(os.listdir(samf)):
        if not trajf.endswith(".xtc"):
            continue
        try:
            c_name = int(trajf[:-4])
        except:
            c_name = int(trajf[4:-4])
            #print("Error for the xtc file:%s"%trajf)
        abs_trajf = os.path.join(samf, trajf)
        hbdict, distances, angles = cal_hb_xtc(abs_trajf, topf)
        angles_dict[c_name] = angles
        distan_dict[c_name] = distances
        hb_dict[c_name] = hbdict
        hb_keys.extend(hbdict["key"].tolist())

        sys.stdout.write("Finished cal %s, %3d/%3d\r"%(trajf,n,nstates))
        sys.stdout.flush(); n += 1

    hb_keys = list(set(hb_keys))
    hb_prob = []
    hb_label= []

    label = lambda hbond : '%s--%s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))

    print("Merge hydrogen bond data")
    n = 0
    for k in hb_keys:
        k_prob = 0
        for i in range(nstates):
            hbdict = hb_dict[i]
            popu   = populations[i]

            tmp_hbcol = hbdict[hbdict['key']==k]
            if tmp_hbcol.empty:
                prob = 0
            else:
                prob = tmp_hbcol["probability"].values[0]
            k_prob += prob * popu
            #sys.stdout.write("del with key: %s,   %5d/%5d \r"%(k, n, len(hb_keys)))
            #sys.stdout.flush(); n += 1

        atom_i = list(map(int, k.split("-")))
        k_label= label(atom_i)

        hb_prob.append(k_prob)
        hb_label.append(k_label)

    hb_msm = pd.DataFrame({"label":hb_label, "frequence":hb_prob, "key":hb_keys})

    hb_msmf = os.path.join(outf, "hb-msm.txt")
    hb_anglesf = os.path.join(outf, "angles-msm.pickl")
    hb_distanf = os.path.join(outf, "distance-msm.pickl")

    hb_msm.to_csv(hb_msmf, sep=" ", index=False)
    pd.to_pickle(angles_dict, hb_anglesf)
    pd.to_pickle(distan_dict, hb_distanf)

if __name__ == "__main__":
    main()
