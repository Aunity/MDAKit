#!/usr/bin/env python
import os
import sys
import argparse

import numpy as np
import mdtraj as md
import pandas as pd

def cal_hb_xtc(trajf, topf, exclude_water=True, periodic=False, sidechain_only=False):
    distance_cutoff = 0.25            # nanometers
    angle_cutoff = 2.0 * np.pi / 3.0  # radians
    distance_indices = [1, 2]
    angle_indices    = [0, 1, 2]

    traj = md.load(trajf, top=topf)
    triplets = md.geometry.hbond._get_bond_triplets(traj.topology,
                                                         exclude_water=exclude_water, sidechain_only=sidechain_only)
    # First we calculate the requested distances
    distances = md.compute_distances(traj, triplets[:, distance_indices], periodic=periodic)

    # Now we discover which triplets meet the distance cutoff often enough
    prevalence = np.mean(distances < distance_cutoff, axis=0)
    mask = prevalence > 0

    # Update data structures to ignore anything that isn't possible anymore
    triplets = triplets.compress(mask, axis=0)
    distances = distances.compress(mask, axis=1)

    # Calculate angles using the law of cosines
    angle_indices = [0,1,2]
    abc_pairs = zip(angle_indices, angle_indices[1:] + angle_indices[:1])
    abc_distances = []

    # Calculate distances (if necessary)
    for abc_pair in abc_pairs:
        if set(abc_pair) == set(distance_indices):
             abc_distances.append(distances)
        else:
            abc_distances.append(md.compute_distances(traj, triplets[:, abc_pair],
                                periodic=periodic))

    # Law of cosines calculation
    a, b, c = abc_distances
    cosines = (a ** 2 + b ** 2 - c ** 2) / (2 * a * b)
    np.clip(cosines, -1, 1, out=cosines) # avoid NaN error
    angles = np.arccos(cosines)

    # Find triplets that meet the criteria
    presence = np.logical_and(distances < distance_cutoff, angles > angle_cutoff)
    mask     = np.mean(presence, axis=0) > 0
    count    = np.sum(presence, axis=0)

    triplets = triplets.compress(mask, axis=0)
    distances= distances.compress(mask, axis=1)
    angles   = angles.compress(mask, axis=1)
    count    = count.compress(mask, axis=0)
    prob     = count/float(traj.n_frames)


    label = lambda hbond : '%s--%s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
    atoml = [ label(hbond) for hbond in triplets ]
    dict_key = [ "%.0f-%.0f-%.0f"%tuple(hbond) for hbond in triplets ]
    hbdict = pd.DataFrame({"labels":atoml, "atom-1":triplets[:,0], "atom-2":triplets[:,1], "atom-3":triplets[:,2], "frequence":count, "probability":prob, "key":dict_key})

    hbdict.n_frame = traj.n_frames

    return hbdict, distances, angles

def parse_args():
    parser = argparse.ArgumentParser(description="Cal hydrogen bond for xtc or pdb file.")
    parser.add_argument("-x", dest="trajf", help="Trajectory file for cal hydrogen bond. xtc, dcd etc.", required=True)
    parser.add_argument("-t", dest="topf", help="Top file for the trajectory file.",required=True)
    parser.add_argument("-o", dest="outf", help="The result floder. default:.", default=".")

    args = parser.parse_args()

    return args.trajf, args.topf, args.outf

def main():
    trajf, topf, outf = parse_args()
    hbdict, distances, angles =  cal_hb_xtc(trajf, topf)

    if outf and not os.path.exists(outf):
        os.mkdir(outf)

    hbdf   = os.path.join(outf, "hbdict.txt")
    anglef = os.path.join(outf, "angle.pickl")
    distf  = os.path.join(outf, "distance.pickl")

    hbdict.to_csv(hbdf, sep=" ", index=False)
    pd.to_pickle(angles, anglef)
    pd.to_pickle(distances, distf)

if __name__ == "__main__":
    main()
