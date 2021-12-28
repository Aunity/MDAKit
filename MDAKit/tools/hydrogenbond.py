#!/usr/bin/env python
import os
import argparse
import logging
import numpy as np
import mdtraj as md
import pandas as pd
from tqdm import tqdm


def cal_hb_PDB(pdbfile:str, outdir=None):
    hbdict, distances, angles = cal_hb_traj(pdbfile, pdbfile, outdir=outdir)

def cal_hb_traj(trajfile:str, topfile:str, outdir=None, exclude_water=True, periodic=False, sidechain_only=False, save=True):
    """[summary]

    Args:
        trajfile (str): [description]
        topfile (str): [description]
        outdir ([type], optional): [description]. Defaults to None.
        exclude_water (bool, optional): [description]. Defaults to True.
        periodic (bool, optional): [description]. Defaults to False.
        sidechain_only (bool, optional): [description]. Defaults to False.

    Returns:
        [type]: [description]
    """
    distance_cutoff = 0.25            # nanometers
    angle_cutoff = 2.0 * np.pi / 3.0  # radians
    distance_indices = [1, 2]
    angle_indices    = [0, 1, 2]

    if trajfile == topfile:
        traj = md.load(topfile)
    else:
        traj = md.load(trajfile, top=topfile)
    triplets = md.geometry.hbond._get_bond_triplets(traj.topology,
                                                         exclude_water=exclude_water, sidechain_only=sidechain_only)
    if len(triplets)==0:
        logging.error(' Make sure your topology file contains Hydrogen atoms.')
        exit(1)
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
    hbdict = pd.DataFrame({"label":atoml, "atom-1":triplets[:,0], "atom-2":triplets[:,1], "atom-3":triplets[:,2], "frequence":count, "probability":prob, "key":dict_key})

    hbdict.n_frame = traj.n_frames

    if save:
        trajname = os.path.split(trajfile)[-1][:-4]
        if outdir is None:
            outdir = trajname
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        hbdf   = os.path.join(outdir, "hbdict.txt")
        anglef = os.path.join(outdir, "angle.pickl")
        distf  = os.path.join(outdir, "distance.pickl")

        hbdict.to_csv(hbdf, sep=" ", index=False)
        pd.to_pickle(angles, anglef)
        pd.to_pickle(distances, distf)

    return hbdict, distances, angles

def cal_hb_trajs(trajdir:str, topfile:str, outdir=None):
    fileFormat = '.xtc'
    traj  = md.load_pdb(topfile)
    trajfiles = [ os.path.join(trajdir, trajfile) for trajfile in sorted(os.listdir(trajdir)) if trajfile.endswith(fileFormat)]
    sep = os.path.sep
    if outdir is None:
        outdir = trajdir.rstrip(sep).split(sep)[-1] + '_hb'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    angles_dict = {}
    distan_dict = {}
    hb_dict = {}
    hb_keys = []

    print("Cal hydrogen bond:")
    pbar = tqdm(trajfiles)
    for i,trajfile in enumerate(pbar):
        pbar.set_description('Processing: %s'%trajfile)
        trajname = os.path.split(trajfile)[-1][:-4]
        hbdict, distances, angles = cal_hb_traj(trajfile, topfile, save=False)
        angles_dict[trajname] = angles
        distan_dict[trajname] = distances
        hb_dict[trajname] = hbdict
        hb_keys.extend(hbdict["key"].tolist())

    hb_keys = list(set(hb_keys))
    hb_prob = []
    hb_label= []

    label = lambda hbond : '%s--%s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))

    print("Merge hydrogen bond data:")
    pbar = tqdm(hb_keys)
    for k1 in pbar:
        k_prob = 0
        pbar.set_description('Processing: %s'%k1)
        for k2 in hb_dict.keys():
            hbdict = hb_dict[k2]
            popu   = 1
            tmp_hbcol = hbdict[hbdict['key']==k1]
            if tmp_hbcol.empty:
                prob = 0
            else:
                prob = tmp_hbcol["probability"].values[0]
            k_prob += prob * popu
            #sys.stdout.write("del with key: %s,   %5d/%5d \r"%(k, n, len(hb_keys)))
            #sys.stdout.flush(); n += 1

        atom_i = list(map(int, k1.split("-")))
        k_label= label(atom_i)

        hb_prob.append(k_prob/(len(trajfiles)*1.0))
        hb_label.append(k_label)

    hb_msm = pd.DataFrame({"label":hb_label, "probability":hb_prob, "key":hb_keys})

    hb_msmf = os.path.join(outdir, "hb-trajs.txt")
    hb_anglesf = os.path.join(outdir, "angles-trajs.pickl")
    hb_distanf = os.path.join(outdir, "distance-trajs.pickl")

    hb_msm.to_csv(hb_msmf, sep=" ", index=False)
    pd.to_pickle(angles_dict, hb_anglesf)
    pd.to_pickle(distan_dict, hb_distanf)



def main():
    parser = argparse.ArgumentParser(description="Calculate the hydrogen bond information.")
    parser.add_argument('-i', dest='INP', help='Input PDB file or Input trajectory file or Input floder contains trajectories.', required=True)
    parser.add_argument('-t', dest='TOP', help='Topology file for trajectory file, need for INP is trajectory.', default=None)
    parser.add_argument('-o', dest='OUP', help='Output directory.', default=None)
    
    args = parser.parse_args()
    if os.path.isdir(args.INP):
        cal_hb_trajs(args.INP, args.TOP, args.OUP)
    elif os.path.isfile(args.INP) and args.TOP:
        cal_hb_traj(args.INP, args.TOP, args.OUP)
    else:
        cal_hb_PDB(args.INP, args.OUP)

if __name__ == "__main__":
    main()