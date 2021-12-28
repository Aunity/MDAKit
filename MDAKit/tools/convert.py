#!/usr/bin/env python
import os
import argparse
import mdtraj as md
import MDAnalysis as mda
from tqdm import tqdm


def convert_traj_to_pdbs(trajfile:str, topologyfile:str, selection='all', outdir=None) -> None:
    """Convert trajectory file to multiple PDB files

    Args:
        trajfile (str): Trajectory file.
        topologyfile (str): topology file for trajectoy file.
        selection (str, optional): Selection to output. Defaults to 'all'.
        outdir ([type], optional): Output directory. Defaults to None.
    """
    u = mda.Universe(topologyfile, trajfile)
    frameN = u.trajectory.n_frames
    format = 10**(len(str(frameN))+1)
    if outdir is None:
        outdir = os.path.split(trajfile)[-1][:-4]
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    pbar = tqdm(u.trajectory)
    for i, ts in enumerate(pbar):
        outfile = os.path.join(outdir, "%d.pdb"%(format+i))
        pbar.set_description('Processing: %s'%outfile)
        SelGroup = u.select_atoms(selection)
        SelGroup.write(outfile)
    print('Output directory: %s'%os.path.abspath(outdir))

def convert_pdbs_to_traj(pdbdir:list, selection='all', outfile=None) -> None:
    """Convert multiple pdb files to a single trajectory file.

    Args:
        pdbdir (list): A floder contains pdb files.
        selection (str, optional): [description]. Defaults to 'all'.
        outfile ([type], optional): [description]. Defaults to None.
    """
    fileformat = '.pdb'
    pdbfiles = sorted([os.path.join(pdbdir, pdbfile) for pdbfile in os.listdir(pdbdir) if pdbfile.endswith(fileformat)])
    u = mda.Universe(pdbfiles[0], pdbfiles)
    SelGroup = u.select_atoms(selection)
    pbar = tqdm(pdbfiles)
    pdbNum = len(pdbfiles)
    if outfile is None:
        pdbname = pdbdir.rstrip(os.path.sep).split(os.path.sep)[-1]
        print(pdbname)
        outfile = '%s_%d.xtc'%(pdbname, pdbNum)
    outTopfile = outfile[:-4] + '.pdb'
    SelGroup.write(outTopfile)
    with mda.Writer(outfile, SelGroup.n_atoms) as W:
        for i, pdbfile in enumerate(pbar):
            pbar.set_description('Processing: %s'%pdbfile)
            W.write(SelGroup)
    print("Output trajectory file: %s"%os.path.abspath(outfile))
    print("Output topology file: %s"%os.path.abspath(outTopfile))

def main():
    parser = argparse.ArgumentParser(description="Convert single trajectory file to multiple PDB file or Convert multiple PDB files to single trajectory file.")
    parser.add_argument('-i', dest='INP', help='Input trajectory file or PDB floder.', required=True)
    parser.add_argument('-t', dest='TOP', help='Topology file for the trajectory file.')
    parser.add_argument('-s', dest='SEL', help='selection str for group select.', default='all')
    parser.add_argument('-o', dest='OUP', help='Output file or output floder.', default=None)

    args = parser.parse_args()
    if os.path.isdir(args.INP):
        convert_pdbs_to_traj(args.INP, args.SEL, args.OUP)
    else:
        convert_traj_to_pdbs(args.INP, args.TOP, args.SEL, args.OUP)

if __name__ == "__main__":
    main()

