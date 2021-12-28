
import os
import argparse
import numpy as np
import mdtraj as md
import tqdm

'''
units: nm
'''

def cal_rg_PDB(pdbfile:str, selection='all'):
    """Calculate the Residue Gyration value for single PDB file.

    Args:
        pdbfile (str): Input PDB files.
        selection (str, optional): Selection group to calculate the Rg value. Defaults to 'all'.
    """
    PDB = md.load(pdbfile)
    AtomIndex = PDB.top.select(selection)
    selGroup = PDB.atom_slice(AtomIndex)
    rg = md.compute_rg(selGroup)
    print("%s: %.4f nm"%(pdbfile, rg))


def cal_rg_traj(topologyfile:str, trajfile:str, selection='all', outfile=None) -> np.array:
    """Calculate the Rg value for single trajecotry file.
    Args:
        topologyfile (str): topology file
        trajfile (str): trajectory file
        selection (str, optional): select action for atoms. Defaults to 'all'.
        outfile (str optional): outfile to save the Rg value.
    Returns:
        np.array: [description]
    """
    traj = md.load(trajfile, top=topologyfile)
    AtomIndex = traj.top.select(selection)
    rgs = md.compute_rg(traj.atom_slice(AtomIndex))
    if outfile:
        np.savetxt(outfile, rgs, fmt='%.4f')
    return rgs

def cal_rg_trajs(topologyfile:str, trajfiles:list, selection='all', outdir=None) ->list:
    """Calculate the Rg value for multiple trajectory files.
    Args:
        topologyfile (str): topology file
        trajfile (str): trajectory file
        selection (str, optional): select action for atoms. Defaults to 'all'.
        outdir (str, optional): 
    Returns:
        list: [description]
    """
    if outdir is None:
        outdir = 'Rg'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    RgDict = {}
    ProcessBar = tqdm.tqdm(trajfiles)
    Rgs = []
    for i, trajfile in enumerate(ProcessBar):
        ProcessBar.set_description("Process: %s "%trajfile)
        trajname = os.path.split(trajfile)[-1][:-4]
        outfile = os.path.join(outdir, trajname+'.rg')
        rg = cal_rg_traj(topologyfile, trajfile, selection, outfile=outfile)
        RgDict[i] = rg
        Rgs.append(rg)
    Rg = np.c_[Rgs]
    return RgDict

def hist_1d():
    """[summary]
    """
    #! TODO: plot 1D png
    pass

def main():
    parser = argparse.ArgumentParser(description='Calculate the RMSF value for ensemble.')
    parser.add_argument('-i', dest='INP', help='A simple pdbfile or a trajectory file or a floder contains trajectory files.', required=True)
    parser.add_argument('-t', dest='TOP', help='Topology file for the input trajectory file.')
    parser.add_argument('-s', dest='SEL', help='Selection atoms for calculate the RMSF.', default='mass >2')
    parser.add_argument('-o', dest='OUP', help='Output file or output directory.', default=None)
    
    trajFormat = '.xtc'
    args = parser.parse_args()
    #print(args)
    if args.TOP is None:
        cal_rg_PDB(args.INP, selection=args.SEL)
    elif os.path.isdir(args.INP):
        trajfiles = sorted([os.path.join(args.INP, trajfile) for trajfile in os.listdir(args.INP) if trajfile.endswith(trajFormat)])
        RgDict = cal_rg_trajs(args.TOP, trajfiles, selection=args.SEL, outdir=args.OUP)
        print(RgDict)
    else:
        Rgs = cal_rg_traj(args.TOP, args.INP, selection=args.SEL, outfile=args.OUP)
        print(Rgs)

if __name__ == "__main__":
    main()
    


