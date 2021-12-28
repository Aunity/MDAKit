import argparse
import os
import sys
from pandas import DataFrame
from mdtraj import load, compute_dssp
from MDAnalysis import Universe
from numpy import sum

# Full mode: 7
comment7 = [
    "#'H' : Alpha helix",
    "#'B' : Residue in isolated beta-bridge",
    "#'E' : Extended strand, participates in beta ladder",
    "#'G' : 3-helix (3/10 helix)",
    "#'I' : 5 helix (pi helix)",
    "#'T' : hydrogen bonded turn",
    "#'S' : bend",
    "#'L' : Loops and irregular elements"
    ]
commet3 = [
    "#‘H’ : Helix. Either of the ‘H’, ‘G’, or ‘I’ codes",
    "#‘E’ : Strand. Either of the ‘E’, or ‘B’ codes",
    "#‘C’ : Coil. Either of the ‘T’, ‘S’ or ‘ ‘ codes"
]

def cal_dssp_PDB(pdbfile:str, mode=3):
    if mode == 7:
        simplified = False
    elif mode == 3:
        simplified = True
    else:
        raise Exception('The mode parameter only accept 3 or 7. %s'%str(mode))
    PDB = load(pdbfile)[0]
    u = Universe(pdbfile)
    chainNames = u.segments.segids
    ssp = compute_dssp(PDB, simplified)[0]
    mask = ssp == ' '
    ssp[mask] = 'L'
    chainIDs = []
    resnames = []
    resids   = []
    for chain in PDB.top.chains:
        chainName = chainNames[chain.index]
        for residue in chain.residues:
            resnames.append(residue.name)
            resids.append(residue.resSeq)
            chainIDs.append(chainName)
    if len(chainIDs) != len(ssp):
        raise Exception('The chainID and ssp results must be same length. %d-%d'%(len(chainIDs), len(ssp)))
    ResultDF = DataFrame({
        'resid'  : resids,
        'resname': resnames,
        'chainID': chainIDs,
        'ssp'    : ssp
    })
    print(ResultDF)
    return ResultDF

def cal_dssp_traj(trajfile:str, topfile:str, mode=3):
    if mode == 7:
        simplified = False
        SStype = [ 'H', 'B', 'E', 'G', 'I', 'T', 'S', 'L']
        comment = comment7
    elif mode == 3:
        simplified = True
        SStype = ['H', 'E', 'C']
        comment = commet3
    else:
        raise Exception('The mode parameter only accept 3 or 7. %s'%str(mode))
    top = load(topfile)[0]
    u = Universe(topfile)
    chainNames = u.segments.segids
    chainIDs = []
    resnames = []
    resids   = []
    for chain in top.top.chains:
        chainName = chainNames[chain.index]
        for residue in chain.residues:
            resnames.append(residue.name)
            resids.append(residue.resSeq)
            chainIDs.append(chainName)
    traj = load(trajfile, top=topfile)
    ssp = compute_dssp(traj, simplified=simplified)
    mask = ssp == ' '
    ssp[mask] == 'L'
    if len(chainIDs) != len(ssp[0]):
        raise Exception('The chainID and ssp results must be same length. %d-%d'%(len(chainIDs), len(ssp[0])))
    ResultDict = {
        'resid'  : resids,
        'resname': resnames,
        'chainID': chainIDs,
    }
    residueSize = ssp.size
    for ss in SStype:
        ResultDict[ss] = sum(ssp==ss, axis=0)
        comment.append("# %s: %.4f"%(ss, sum(ssp==ss)/residueSize))
    ResultDF = DataFrame(ResultDict)
    print("\n".join(comment))
    print(ResultDF)
    return ResultDict, ssp
    
def main():
    parser = argparse.ArgumentParser(description='Calculate the secondary information for a PDB file.')
    parser.add_argument('-i', dest='PDB', help='The input PDB file.', required=True)
    parser.add_argument('-t', dest='TRAJ', help='Trajectory file to calculate the ssp.', default=None)
    parser.add_argument('-m', dest='MODE', help='Simple mode or complex mode.', default=3, choices=[3,7], type=int)
    
    args = parser.parse_args()
    if args.TRAJ:
        cal_dssp_traj(args.TRAJ, args.PDB, args.MODE)
    else:
        cal_dssp_PDB(args.PDB, args.MODE)

if __name__ == "__main__":
    main()
