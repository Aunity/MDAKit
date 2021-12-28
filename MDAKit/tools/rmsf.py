
import os
import argparse

import tqdm
import pandas as pd
import MDAnalysis as mda

'''
Calculat RMSF value for xtc
'''

import mdtraj as md

def cal_rmsf_traj(topologyfile:str, trajfile:str, selection='mass >= 2', mode='residue', outfile=None) -> pd.DataFrame:
    """Calculate the root mean square fluctuation for the trajfile.

    Args:
        topologyfile (str): topology file for the trajectory
        trajfile (str): molecular dynamic simulation trajectory file.
        selection (str, optional): atom index select to calculate the RMSF. Defaults to 'mass >= 2'.
        mode (str, optional): calculation mode, atom or residue. Defaults to 'residue'.
        outfile (str, optional): default is None, outfile to save result.
    Raises:
        Exception: mode not reisude or atom
    Returns:
        pd.DataFrame: RMSF value with data frame
    """
    univer = mda.Universe(topologyfile)
    ChainNames = univer.segments.segids
    traj = md.load(trajfile, top=topologyfile)
    AtomIndex = traj.topology.select(selection)
    TrajSelect = traj.atom_slice(AtomIndex)
    # target, reference, frame of reference
    rmsf = md.rmsf(TrajSelect, TrajSelect, 0)
    topology = TrajSelect.topology
    ColumnNames = ['chain', 'resid', 'resname', 'atomid', 'AtomName', 'RMSF']
    modes = ['residue', 'atom']
    if mode not in modes:
        raise Exception('Unsupport mode: %s, only accept residue or atom.'%mode)
    
    index = 0
    records = []
    for chain in topology.chains:
        for atom in chain.atoms:
            record = (ChainNames[chain.index], atom.residue.resSeq, atom.residue.name, atom.serial, atom.name, rmsf[index])
            records.append(record)
            index += 1
            #RMSFdict['resid'].append(atom.residue.resSeq)
            #RMSFdict['atom']
    RMSFdf = pd.DataFrame(records, columns=ColumnNames)
    if mode == 'residue':
        RMSFdf = RMSFdf.groupby(by=['chain','resid', 'resname'], as_index=False).agg({'RMSF':'mean'})
    if outfile:
        RMSFdf.to_csv(outfile, index=False)
    return RMSFdf

def cal_rmsf_trajs(topologyfile:str, trajfiles:list, selection='mass>=2', mode='residue', outdir=None) -> pd.DataFrame:
    """Calculate the root mean square fluctuation for multiple trajfile.
    Args:
        topologyfile (str): topology file for the trajectory
        trajfiles (list): molecular dynamic simulation trajectory files.
        selection (str, optional): atom index select to calculate the RMSF.. Defaults to 'mass>=2'.
        mode (str, optional): calculation mode, atom or residue. Defaults to 'residue'.
        outdir (str, optional): default is None, output result directory.
    Returns:
        pd.DataFrame: [description]
    """
    trajRMSF = {}
    RMSFdict = None
    if outdir is not None:
        os.mkdir(outdir)
    ProcessBar = tqdm.tqdm(trajfiles)
    for i, trajfile in enumerate(ProcessBar):
        ProcessBar.set_description("Process: %s "%trajfile)
        trajname = os.path.split(trajfile)[-1][:-4]
        outfile = os.path.join(outdir, trajname+'.rmsf')
        trajRMSF[i] = cal_rmsf_traj(topologyfile, trajfile, selection, mode, outfile)
        if RMSFdict is not None:
            RMSFdict['RMSF'] = RMSFdict['RMSF'].values + trajRMSF[i]['RMSF'].values
        else:
            RMSFdict = trajRMSF[i]
    if RMSFdict is not None:
        RMSFdict['RMSF'] = RMSFdict['RMSF'].values/len(trajfiles)
    return RMSFdict, trajRMSF

def cal_rmsf_MSM(topologyfile:str, stateSamples:str, populations:list, selection='mass >2', mode='residue') -> pd.DataFrame:
    """Calculate the root mean square fluctuation for multiple trajfile.
    Args:
        topologyfile (str): topology file for the trajectory
        stateSamples (str): state sample trajectory files generated from makov state model.
        populations (list): markov state model population for each state.
        selection (str, optional): atom index select to calculate the RMSF.. Defaults to 'mass>=2'.
        mode (str, optional): calculation mode, atom or residue. Defaults to 'residue'.
    Returns:
        pd.DataFrame: [description]
    """
    stateRMSF = {}
    RMSFdict = None
    ProccessBar = tqdm.tqdm(stateSamples)
    for i, stateSample in enumerate(ProccessBar):
        ProccessBar.set_description("Process: %s " % stateSample)
        stateRMSF[i] = cal_rmsf_traj(topologyfile, stateSample, selection, mode)
        if RMSFdict is not None:
            pi = populations[i]
            RMSFdict['RMSF'] = RMSFdict['RMSF'].values + stateRMSF[i]['RMSF'].values * pi
        else:
            RMSFdict['RMSF'] = stateRMSF[i] * pi
    return RMSFdict, stateRMSF

def main():
    parser = argparse.ArgumentParser(description='Calculate the RMSF value for ensemble.')
    parser.add_argument('-i', dest='INP', help='Input trajectory file or a floder contains trajectory files.', required=True)
    parser.add_argument('-t', dest='TOP', help='Topology file for the input trajectory file.', required=True)
    parser.add_argument('-s', dest='SEL', help='Selection atoms for calculate the RMSF.', default='mass >2')
    parser.add_argument('-m', dest='MODE', help='Calculation mode: residue or atom, default is residue.', choices=['residue','atom'], default='residue')
    parser.add_argument('-o', dest='OUP', help='Output file or output directory.', default=None)
    
    trajFormat = '.xtc'
    args = parser.parse_args()
    if os.path.isdir(args.INP):
        trajfiles = sorted([os.path.join(args.INP, trajfile) for trajfile in os.listdir(args.INP) if trajfile.endswith(trajFormat)])
        RMSFdict, trajRMSF = cal_rmsf_trajs(args.TOP, trajfiles, mode=args.MODE, selection=args.SEL, outdir=args.OUP)
        print(RMSFdict)
    else:
        RMSFdf = cal_rmsf_traj(args.TOP, args.INP, selection=args.SEL, mode=args.MODE, outfile=args.OUP)
        print(RMSFdf)

if __name__ == "__main__":
    main()


