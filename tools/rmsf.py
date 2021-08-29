
from .base import BaseCalculator

import tqdm
import pandas as pd
import MDAnalysis as mda

'''
Calculat RMSD value for xtc
'''

import mdtraj as md

class RMSFCalculator(BaseCalculator):

    def __init__(self):
        super().__init__()

    def cal_rmsf_xingle_traj(self, topologyfile:str, trajfile:str, selection='mass >= 2', mode='residue') -> pd.DataFrame:
        """Calculate the root mean square fluctuation for the trajfile.

        Args:
            topologyfile (str): topology file for the trajectory
            trajfile (str): molecular dynamic simulation trajectory file.
            selection (str, optional): atom index select to calculate the RMSF. Defaults to 'mass >= 2'.
            mode (str, optional): calculation mode, atom or residue. Defaults to 'residue'.

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
            
        return RMSFdf

    def cal_rmsf_mul_traj(self, topologyfile:str, trajfiles:list, selection='mass>=2', mode='residue') -> pd.DataFrame:
        """Calculate the root mean square fluctuation for multiple trajfile.

        Args:
            topologyfile (str): topology file for the trajectory
            trajfiles (list): molecular dynamic simulation trajectory files.
            selection (str, optional): atom index select to calculate the RMSF.. Defaults to 'mass>=2'.
            mode (str, optional): calculation mode, atom or residue. Defaults to 'residue'.

        Returns:
            pd.DataFrame: [description]
        """
        trajRMSF = {}
        RMSFdict = None
        ProcessBar = tqdm.tqdm(trajfiles)
        for i, trajfile in enumerate(ProcessBar):
            ProcessBar.set_description("Process: %s "%trajfile)
            trajRMSF[i] = self.cal_rmsf_xingle_traj(topologyfile, trajfile, selection, mode)
            if RMSFdict is not None:
                RMSFdict['RMSF'] = RMSFdict['RMSF'].values + trajRMSF[i]['RMSF'].values
            else:
                RMSFdict = trajRMSF[i]
        if RMSFdict is not None:
            RMSFdict['RMSF'] = RMSFdict['RMSF'].values/len(trajfiles)
        return RMSFdict, trajRMSF

    def cal_rmsf_MSM(self, topologyfile:str, stateSamples:str, population:list, selection='mass >2', mode='residue') -> pd.DataFrame:
        """Calculate the root mean square fluctuation for multiple trajfile.

        Args:
            topologyfile (str): topology file for the trajectory
            stateSamples (str): state sample trajectory files generated from makov state model.
            population (list): markov state model population for each state.
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
            stateRMSF[i] = self.cal_rmsf_xingle_traj(topologyfile, stateSample, selection, mode)
            if RMSFdict is not None:
                pi = populations[i]
                RMSFdict['RMSF'] = RMSFdict['RMSF'].values + stateRMSF[i]['RMSF'].values * pi
            else:
                RMSFdict['RMSF'] = stateRMSF[i] * pi
        return RMSFdict, stateRMSF

