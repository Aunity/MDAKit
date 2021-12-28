
from .base import BaseCalculator

import tqdm
import numpy as np
import pandas as pd

'''
Calculat RMSD value for xtc
'''

import mdtraj as md
import MDAnalysis as mda


class RMSDCalculator(BaseCalculator):

    def __init__(self):
        super().__init__()

    def cal_rmsd_xingle_traj(self, topologyfile:str, trajfile:str, reffile=None, selection='mass >= 2') -> np.array:
        """Calculate the root mean square deviation for the trajfile.

        Args:
            topologyfile (str): topology file for the trajectory
            trajfile (str): molecular dynamic simulation trajectory file.
            reffile ([type], optional): [description]. Defaults to None.
            selection (str, optional): atom index select to calculate the RMSF. Defaults to 'mass >= 2'.


        Returns:
            np.array: RMSD value with data frame
        """
        traj = md.load(trajfile, top=topologyfile)
        if reffile is not None:
            ref = md.load(reffile)
        else:
            ref = traj[0]
        AtomIndex = traj.topology.select(selection)
        # target, reference, frame of reference
        rmsd = md.rmsd(traj, ref, atom_indices=AtomIndex)
            
        return rmsd

    def cal_rmsd_mul_traj(self, topologyfile:str, trajfiles:list, reffile=None, selection='mass>=2') -> pd.DataFrame:
        """Calculate the root mean square fluctuation for multiple trajfile.

        Args:
            topologyfile (str): topology file for the trajectory
            trajfiles (list): molecular dynamic simulation trajectory files.
            selection (str, optional): atom index select to calculate the RMSF.. Defaults to 'mass>=2'.
            mode (str, optional): calculation mode, atom or residue. Defaults to 'residue'.

        Returns:
            pd.DataFrame: [description]
        """
        trajRMSD = {}
        ProcessBar = tqdm.tqdm(trajfiles)
        for i, trajfile in enumerate(ProcessBar):
            ProcessBar.set_description("Process: %s "%trajfile)
            trajRMSD[i] = self.cal_rmsd_xingle_traj(topologyfile, trajfile, reffile, selection)
        return trajRMSD

    def cal_rmsd_ensemble(self, topologyfile:str, stateSamples:str, reffile:str, populations:list, selection='mass >2', mode='residue') -> pd.DataFrame:
        """Calculate the root mean square fluctuation for multiple trajfile.

        Args:
            topologyfile (str): topology file for the trajectory
            stateSamples (str): state sample trajectory files generated from makov state model.
            populations (list): markov state model population for each state.
            selection (str, optional): atom index select to calculate the RMSD. Defaults to 'mass>=2'.
            mode (str, optional): calculation mode, atom or residue. Defaults to 'residue'.

        Returns:
            pd.DataFrame: [description]
        """
        stateRMSD = {}
        RMSDList = []
        weights = []
        ProccessBar = tqdm.tqdm(stateSamples)
        for i, stateSample in enumerate(ProccessBar):
            ProccessBar.set_description("Process: %s " % stateSample)
            stateRMSD[i] = self.cal_rmsd_xingle_traj(topologyfile, stateSample, reffile=reffile, selection=selection)
            RMSDList.append(stateRMSD[i])
            weights.append([populations[i]]*len(stateRMSD[i]))
        RMSDList = np.concatenate(RMSDList)
        weights = np.array(weights)
        return RMSDList, populations, stateRMSD

