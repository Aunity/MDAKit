
from .base import BaseCalculator
import numpy as np
import MDAnalysis as mda

class RgCalculator(BaseCalculator):

    def __init__(self) -> None:
        super().__init__()

    def cal_rg_single_traj(self, topologyfile:str, trajfile:str, selection='all') -> np.array:

        universe = mda.Universe(trajfile, topologyfile)
        rgs = []
        for ts in universe.trajectory:
            atoms = universe.select_atoms(selection)
            rgs.append(atoms.radius_of_gyration())
        rgs = np.array(rgs)

        return rgs
        

