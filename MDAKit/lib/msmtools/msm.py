

import os

import MDAnalysis as mda
from collections import Counter
from pyemma import coordinates
from pyemma.coordinates import source
class AutoBuildMarkovStateModel(object):
    '''
    Auto build markov state model
    '''
    def __init__(self, topologyfile:str, trajfiles:list, workdir='automsm') -> None:
        super().__init__()
        self.topologyfile = os.path.abspath(topologyfile)
        self.trajfiles = [ os.path.abspath(trajfile) for trajfile in trajfiles ]
        self.workdir = workdir
        if not os.path.exists(workdir):
            os.mkdir(workdir)
        self.src = source(self.trajfiles, top=self.topologyfile)
    
    @property
    def units(self) -> dict:
        """Obtain the units of the trajectories.

        Returns:
            dict: a dict of units
              example:
              {
                  'time'  : 'ps',
                  'length': 'nm'
              }
        """
        if hasattr(self, "_units"):
            units = self._units
        else:
            universe = mda.Universe(self.topologyfile, self.trajfiles[0])
            units = universe.trajectory.units
            self._units = units
        return units

    @property
    def frame_step(self) -> float:
        """Obtain the interval of every frame.

        Returns:
            float: interval step of two frame.
        """
        if hasattr(self, "_dt"):
            dt = self._dt
        else:
            universe = mda.Universe(self.topologyfile, self.trajfiles[0])
            dt = universe.trajectory.dt
            self._dt = dt
        return dt
    
    @property
    def n_traj(self) -> int:
        """Number of trajecotry

        Returns:
            int: the number of trajectory
        """
        return len(self.trajfiles)

    @property
    def n_frames(self) -> int:
        """Number of frames

        Returns:
            int: number of frame
        """
        if hasattr(self, '_n_frame'):
            n_frames = self._n_frames
        else:
            n_frames = 0
            for trajfile in self.trajfiles:
                universe = mda.Universe(self.topologyfile, trajfile)
                trajectory = universe.trajectory
                n_frames += trajectory.n_frames
        self._n_frames = n_frames
        return n_frames

    def traj_stat(self, timeunit='ns') -> None:
        """Statistic trajectory informations

        Args:
            timeunit (str, optional): time units. Defaults to 'ns'.

        Raises:
            Exception: time step for each frame should be same

        """
        UnitDict = {
            'ps' : 1,
            'ns' : 2,
            'us' : 3,
            'ms' : 4
        }
        TotalTime = 0
        if hasattr(self, '_timeunit') and self._timeunit == timeunit:
            TotalTime = self.total_times
            TrajTimes = self.traj_times
        else:
            UnitConvert = 1000**(UnitDict[self.units['time']] - UnitDict[timeunit])
            TrajTimes = []
            for i,trajfile in enumerate(self.trajfiles):
                universe = mda.Universe(self.topologyfile, trajfile)
                trajectory = universe.trajectory
                if trajectory.dt != self.frame_step:
                    raise Exception('The time step should be sample for all the trajectory, please check %s'%trajfile)
                TrajTime = trajectory.totaltime*UnitConvert
                TrajTimes.append(TrajTime)
                TotalTime += TrajTime
            self.traj_times = TrajTimes
            self._timeunit = timeunit
            self.total_times = TotalTime
        counter = Counter(TrajTimes)
        print('Trajectory total time is: %.4f (%s)' % (TotalTime, timeunit))
        print('Trajecotry total number frames: %d ' % self.n_frames)
        print('Trajectory length statistic:')
        print('   TrajLen(%s)\t Count'%timeunit)
        for TrajLen, count in counter.items():
            print('%11.4f\t% 5d' % (TrajLen, count))
            
    def featurelize(self):
         feature = coordinates.featurizer(self.topologyfile)
         feature.add_backbone_torsions()
         self.src.featurizer = feature

        
