
import numpy as np
#! TODO

class BaseCalculator(object):

    def __init__(self) -> None:
        pass

    def cal_msm(self) -> None:
        """calculate for msm model
        """

    def cal_xingle_traj(self) -> None:
        """calculate for single trjectory file
        """
    def cal_mul_traj(self) -> None:
        """calculate for multi trjectory file
        """
    def hist_1D(self, dat: np.array, bins=30, range=None, normed=None, weights=None):
        """Compute the histogram of a dataset.

        Args:
            dat (np.array): Input data. The histogram is computed over the flattened array.
            bins (int, optional): [description]. Defaults to 30.
            range ([type], optional): [description]. Defaults to None.
            normed ([type], optional): [description]. Defaults to None.
            weights ([type], optional): [description]. Defaults to None.
        """
        hist, edges = np.histogram(dat, bins=bins, range=range, normed=normed, weights=weights)
        return hist, edges

