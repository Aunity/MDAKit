
import os
from tools.rg import RgCalculator
from tools.rmsf import RMSFCalculator

from lib.msmtools.msm import AutoBuildMarkovStateModel


#### MSM module test
def test_for_msm():
    topologyfile = './example/pentapeptide/pentapeptide-impl-solv.pdb'
    xtcfiles = sorted([ os.path.join('./example/pentapeptide/', _) for _ in os.listdir('./example/pentapeptide/') if _.endswith('impl-solv.xtc')])
    AutoBuildMSM = AutoBuildMarkovStateModel(topologyfile, xtcfiles, 'automsm')
    print(AutoBuildMSM.units)
    AutoBuildMSM.traj_stat(timeunit='us')
    AutoBuildMSM.traj_stat(timeunit='us')
    print(AutoBuildMSM.n_traj)
    print(AutoBuildMSM.total_times)


#### RMSF module test
def test_for_single():
    RMSFCal = RMSFCalculator()
    #topologyfile = './example/pentapeptide/pentapeptide-impl-solv.pdb'
    #xtcfile = './example/pentapeptide/pentapeptide-00-500ns-impl-solv.xtc'
    topologyfile = './example/1ycr.pdb'
    xtcfile = './example/1ycr.pdb'
    rmsf = RMSFCal.cal_rmsf_xingle_traj(topologyfile, xtcfile)
    print(rmsf)

def test_for_multiple():
    RMSFCal = RMSFCalculator()
    topologyfile = './example/pentapeptide/pentapeptide-impl-solv.pdb'
    xtcfiles = sorted([ os.path.join('./example/pentapeptide/', _) for _ in os.listdir('./example/pentapeptide/') if _.endswith('impl-solv.xtc')])
    rmsf,_ = RMSFCal.cal_rmsf_mul_traj(topologyfile, xtcfiles, selection='name CA', mode='residue')
    print(rmsf)

def test_for_rg():
    topologyfile = './example/1ycr.pdb'
    xtcfile = './example/1ycr.pdb'
    topologyfile = './example/pentapeptide/pentapeptide-impl-solv.pdb'
    xtcfiles = sorted([ os.path.join('./example/pentapeptide/', _) for _ in os.listdir('./example/pentapeptide/') if _.endswith('impl-solv.xtc')])

    RgCal = RgCalculator()
    rg = RgCal.cal_rg_single_traj(topologyfile, xtcfiles[0], selection='all')
    print(rg)

if __name__ == "__main__":
    # test_for_multiple()
    #test_for_msm()
    test_for_rg()