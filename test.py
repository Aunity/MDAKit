
import os
from tools.rmsf import RMSFCalculator


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
if __name__ == "__main__":
    test_for_multiple()