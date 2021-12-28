import os
import sys
import argparse
# import argcomplete
import readline

CMD = ['DSSP', 'RMSF', 'Rg', 'XTC2PDB', 'PDB2XTC', 'catTraj']
def completer(text, state):
    options = [cmd for cmd in CMD if cmd.startswith(text)]
    if state < len(options):
        return options[state]
    else:
        return None



def parse_dssp(args):
    from .tools.dssp import cal_dssp_PDB, cal_dssp_traj
    if args.TRAJ:
        cal_dssp_traj(args.TRAJ, args.PDB, args.MODE)
    else:
        cal_dssp_PDB(args.PDB, args.MODE)

def parse_rmsf(args):
    from .tools.rmsf import cal_rmsf_traj, cal_rmsf_trajs
    trajFormat = '.xtc'
    if os.path.isdir(args.INP):
        trajfiles = sorted([os.path.join(args.INP, trajfile) for trajfile in os.path.listdir(args.INP) if trajfile.endswith(trajFormat)])
        RMSFdict, trajRMSF = cal_rmsf_trajs(args.TOP, trajfiles, mode=args.MODE, selection=args.SEL, outdir=args.OUP)
    else:
        RMSFdf = cal_rmsf_traj(args.TOP, args.INP, selection=args.SEL, mode=args.MODE, outfile=args.OUP)

def parse_rg(args):
    from .tools.rg import cal_rg_PDB, cal_rg_trajs, cal_rg_traj
    if args.TOP is None:
        cal_rg_PDB(args.INP, selection=args.SEL)
    elif os.path.isdir(args.INP):
        trajFormat = '.xtc'
        trajfiles = sorted([os.path.join(args.INP, trajfile) for trajfile in os.listdir(args.INP) if trajfile.endswith(trajFormat)])
        RgDict = cal_rg_trajs(args.TOP, trajfiles, selection=args.SEL, outdir=args.OUP)
        print(RgDict)
    else:
        Rgs = cal_rg_traj(args.TOP, args.INP, selection=args.SEL, outfile=args.OUP)
        print(Rgs)

def parse_xtc2pdb(args):
    from .tools.convert import convert_traj_to_pdbs
    convert_traj_to_pdbs(args.INP, args.TOP, args.SEL, args.OUP)

def parse_pdb2xtc(args):
    from .tools.convert import convert_pdbs_to_traj
    convert_pdbs_to_traj(args.INP, args.SEL, args.OUP)

def parse_cattraj(args):
    from .tools.catTraj import cat_traj
    cat_traj(args.INP, args.TOP)

def parse_hb(args):
    from .tools.hydrogenbond import cal_hb_PDB, cal_hb_traj, cal_hb_trajs
    if os.path.isdir(args.INP):
        cal_hb_trajs(args.INP, args.TOP, args.OUP)
    elif os.path.isfile(args.INP) and args.TOP:
        cal_hb_traj(args.INP, args.TOP, args.OUP)
    else:
        cal_hb_PDB(args.INP, args.OUP)

def parse_hbcomp(args):
    from .tools.hydrogenbondCompare import compare_hydrogenbonds
    compare_hydrogenbonds(args.h1f, args.h2f, args.outf)

def parse_args():
    readline.parse_and_bind("tab: complete")
    readline.set_completer(completer)

    MainParser = argparse.ArgumentParser(description='MD tools kits', prog='mdtools')
    subparsers = MainParser.add_subparsers()

    # dssp
    dsspParser = subparsers.add_parser('DSSP', help='Calculate the secondary information.')
    dsspParser.add_argument('-i', dest='PDB', help='The input PDB file.', required=True)
    dsspParser.add_argument('-t', dest='TRAJ', help='Trajectory file to calculate the ssp.', default=None)
    dsspParser.add_argument('-m', dest='MODE', help='Simple mode or complex mode.', default=3, choices=[3,7], type=int)
    dsspParser.set_defaults(func=parse_dssp)

    # rmsf
    rmsfParser = subparsers.add_parser('RMSF', help='Calculate the RMSF value.')
    rmsfParser.add_argument('-i', dest='INP', help='Input trajectory file or a floder contains trajectory files.', required=True)
    rmsfParser.add_argument('-t', dest='TOP', help='Topology file for the input trajectory file.', required=True)
    rmsfParser.add_argument('-s', dest='SEL', help='Selection atoms for calculate the RMSF.', default='mass >2')
    rmsfParser.add_argument('-m', dest='MODE', help='Calculation mode: residue or atom, default is residue.', choices=['residue','atom'], default='residue')
    rmsfParser.add_argument('-o', dest='OUP', help='Output file or output directory.', default=None)
    rmsfParser.set_defaults(func=parse_rmsf)

    # Rg
    rgParser = subparsers.add_parser('Rg', help='Calculate the Residue gyration value for ensemble.')
    rgParser.add_argument('-i', dest='INP', help='A simple pdbfile or a trajectory file or a floder contains trajectory files.', required=True)
    rgParser.add_argument('-t', dest='TOP', help='Topology file for the input trajectory file.')
    rgParser.add_argument('-s', dest='SEL', help='Selection atoms for calculate the RMSF.', default='mass >2')
    rgParser.add_argument('-o', dest='OUP', help='Output file or output directory.', default=None)
    rgParser.set_defaults(func=parse_rg)

    # xtc2pdb
    xtc2pdbParser = subparsers.add_parser('XTC2PDB',help="Convert single trajectory file to multiple PDB file.")
    xtc2pdbParser.add_argument('-i', dest='INP', help='Input trajectory file.', required=True)
    xtc2pdbParser.add_argument('-t', dest='TOP', help='Topology file for the trajectory file.', required=True)
    xtc2pdbParser.add_argument('-s', dest='SEL', help='selection str for group select.', default='all')
    xtc2pdbParser.add_argument('-o', dest='OUP', help='Output file or output floder.', default=None)
    xtc2pdbParser.set_defaults(func=parse_xtc2pdb)

    # pdb2xtc
    pdb2xtcParser = subparsers.add_parser('PDB2XTC', help="Convert multiple PDB files to single trajectory file.")
    pdb2xtcParser.add_argument('-i', dest='INP', help='Input PDB file floder.', required=True)
    pdb2xtcParser.add_argument('-s', dest='SEL', help='selection str for group select.', default='all')
    pdb2xtcParser.add_argument('-o', dest='OUP', help='Output file or output floder.', default=None)
    pdb2xtcParser.set_defaults(func=parse_pdb2xtc)

    # catTraj
    cattrajParser = subparsers.add_parser('catTraj', help='Read the trajectory file and obtain the information.')
    cattrajParser.add_argument('-i', dest='INP', help='Trajectory file to cat.', required=True)
    cattrajParser.add_argument('-t', dest='TOP', help='Topology file for the trajectory file.', required=True)
    cattrajParser.set_defaults(func=parse_cattraj)
    
    # HB
    hbParser = subparsers.add_parser('HB', help="Calculate the hydrogen bond information.")
    hbParser.add_argument('-i', dest='INP', help='Input PDB file or Input trajectory file or Input floder contains trajectories.', required=True)
    hbParser.add_argument('-t', dest='TOP', help='Topology file for trajectory file, need for INP is trajectory.', default=None)
    hbParser.add_argument('-o', dest='OUP', help='Output directory.', default=None)
    hbParser.set_defaults(func=parse_hb)

    #HBComp
    hbCompParser = subparsers.add_parser('HBComp', help="Compare two msm ensemble hydrogen bond. h1 - h2")
    hbCompParser.add_argument("-h1", dest="h1f", help="The ensemble hydrogen bond file.",required=True)
    hbCompParser.add_argument("-h2", dest="h2f", help="As the h1",required=True)
    hbCompParser.add_argument("-o", dest="outf", help="The result file name.", default="hb-msm-comp.txt")
    hbCompParser.set_defaults(func=parse_hbcomp)

    
    
    # argcomplete.autocomplete(MainParser)
    args = MainParser.parse_args()
    if not hasattr(args, 'func'):
        args = MainParser.parse_args(['-h'])
    args.func(args)
    #print(args)


def main():
    parse_args()

if __name__ == "__main__":
    main()