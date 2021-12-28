usage: mdtools [-h] {DSSP,RMSF,Rg,XTC2PDB,PDB2XTC,catTraj,HB,HBComp} ...

MD tools kits

positional arguments:
  {DSSP,RMSF,Rg,XTC2PDB,PDB2XTC,catTraj,HB,HBComp}
    DSSP                Calculate the secondary information.
    RMSF                Calculate the RMSF value.
    Rg                  Calculate the Residue gyration value for ensemble.
    XTC2PDB             Convert single trajectory file to multiple PDB file.
    PDB2XTC             Convert multiple PDB files to single trajectory file.
    catTraj             Read the trajectory file and obtain the information.
    HB                  Calculate the hydrogen bond information.
    HBComp              Compare two msm ensemble hydrogen bond. h1 - h2

optional arguments:
  -h, --help            show this help message and exit
