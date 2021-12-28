'''

Author: maohuay
Email: maohuay@hotmail.com
'''

__all__ = [
    'obtain_sequence', 
    'convert_residue_nonstandard_to_standard',
    'convert_three_letters_to_one_letters',
    'convert_one_letters_to_three_letters',
    'print_sequence',
    'print_sequence_dict'
    ]

import argparse
import logging
from MDAnalysis import Universe
from aminoacids import nonstandardResideuMap, ThreeLettersMap, OneLettersMap


def obtain_sequence(pdbfile:str, mode=3, standard=False) -> dict:
    """Obtain the sequence from a input pdbfile.

    Args:
        pdbfile (str): A input file. support format: pdb, gro, pdbmx
        mode (int, optional): sequence letters only support input 3 or 1. Defaults to 3.
        standard (bool, optional): Convert the resname to standard resname. Defaults to False.

    Raises:
        Exception: The input mode only support a int of 3 or 1

    Returns:
        dict: A dict with sequence list for each chain.
           {
               chainName: sequenceList
           }
    """
    if mode not in [1,3]:
        raise Exception("Only support 1 or 3 for the input of mode. %s"%str(mode))
    u = Universe(pdbfile)
    sequenceDict = {}
    for segment in u.segments:
        ChainID = segment.segid
        segmentSequence = []
        for residue in segment.residues:
            segmentSequence.append(residue.resname)
        if mode == 1:
            segmentSequence = convert_three_letters_to_one_letters(segmentSequence)
        if mode == 3 and standard:
            segmentSequence = convert_residue_nonstandard_to_standard(segmentSequence)
        sequenceDict[ChainID] = segmentSequence
    return sequenceDict

def convert_residue_nonstandard_to_standard(sequenceList: list) -> list:
    """Convert a list of three letter resname to standard protein residue name.

    Args:
        sequenceList (list): A list contains resnames

    Returns:
        list: A list of resname which were converted to standard resname.
    
    Warnings:
        Will the resname not found in the nonstandardResideuMap
    """
    standardSequenceList = []
    for resname in sequenceList:
        if resname not in nonstandardResideuMap:
            logging.warn("Can't convert the resname %s to a standard residue."%resname)
        else:
            resname = nonstandardResideuMap[resname]
        standardSequenceList.append(resname)
    return standardSequenceList

def convert_three_letters_to_one_letters(sequenceList: list) -> list:
    """Convert a list of three letter resname to one letter resname.

    Args:
        sequenceList (list): A list contains resnames

    Returns:
        list: A list of one letter resname.
    """
    OneLettersList = []
    standardSequenceList = convert_residue_nonstandard_to_standard(sequenceList)
    for resname in standardSequenceList:
        if resname not in ThreeLettersMap:
            OneLetter = 'X'
        else:
            OneLetter = ThreeLettersMap[resname]
        OneLettersList.append(OneLetter)
    return OneLettersList

def convert_one_letters_to_three_letters(sequenceList: list) -> list:
    """Convert a list of one letter resname to three letter resname.

    Args:
        sequenceList (list): A list contains one letter resname

    Returns:
        list: A list of three letter resname
    
    Warnings:
        The UNK resname in the result was means a unkonw resname.
    """
    ThreeLettersList = []
    for OneLetter in sequenceList:
        if OneLetter not in OneLettersMap:
            resname = 'UNK'
        else:
            resname = OneLettersMap[OneLetter]
        ThreeLettersList.append(resname)
    return ThreeLettersList

def print_sequence(sequenceList:list, lineLen=80) -> str:
    """Convert a list of sequence resname to a format print.

    Args:
        sequenceList (list): A list of sequence resnames.
        lineLen (int, optional): line length. Defaults to 80.

    Raises:
        Exception: While the length of the first resname in the sequence list
                   was not equation to 1 or 3. 

    Returns:
        str: A str line of format sequence.
    """
    if len(sequenceList)<1:
        return ''
    resLen = len(sequenceList[0])
    if resLen == 3:
        sep = ' '
    elif resLen == 1: 
        sep = ''
    else:
        raise Exception("Found unkonw resname. Resname must be length 3 or 1. %s"%sequenceList[0])
    lines = []
    for i, resname in enumerate(sequenceList):
        if i!=0 and i%lineLen==0:
            lines.append('\n')
        else:
            lines.append(resname)
    linestr = sep.join(lines)
    return linestr

def print_sequence_dict(sequenceDict: dict, lineLen=80):
    """Print a sequence dict like:
        >Chain A length: 85
        GLU THR LEU VAL ARG PRO LYS PRO LEU LEU LEU LYS LEU LEU LYS SER VAL GLY ALA GLN LYS ASP THR TYR THR MET LYS GLU VAL LEU PHE TYR LEU GLY GLN TYR ILE MET THR LYS ARG LEU TYR ASP GLU LYS GLN GLN HIS ILE VAL TYR CYS SER ASN ASP LEU LEU GLY ASP 
        PHE GLY VAL PRO SER PHE SER VAL LYS GLU HIS ARG LYS ILE TYR THR MET ILE TYR ARG ASN LEU VAL VAL
        >Chain B length: 13
        GLU THR PHE SER ASP LEU TRP LYS LEU LEU PRO GLU ASN
        
        >Chain A length: 85
        ETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDEKQQHIVYCSNDLLGD
        FGVPSFSVKEHRKIYTMIYRNLVV
        >Chain B length: 13
        ETFSDLWKLLPEN

    Args:
        dict: A dict with sequence list for each chain.
           {
               chainName: sequenceList
           }
        lineLen (int, optional): line length. Defaults to 80.
    """
    for chainID, sequenceList in sequenceDict.items():
        print('>Chain %s length: %d' % (chainID, len(sequenceList)))
        linestr = print_sequence(sequenceList, lineLen)
        print(linestr)

def main():
    parser = argparse.ArgumentParser(description='Sequence calculation tools.')
    parser.add_argument('-i', dest='PDB', help='A input file. pdb, gro ', required=True)
    parser.add_argument('-m', dest='MODE', help='Resname mode, three letters or one letters.', 
                            choices=[1,3], default=1, type=int)
    parser.add_argument('-n', dest='NORM', help='Convert the resname to standart resname.', 
                            default=True, action='store_false')

    args = parser.parse_args()
    print(args)
    sequenceDict = obtain_sequence(args.PDB, mode=args.MODE, standard=args.NORM)
    print_sequence_dict(sequenceDict, lineLen=60)

if __name__ == "__main__":
    main()