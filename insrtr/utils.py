"""
Utility functions
"""

resname_3to1_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER':'*',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA':'X', 'MSE': 'M'}

def resname_3to1(resname3: list) -> str:
    """Takes an array of 3 letter names and returns one letter names"""
    #TODO: check for unknown residues
    result = [resname_3to1_dict[resname] for resname in resname3]
    return result
