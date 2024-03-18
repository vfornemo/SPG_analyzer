"""
This module contains functions to parse .mol files and convert them to Molecule objects. 
Also contains functions to convert Molecule objects to .xyz and .mol files.

Functions:
- `from_mol(mol)`: Convert a .mol V2000 file to a Molecule object.
- `to_xyz(mole, filename)`: Convert a Molecule object to a .xyz file.
- `to_mol(mole, filename)`: Convert a Molecule object to a .mol file.

"""


from genericpath import exists
from spg import Molecule
'''
example mol file

CT1000292221

Created by GaussView 6.0.16
  3  2  0  0  0  0  0  0  0  0  0    0
    0.0021   -0.0041    0.0020 H   0  0  0  0  0  0  0  0  0  0  0  0    # atom block must be aligned perfectly,
   -0.0110    0.9628    0.0073 O   0  0  0  0  0  0  0  0  0  0  0  0    # and the space between coords must not be less than 3, don't know why
    0.8669    1.3681    0.0011 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END

'''

def from_mol(molfile):
    """
    Convert a .mol V2000 file to a Molecule object.

    Args:
        mol: .mol file.

    Returns:
        Molecule object 
    """
    # .mol file contains several blocks, first three lines are header, then atom block, then bond block, then properties block.
    # We only need the atom block.
    
    if exists(molfile) == False:
        raise FileNotFoundError(f'File does not exist in {molfile}.')
    
    f = open(molfile, 'r')
    lines = f.readlines()
    # read atom block
    natm = int(lines[3].split()[0])
    atm_coord = []
    atm_name = []
    for i in range(4, 4 + natm):
        line = lines[i].split()
        atm_coord.append([float(line[0]), float(line[1]), float(line[2])])
        atm_name.append(line[3])
    f.close()

    mole = [[atm_name[i], atm_coord[i]] for i in range(natm)]
    molec = Molecule(mole)
    molec.name = molfile.split('/')[-1].split('.')[0]
    molec.headblock = lines[:4]
    molec.atmblock = lines[4:4 + natm]
    molec.atmblock = [x.split()[3:] for x in molec.atmblock]
    molec.bondblocktoend = lines[4 + natm:]
    return molec


def to_xyz(mole, filename):
    """
    Convert a Molecule object to a .xyz file.

    Args:
        mole: Molecule object.
        filename: name of the .xyz file.

    Returns:
        .xyz file.
    """
    natm = mole.natm
    f = open('../tests/' + filename + '.xyz', 'w')
    f.write(f'{natm}\n\n')
    for i in range(natm):
        f.write(f'{mole.atm_name[i]} {mole.coordinates[i][0]} {mole.coordinates[i][1]} {mole.coordinates[i][2]}\n')
    f.close()
    return

def to_mol(mole, filename):
    """
    Convert a Molecule object to a .mol file.

    Args:
        mole: Molecule object.
        filename: name of the .mol file.

    Returns:
        .mol file.
    """
    f = open('../tests/' + filename + '.mol', 'w')
    if not hasattr(mole, 'headblock'):
        RuntimeError('Molecule object does not have .mol file information.')
    f.writelines(mole.headblock)
    # mole.atm_block_updater()
    for i in range(mole.natm):
        f.write(f'   {mole.coordinates[i][0]: .4f}   {mole.coordinates[i][1]: .4f}   {mole.coordinates[i][2]: .4f}'+ ' ' +'  '.join(mole.atmblock[i])+'\n')
    f.writelines(mole.bondblocktoend)
    f.close()
    return


# test to check if the function works
def test_from_mol():
    mol = '../tests/water.mol'
    mole = from_mol(mol)
    print(mole)
    return

if __name__ == "__main__":
    test_from_mol()

