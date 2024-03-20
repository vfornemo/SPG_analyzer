# This is a script to rotate molecule for test

import os, sys
sys.path.append(os.getcwd())
sys.path.append(os.getcwd()+"/spg_analyzer")

# import numpy as np
from spg_analyzer.utils import rotate_molecule
from spg_analyzer.mol_parser import from_mol, to_mol

mol = from_mol('tests/testset/D8h.mol')
# mol = mol_parser.from_mol('../tests/mol_samples/C40.mol')

# Build the molecule object
mol.build()

# Rotate the molecule
mol.coordinates = rotate_molecule(mol.coordinates)

# Save the rotated molecule to .mol file
to_mol(mol, mol.name + '_rotated', 'tests/mol_rotated/')


