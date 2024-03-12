import numpy as np
from pyscf import gto
from pyscf.symm import *

# input molecule in pyscf format
mol = gto.Mole()
mol.atom = '''
C    1.3940    0.0000    0.0000 
C    0.6970   -1.2072    0.0000 
C   -0.6970   -1.2072    0.0000 
C   -1.3940    0.0000    0.0000 
C   -0.6970    1.2072    0.0000 
C    0.6970    1.2072    0.0000 
H    2.4784    0.0000    0.0000 
H    1.2392   -2.1464    0.0000 
H   -1.2392   -2.1464    0.0000 
H   -2.4784    0.0000    0.0000 
H   -1.2392    2.1464    0.0000 
H    1.2392    2.1464    0.0000 
'''
mol.basis = 'sto-3g'
mol.build()

atom = [['C', (1.3940, 0.0000, 0.0000)],
    ['C', (0.6970, -1.2072, 0.0000)],
    ['C', (-0.6970, -1.2072, 0.0000)],
    ['C', (-1.3940, 0.0000, 0.0000)],
    ['C', (-0.6970, 1.2072, 0.0000)],
    ['C', (0.6970, 1.2072, 0.0000)],
    ['H', (2.4784, 0.0000, 0.0000)],
    ['H', (1.2392, -2.1464, 0.0000)],
    ['H', (-1.2392, -2.1464, 0.0000)],
    ['H', (-2.4784, 0.0000, 0.0000)],
    ['H', (-1.2392, 2.1464, 0.0000)],
    ['H', (1.2392, 2.1464, 0.0000)]]

gpname, orig, axes = detect_symm(atom)
print(gpname, orig, axes)
atom = shift_atom(atom, orig, axes)
print(gpname, symm_identical_atoms(subgroup(gpname, axes)[0], atom))

a = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])	
b = np.array([[4,7,9]])
# dot
print(np.dot(a, b.T))
print(np.dot(b, a))
