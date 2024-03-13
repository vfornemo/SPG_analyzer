import mol_parser
from spg import Molecule, SPG
from utils import sort_coords, thre_cut
import sym_op
import numpy as np

# test to check if the function works

# One way to create a molecule
# h2o = Molecule([['O', [0.0, 0.0, 0.0]], ['H', [0.0, 0.0, 1.0]], ['H', [0.0, 1.0, 0.0]]])

# Another way to create a molecule using mol_parser
water = mol_parser.from_mol('../tests/testset/12_D6h.mol')
# water = mol_parser.from_mol('../tests/mol_samples/H2O2.mol')

# Build the molecule after creating it
water.build()

# Create a SPG object from the molecule
water_spg = SPG(water)
mol_parser.to_mol(water, water.name + '_original')

# Build the SPG object after creating it
water_spg.build()
mol_parser.to_mol(water_spg.mol, water.name + '_aligned')

# Check rotation

print(water_spg.mol.coordinates)

rot_z = thre_cut(sym_op.Cn_rot(water_spg.mol.coordinates, 2, 'z'))
rot_y = thre_cut(sym_op.Cn_rot(water_spg.mol.coordinates, 2, 'y'))
rot_x = thre_cut(sym_op.Cn_rot(water_spg.mol.coordinates, 2, 'x'))


print("rotation_z", rot_z, '\n', sort_coords(rot_z))
print("rotation_y", rot_y, '\n', sort_coords(rot_y))
print("rotation_x", rot_x, '\n', sort_coords(rot_x))



if np.array_equal(sort_coords(rot_z), sort_coords(water_spg.mol.coordinates)):
    print("C2 rotation about z-axis is a symmetry operation for the molecule.")
elif np.array_equal(sort_coords(rot_y), sort_coords(water_spg.mol.coordinates)):
    print("C2 rotation about y-axis is a symmetry operation for the molecule.")
elif np.array_equal(sort_coords(rot_x), sort_coords(water_spg.mol.coordinates)):
    print("C2 rotation about x-axis is a symmetry operation for the molecule.")
else:
    print("C2 rotation is not a symmetry operation for the molecule.")
