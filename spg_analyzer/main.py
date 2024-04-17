import mol_parser
from spg import SPG

# test to check if the function works

# One way to create a molecule
# h2o = Molecule([['O', [0.0, 0.0, 0.0]], ['H', [0.0, 0.0, 1.0]], ['H', [0.0, 1.0, 0.0]]])

# Another way to create a molecule using mol_parser
mol = mol_parser.from_mol('../tests/testset/Ih_2.mol')
# mol = mol_parser.from_mol('../tests/mol_rotated/D8h_rotated.mol')
# mol = mol_parser.from_mol('../tests/mol_samples/C40.mol')

# Build the molecule object
mol.build()

# Create a SPG object from the molecule
mol_spg = SPG(mol)
# mol_spg.mode = 'loose'
# mol_spg.mode = 'very_loose'
# Save the original molecule to .mol file
# mol_parser.to_mol(mol, mol.name + '_original')

# Build the SPG object
mol_spg.build()
# Save the aligned molecule to .mol file
# mol_parser.to_mol(mol_spg.mol, mol.name + '_aligned', '../tests/mol_aligned/')

# print("Coordinates of molecule", mol_spg.mol.coordinates)

# Check the C2 axis
# sym_op.check_C2_perp_Cn(mol_spg.mol.coordinates, mol_spg.Cn_axis)
# sym_op.check_C2_full(mol_spg.mol.coordinates)

# Print the point group of the molecule
print(f"Point group of the molecule {mol_spg.mol.name} is {mol_spg.spg}")
print("Distinguishable symmetry operations are:", mol_spg.so)

