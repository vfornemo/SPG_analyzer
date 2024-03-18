import mol_parser
from spg import Molecule, SPG
import sym_op

# test to check if the function works

# One way to create a molecule
# h2o = Molecule([['O', [0.0, 0.0, 0.0]], ['H', [0.0, 0.0, 1.0]], ['H', [0.0, 1.0, 0.0]]])

# Another way to create a molecule using mol_parser
water = mol_parser.from_mol('../tests/testset/R3.mol')
# water = mol_parser.from_mol('../tests/mol_samples/C40.mol')

# Build the molecule object
water.build()

# Create a SPG object from the molecule
water_spg = SPG(water)
# Save the original molecule to .mol file
# mol_parser.to_mol(water, water.name + '_original')

# Build the SPG object
water_spg.build()
# Save the aligned molecule to .mol file
# mol_parser.to_mol(water_spg.mol, water.name + '_aligned')

# print("Coordinates of molecule", water_spg.mol.coordinates)

# Print the point group of the molecule
print(f"Point group of the molecule {water_spg.mol.name} is {water_spg.spg}")