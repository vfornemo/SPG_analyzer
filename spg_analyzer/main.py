import mol_parser
from spg import Molecule, SPG

# test to check if the function works

# One way to create a molecule
# h2o = Molecule([['O', [0.0, 0.0, 0.0]], ['H', [0.0, 0.0, 1.0]], ['H', [0.0, 1.0, 0.0]]])

# Another way to create a molecule using mol_parser
water = mol_parser.from_mol('../tests/testset/4_Td.mol')
# water = mol_parser.from_mol('../tests/mol_samples/H2O2.mol')

# Build the molecule after creating it
water.build()

# Create a SPG object from the molecule
water_spg = SPG(water)
mol_parser.to_mol(water, water.name + '_original')

# Build the SPG object after creating it
water_spg.build()
mol_parser.to_mol(water_spg.mol, water.name + '_aligned')
