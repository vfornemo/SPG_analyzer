# This is a script to test the functioning of the spg_analyzer package. 
# It reads the .mol files from the testset folder and checks if the point group of the molecule is correctly identified.

import os, sys
sys.path.append(os.getcwd())
sys.path.append(os.getcwd()+"/spg_analyzer")

from spg_analyzer.mol_parser import from_mol
from spg_analyzer.spg import SPG


# test to check if the function works

test_path = 'tests/testset2/'

f = open('tests/bulk_test.log', 'w')
f.write("Type  ,  Point Group  ,  Result  ,  Pass/Fail\n")


for file in os.listdir(test_path):
    if file.endswith('.mol'):
        mol = from_mol(test_path + file)
        mol.build()
        mol_spg = SPG(mol)
        mol_spg.build()
        if {mol_spg.mol.name} == {mol_spg.spg}:
            pass_test = "Pass"
        else:
            pass_test = "Fail"
        f.write(f"{mol_spg.sym_type}  ,  {mol_spg.mol.name}  ,  {mol_spg.spg}  ,  {pass_test}\n")

