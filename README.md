# SPG_analyzer
Analyze the symmetry point group (SPG) of a given molecule. 

## Usage
```
from spg_analyzer import mol_parser
from spg_analyzer.spg import SPG

# set path to .mol file
mol = mol_parser.from_mol('../tests/testset/Ih.mol')

# Build the molecule object
mol.build()

# Create a SPG object from the molecule
mol_spg = SPG(mol)

# Build the SPG object
mol_spg.build()

# Print the point group of the molecule
print(f"Point group of the molecule {mol_spg.mol.name} is {mol_spg.spg}")

```

## Parameters

### Precision mode
To adjust precision mode, use:
```
mol_spg.mode = "loose"
```
The default settings for precision mode is "medium". Other modes are `"very_tight"`, `"tight"`, `"very_loose"`, `"super_loose"`, `"ultra_loose"`. 

Warning: excessively loose or tight precision brings inaccuracy. Using default 
precision is recommended.

### Aligning
Auto aligning of the molecule may cause problems, especially when the molecule 
is already aligned perfectly. To turn off aligning, use: 
```
mol_spg.align = False
```