---
author: Yufeng Tong
date: 2021-01-31
...
# A PyMOL script to calculate and display conformational changes

This is a `PyMOL` script to calculate the r.m.s.d. of two aligned structures with same residue numbering, for comparing the conformations of one protein in two different structures or two proteins (e.g. one wild type and one mutant) in two structures.

The script will create a copy of the target protein, replace the B-factor with the r.m.s.d. and generate a sausage display of the target protein.

## Installation
Save the script to any folder accessible in `PyMOL`.

## Usage
This script has been tested on `PyMOL` version 2.4.

1. Load two structures in the `PyMOL` program. Select polypeptide chains of interest, align them using PyMOL `align` or `super` command, and give them names, say `mol1`, and `mol2`
2. `File` → `Run Script` and load this script.
3. Run `rmsdCA mol1, mol2`, where `mol1` is the reference structure, and `mol2` is the target structure that will be colored.
4. A `tgt_gzt` object will be created and displayed in cartoon-sausage mode with a color ramp.
5. A pdb file `rmsdBFactor_mol2.pdb` will be generated and saved where the B-factor is the r.m.s.d. of the Cα atoms.
6. A csv file `rmsdCA_mol2.csv` will be generated and save where it tabulate the r.m.s.d per residue.


## TODO

