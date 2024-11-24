from ase.io import read
from ase.visualize import view

# Load the .xyz file
atoms = read('../op_text/temp.xyz')

# Visualize the structure
view(atoms)
