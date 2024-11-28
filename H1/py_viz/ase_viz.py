from ase.io import read
from ase.visualize import view

# Load the .xyz file
atoms = read('../op_text/equilibration_pos.xyz')

# Visualize the structure
view(atoms)
