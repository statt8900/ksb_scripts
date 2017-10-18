import sys
import numpy as np
from ase import io
from copy import deepcopy

"""Arguments:
	- filename of traj file to read
	- filename of output cif file
	- approximate x dimension (nm)
	- approximate y dimension (nm)
	- approximate z dimension (nm)

Script identifies which layer each atom is in via the "tag" field
This can be modified in ase gui (highlight atoms, ctrl+y, then set the tag)
The tag convention for surfaces is that adsorbates are 0, first layer is 1, etc.
This script just requires that fixed layers are tagged
"""
# Process Input
x,y,z = float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]) #dimensions, nm
atoms = io.read(sys.argv[1])
cell  = atoms.get_cell()
dx    = np.linalg.norm(cell[0])
dy    = np.linalg.norm(cell[1])

# Identify tags and layers
maxTag       = max([atom.tag for atom in atoms])
constrained  = [deepcopy(atoms[i]) for i in atoms.constraints[0].get_indices()]
nFixedLayers = len(set([a.tag for a in constrained]))
# Identify layer spacing
bottomLayerHeight = min([atom.z for atom in atoms if atom.tag == maxTag])
nextToBottomHeight = min([atom.z for atom in atoms if atom.tag == (maxTag - 1)])
spacing = nFixedLayers * (nextToBottomHeight - bottomLayerHeight)
if spacing < 1: print "Don't forget to tag fixed layers!"
# Get number of replications needed to exceed specified dimensions
nx = int(x / (dx/10)) # how many dx's (converted to nm) can fit in x?
ny = int(y / (dy/10))
nz = int(z / (spacing/10))
# Add atoms to get desired thickness
for i in range(nz):
	atoms.set_cell(atoms.get_cell()+[[0,0,0],[0,0,0],[0,0,spacing]])
	for atom in atoms:
		atom.z+=spacing
	for a in constrained: atoms.append(a)
# Scale x,y directions and write output
atoms = atoms*(nx,ny,1)
io.write(sys.argv[2],atoms)

if True:
	eleDict = {'Cu':1,'O':2}
	with open('feff.inp','w') as f:
		f.write(""" POTENTIALS
	 *   ipot   Z  element        l_scmt  l_fms   stoichiometry
	     0     29     Cu           3       3       0.01
	     1     29     Cu           3       3       0.98
	     1     8     O            3       3       0.01\n\n""")

		for a in atoms:
			pos = a.position
			f.write('\n{0} {1} {2} {3}'.format(pos[0],pos[1],pos[2],eleDict[a.symbol]))
			


		
	
	