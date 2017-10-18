import sys
from ase.lattice.surface import surface
from ase.io import read
from ase.visualize import view
# facet = sys.argv[1]
try:
    layers = int(sys.argv[2])
except:
    layers = 1

try:
    vacuum = float(sys.argv[3])
except:
    vacuum = 9


# facet = [int(i) for i in facet]
facet = [1,1,1]


atoms = read(sys.argv[1])  #Bulk traj file (no vacuum)
atoms = surface(atoms,facet,layers,vacuum)
atoms.set_pbc(True)
view(atoms)

