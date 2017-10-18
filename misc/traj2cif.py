from ase.all import *
import sys
atoms = read(sys.argv[1])
newname = (sys.argv[1]).split('.')[0]+".cif"
write(newname,atoms)
