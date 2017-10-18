from pymatgen import Molecule

"""Database for molecules. Consider merging with miniDatabase"""

CO = Molecule(["C", "O"], [[0, 0, 0],[0, 0, 1.3]])
H  = Molecule(["H"],[[0,0,0]])


molDict = {'CO':CO,'H':H}