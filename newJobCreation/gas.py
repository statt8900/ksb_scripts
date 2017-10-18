from pymatgen 						import Molecule
from pymatgen.io.ase 				import AseAtomsAdaptor

##############
# Adsorbates #
##############
class Adsorbate(object):
	def __init__(self,pmg,vector): 
		self.name  	= name
		self.pmg  	= pmg #should be centered around 0,0,0
		self.vector = vector
		self.offset = offset

H2O = Molecule(['O','H','H'],[[.000000, .000000, .119262]
							,[.000000, .763239,-.477047]
							,[.000000,-.763239,-.477047]])

O2  = Molecule(["O", "O"], 	[[0, 0, 0]
							,[0, 0, 1.25]])

CH4 = Molecule(['C','H','H','H','H'],	[[0,0,0]
										,[ .629118, .629118, .629118]
										,[ -.629118,-.629118, .629118]
										,[.629118,-.629118,-.629118]
										,[ -.629118, .629118,-.629118]])

CO  = Molecule(["C", "O"], 	[[0, 0, 0]
							,[0, 0, 1.2]])

CO2 = Molecule(["C", "O",'O'], 	[[0, 0, 0]
								,[0, 0, 1.15]
								,[0, 0, -1.15]])

H   = Molecule(["H"],[[0,0,0]])

H2  = Molecule(["H",'H'],	[[0,0,0]
							,[0,0,0.75]])

N2  = Molecule(["N",'N'],	[[0,0,0]
							,[0,0,1.1]])

Cl2 = Molecule(["Cl",'Cl'],	[[0,0,0]
							,[0,0,2]])

Br2 = Molecule(["Br",'Br'],	[[0,0,0]
							,[0,0,2.3]])

F2  = Molecule(["F",'F'],[[0,0,0]
						,[0,0,1.4]])

molDict = 	{'CO':CO
			,'H':H
			,'H2':H2
			,'Cl2':Cl2
			,'N2':N2
			,'Br2':Br2
			,'F2':F2
			,'H2O':H2O
			,'O2':O2
			,'CH4':CH4
			,'CO2':CO2}

aseMolDict = {k:AseAtomsAdaptor.get_atoms(v.get_boxed_structure(20,20,20)) for k,v in molDict.items()}
