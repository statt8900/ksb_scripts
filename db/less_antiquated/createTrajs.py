import warnings
warnings.filterwarnings("ignore", message="Moved to ase.neighborlist")

from emt  					import EMT
from printParse 			import parseChemicalFormula

from ase 					import io, Atoms
from ase.lattice.cubic 		import FaceCenteredCubic,BodyCenteredCubic,Diamond
from ase.lattice.hexagonal 	import HexagonalClosedPacked

import data_solids_wPBE as data

"""
A systematic way to create Atoms objects to initialize an ASE Database ($SCRATCH/db/ase.db)

This is fairly dependent on the data_solids_wPBE dataset from Keld

Initial magnetic moments are set to 3, but could be changed
"""

##################
# INPUT parameters
##################
magmomInit = 3 # If magnetic, initialize with this magmom

################
# Database Setup
################
from ase.db import connect
db = connect('/scratch/users/ksb/db/ase.db')

####################
#Auxillary Functions
####################
def magLookup(symbs):
	for d in data.data.values():
		if d['symbols'] == symbs: 
			return d['magmom'] is not None
	raise KeyError, 'Symbols {0} were not found in data_solids_wPBE'.format(symbs)

def dataLookup(symbs):
	for d in data.data.values():
		if d['symbols'] == symbs: 
			return d['lattice parameter']
	raise KeyError, 'Symbols {0} were not found in data_solids_wPBE'.format(symbs)

###############
# Mini Database
###############
bccs = 		[('Li','bcc'),   ('Na','bcc'), ('K','bcc'), ('Rb','bcc'),
             ('Ba','bcc'),   ('Mo','bcc'), ('W','bcc'), ('Fe','bcc'),
             ('Nb','bcc'),   ('Ta','bcc')]
			
fccs =      [('Rh','fcc'),   ('Ir','fcc'), ('Al','fcc'),('Pb','fcc'),
             ('Ni','fcc'),   ('Pd','fcc'), ('Pt','fcc'),('Au','fcc'),
			 ('Cu','fcc'),   ('Ag','fcc'), ('Ca','fcc'),('Sr','fcc')]

hcps = 		[('Cd','hcp'),   ('Co','hcp'), ('Os','hcp'), ('Ru','hcp'), 
			 ('Zn','hcp'),   ('Zr','hcp'), ('Sc','hcp'), ('Be','hcp'), 
			 ('Mg','hcp'),   ('Ti','hcp')]

diamonds =  [('Sn','diamond'),('C','diamond'),('Si','diamond'),('Ge','diamond')]

rocksalts = [	('LiH','rocksalt'),	('LiF','rocksalt'),
				('LiCl','rocksalt'),('NaF','rocksalt'),
				('NaCl','rocksalt'),('MgO','rocksalt'),
				('MgS','rocksalt'),	('CaO','rocksalt'),
				('TiC','rocksalt'),	('TiN','rocksalt'),
				('ZrC','rocksalt'),	('ZrN','rocksalt'),
				('VC','rocksalt'),	('VN','rocksalt'),
				('NbC','rocksalt'),	('NbN','rocksalt'),
				('KBr','rocksalt'),	('CaSe','rocksalt'),
				('SeAs','rocksalt'),('RbI','rocksalt'),
				('LiI','rocksalt'),	('CsF','rocksalt'),
				('AgF','rocksalt'),	('AgCl','rocksalt'),
				('AgBr','rocksalt'),('CaS','rocksalt'),
				('BaO','rocksalt'),	('BaSe','rocksalt'),
				('CdO','rocksalt'),	('MnO','rocksalt'),
				('MnS','rocksalt'),	('ScC','rocksalt'),
				('CrC','rocksalt'),	('MnC','rocksalt'),
				('FeC','rocksalt'),	('CoC','rocksalt'),
				('NiC','rocksalt'),	('ScN','rocksalt'),
				('CrN','rocksalt'),	('MnN','rocksalt'),
				('FeN','rocksalt'),	('CoN','rocksalt'),
				('NiN','rocksalt'),	('MoC','rocksalt'),
				('RuC','rocksalt'),	('RhC','rocksalt'),
				('PdC','rocksalt'),	('MoN','rocksalt'),
				('RuN','rocksalt'),	('RhN','rocksalt'),
				('PdN','rocksalt'),	('LaC','rocksalt'),
				('TaC','rocksalt'),	('WC','rocksalt'),
				('OsC','rocksalt'),	('IrC','rocksalt'),
				('PtC','rocksalt'),	('LaN','rocksalt'),
				('TaN','rocksalt'),	('WN','rocksalt'),
				('OsN','rocksalt'),	('IrN','rocksalt'),
				('PtN','rocksalt')]

cscls = 		[('CsI','cscl'),	    ('NiAl','cscl'),
				('FeAl','cscl'),	    ('CoAl','cscl')]

zincblendes	=	[('BN','zincblende'),		('BP','zincblende'),
				('BAs','zincblende'),		('AlN','zincblende'),
				('AlP','zincblende'),		('AlAs','zincblende'),
				('GaN','zincblende'),		('GaP','zincblende'),
				('GaAs','zincblende'),		('InP','zincblende'),
				('InAs','zincblende'),		('SiC','zincblende')]

magElems = ['Fe','Mn','Cr','Co','Ni']

#######################################
# ATOMIC PROPERTIES
#######################################
radDict   = {'Sn':4.3}

###############
# Main Function
###############

if __name__ == '__main__':
	main()

def main():
	for x in fccs+bccs+hcps+diamonds+rocksalts+zincblendes+cscls:
		f = x in fccs
		b = x in bccs
		h = x in hcps
		d = x in diamonds
		r = x in rocksalts
		z = x in zincblendes
		c = x in cscls
		
		try: 
			if   f: a = FaceCenteredCubic(x[0]).get_cell()[0][0]/2
			elif b: a = BodyCenteredCubic(x[0]).get_cell()[0][0]/2
			elif d: a = Diamond(x[0]).get_cell()[0][0]/2
			elif h: 
					cell = HexagonalClosedPacked(x[0]).get_cell()
					a,c = cell[0][0],cell[2][2]

			elif r | z | c: a = dataLookup(x[0])
			else: raise NotImplementedError

		except ValueError: 
			a = sum([radDict[e] for e in elems])/len(elems)
			print "Had to guess lattice constant of "+x[0]

		if   f: name,struc,pos,cell,n = '-fcc',       'fcc',        [[0,0,0]],                  [[0,a,a],[a,0,a],[a,a,0]],   1
		elif b: name,struc,pos,cell,n = '-bcc',       'bcc',        [[0,0,0]],                  [[a,a,-a],[-a,a,a],[a,-a,a]],1
		elif h: name,struc,pos,cell,n = '-hcp',       'hexagonal',  [[0,0,0],[2./3,1./3,1./2]], [a,a,c,90,90,120],           2
		elif d: name,struc,pos,cell,n = '-diamond',   'diamond',    [[0,0,0],[0.25,0.25,0.25]], [[0,a,a],[a,0,a],[a,a,0]],   2
		elif z: name,struc,pos,cell,n = '-zincblende','zincblende', [[0,0,0],[0.25,0.25,0.25]], [[0,a,a],[a,0,a],[a,a,0]],   1
		elif r: name,struc,pos,cell,n = '-rocksalt',  'rocksalt',   [[0,0,0],[0.5,0.5,0.5]],    [[0,a,a],[a,0,a],[a,a,0]],   1
		elif c: name,struc,pos,cell,n = '-cscl',      'cubic',      [[0,0,0],[0.5,0.5,0.5]],    [a,a,a,90,90,90],            1

		mag     = magLookup(x[0])
		elems 	= parseChemicalFormula(x[0]).keys()*n
		magmoms = [magmomInit if e in magElems else 0 for e in elems] 

		atoms = Atoms(elems,scaled_positions=pos,cell=cell,pbc=[1,1,1],magmoms=magmoms,calculator=EMT())
		
		info = 	{'name': 		x[0]+name
				,'relaxed': 	False
				,'emt': 		atoms.get_potential_energy()/len(elems)	#normalized to per-atom basis
				,'comments':	'Autogenerated by createTrajs'
				,'kind': 		'bulk' # vs surface/molecules
				### Stuff for bulk
				,'structure': 	struc
				### Stuff for surfaces
				,'parent': 		None 
				,'sites': 		None
				,'facet': 		None
				,'xy': 			None
				,'layers': 		None
				,'constrained': None
				,'symmetric': 	None
				,'vacuum': 		None
				,'adsorbates':	None}

		db.write(atoms,key_value_pairs=info)







