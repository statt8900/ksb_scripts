
"""
This module contains data which forms a backbone for other modules.
The primary data stored here are:
	- domains of interest for all DFT parameters
	- data about traj files to be considered (this could be automated if a designated folder were to be used)
	- Which stoichiometries are known to be magnetic 
	- Ionic radii for use in initial geometry guesses
"""

__all__ = []
#######################################
# RELEVANT DOMAINS 
#######################################

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

pureMetals = bccs+fccs+hcps+diamonds

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

binaryAlloys = rocksalts + cscls + zincblendes

otherStoichs=[] # what structure for PdAu? Need to make a new one for FCC alloy.

stoichStructDomain= pureMetals + binaryAlloys+otherStoichs

__all__+=['bccs','fccs','hcps','diamonds','pureMetals','rocksalts','cscls','zincblendes'
			,'binaryAlloys','otherStoichs','stoichStructDomain']

# Calc Domains
dftDomain     = ['gpaw','quantumEspresso']
xcDomain      = ['PBE','BEEF','mBEEF','RPBE']
pspDomain     = ['gbrv15pbe','sg15']
pwDomain      = [400,700,900,1100,1300,1500,1700]
kptDomain     = [2,4,6] # pts / A^-1
eConvDomain   = [1e-5] # 1e-5 eV
mixDomain     = [0.1] #0.05, #0.1
nMixDomain    = [10]
maxStepDomain = [500]
nBandDomain   = [12]
sigmaDomain   = [0.1]   
fMaxDomain    = [0.05]  # eV/A
xTolDomain    = [0.005] # A
strainDomain  = [0.03]  

__all__ += ['dftDomain','xcDomain','pwDomain','kptDomain','eConvDomain','mixDomain','nMixDomain'
			,'maxStepDomain','nBandDomain','sigmaDomain','xTolDomain','fMaxDomain','strainDomain'
			,'pspDomain']
"""
# Surface Domains
facetDomain   = [(1,1,1),(1,0,0),(0,0,1)]
scaleDomain   = [(1,1,4),(2,2,3),(2,2,4)]
fixedDomain   = [1,2,3]
vacDomain     = [5,7,9]
adsDomain     = [{}]

# And...
clusterDomain = ['sherlock','suncat2']

__all__ += ['facetDomain','scaleDomain','fixedDomain','vacDomain','adsDomain','clusterDomain']
##########################################
# TRAJ DATA
##########################################
"""
#######################################
# MAGNETISM
#######################################

magList = [('Fe2','bcc'),('Ni','fcc'),('Co2','hcp'),('MnO','rocksalt')
			,('MnS','rocksalt'),('MnN','rocksalt'),('MnC','rocksalt')
			,('FeC','rocksalt'),('FeAl','cscl'),('CrC','rocksalt'),('CrN','rocksalt')]
magElems = ['Fe','Mn','Cr','Co','Ni']

__all__+=['magList','magElems']
#######################################
# ATOMIC PROPERTIES
#######################################
radDict   = {'H':0.31	,'He':0.28	,'Li':1.13	,'Be':0.76	,'B':0.84	,'C':1.05,'N':0.71	,'O':0.66	,'F':0.57	,'Ne':0.58	,'Na':1.66	,'Mg':1.41	,'Al':1.26,'Si':1.11	,'P':1.07	,'S':1.05	,'Cl':1.02	,'Ar':1.06	,'K':1.93	,'Ca':1.76
	,'Sc':1.7	,'Ti':1.6	,'V':1.53	,'Cr':1.39	,'Mn':1.553504	,'Fe':1.48417	,'Co':1.2975	,'Ni':1.24	,'Cu':1.22	,'Zn':1.22	,'Ga':1.22	,'Ge':1.2	,'As':1.19	,'Se':1.2	,'Br':1.2	,'Kr':1.16	,'Rb':2.2	,'Sr':1.95	,'Y':1.9	,'Zr':1.75	,'Nb':1.64	,'Mo':1.54	,'Tc':1.47	,'Ru':1.46
	,'Rh':1.42	,'Pd':1.39	,'Ag':1.45	,'Cd':1.44	,'In':1.42	,'Sn':2.29	,'Sb':1.39	,'Te':1.38	,'I':1.39	,'Xe':1.4	,'Cs':2.44	,'Ba':2.15	,'La':2.07	,'Ce':2.04	,'Pr':2.03	,'Nd':2.01	,'Pm':1.99	,'Sm':1.98	,'Eu':1.98	,'Gd':1.96	,'Tb':1.94	,'Dy':1.92	,'Ho':1.92	,'Er':1.89	,'Tm':1.9	,'Yb':1.87	,'Lu':1.87	,'Hf':1.75	,'Ta':1.7	,'W':1.62	,'Re':1.51
	,'Os':1.44	,'Ir':1.41	,'Pt':1.36	,'Au':1.5	,'Hg':1.32	,'Tl':1.45	,'Pb':1.46	,'Bi':1.48	,'Po':1.4	,'At':1.5	,'Rn':1.5	,'Fr':2.6	,'Ra':2.21	,'Ac':2.15	,'Th':2.06	,'Pa':2	,'U':1.96	,'Np':1.9	,'Pu':1.87	,'Am':1.8	,'Cm':1.69}

__all__+=['radDict']

## # of electrons in gbrv15
# gbrv15pbe = PSP('gbrv15pbe','/home/vossj/suncat/psp/gbrv1.5pbe',[1,0,3,4,3,4,5,6,7,0,9,10,3,4,5,6,7,0,9,10,11,12,13,14,15,16,17,18,19,20,19,12,5,6,7,0,9,10,11,12,13,14,15,16,15,16,19,12,13,14,15,6,7,0,9,10,11,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,13,11,15,16,15,16,11,12,13,14,8,0,0,0,0,0,0,0,0,0,0,8,8])
