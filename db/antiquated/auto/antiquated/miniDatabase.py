from traj import TrajDataBulk, TrajDataSurf

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

bccs = 		[('Li2','bcc'),   ('Na2','bcc'), ('K2','bcc'), ('Rb2','bcc'),
             ('Ba2','bcc'),   ('Mo2','bcc'), ('W2','bcc'), ('Fe2','bcc'),
             ('Nb2','bcc'),   ('Ta2','bcc')]
			
fccs =      [('Rh','fcc'),   ('Ir','fcc'), ('Al','fcc'),('Pb','fcc'),
             ('Ni','fcc'),   ('Pd','fcc'), ('Pt','fcc'),('Au','fcc'),
			 ('Cu','fcc'),   ('Ag','fcc'), ('Ca','fcc'),('Sr','fcc')]

hcps = 		[('Cd2','hcp'),   ('Co2','hcp'), ('Os2','hcp'), ('Ru2','hcp'), 
			 ('Zn2','hcp'),   ('Zr2','hcp'), ('Sc2','hcp'), ('Be2','hcp'), 
			 ('Mg2','hcp'),   ('Ti2','hcp')]

diamonds =  [('Sn2','diamond'),('C2','diamond'),('Si2','diamond'),('Ge2','diamond')]

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
xcDomain      = ['PBE','BEEF','mBEEF']
pwDomain      = [500,700,900,1100,1300,1500,1700]
kptDomain     = [(5,5,5),(8,8,8),(11,11,11),(15,15,15),(4,4,1)]
eConvDomain   = [5] # 1e-5 eV
mixDomain     = [5,10] #0.05, #0.1
nMixDomain    = [5]
maxStepDomain = [500]
nBandDomain   = [12]
sigmaDomain   = [1]   #0.1
xTolDomain    = [5,20] #1e-3 Angstrom
magDomain     = [0,3]
pspDomain     = ['gbrv15pbe','sgs15']


__all__ += ['dftDomain','xcDomain','pwDomain','kptDomain','eConvDomain','mixDomain','nMixDomain'
			,'maxStepDomain','nBandDomain','sigmaDomain','xTolDomain','magDomain'
			,'pspDomain']
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
ni  = TrajDataBulk('NiTraj', '/home/ksb/autoTraj/bulk.traj','trigonal')
nis = TrajDataSurf('NiSurf', '/home/ksb/autoTraj/Nisurf.traj')
bTrajDomain    = [ni]
sTrajDomain    = [nis]

bulkDomain = bTrajDomain + stoichStructDomain

trajDict = {t.name:t for t in bTrajDomain+sTrajDomain}


__all__ += ['trajDict','bulkDomain']
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
	,'Rh':1.42	,'Pd':1.39	,'Ag':1.45	,'Cd':1.44	,'In':1.42	,'Sn':1.69	,'Sb':1.39	,'Te':1.38	,'I':1.39	,'Xe':1.4	,'Cs':2.44	,'Ba':2.15	,'La':2.07	,'Ce':2.04	,'Pr':2.03	,'Nd':2.01	,'Pm':1.99	,'Sm':1.98	,'Eu':1.98	,'Gd':1.96	,'Tb':1.94	,'Dy':1.92	,'Ho':1.92	,'Er':1.89	,'Tm':1.9	,'Yb':1.87	,'Lu':1.87	,'Hf':1.75	,'Ta':1.7	,'W':1.62	,'Re':1.51
	,'Os':1.44	,'Ir':1.41	,'Pt':1.36	,'Au':1.5	,'Hg':1.32	,'Tl':1.45	,'Pb':1.46	,'Bi':1.48	,'Po':1.4	,'At':1.5	,'Rn':1.5	,'Fr':2.6	,'Ra':2.21	,'Ac':2.15	,'Th':2.06	,'Pa':2	,'U':1.96	,'Np':1.9	,'Pu':1.87	,'Am':1.8	,'Cm':1.69}

__all__+=['radDict']