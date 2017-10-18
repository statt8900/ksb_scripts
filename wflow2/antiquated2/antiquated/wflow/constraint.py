"""
- Defines a 'constraint', which is used to filter database entries of jobs (data.db) 
- Rather than write out a string each time you want to select, db.plotQuery() will return columns that satisfy any list of Constraints
- Basic Constraints are conditions imposed on some column value. They can be converted to strings for Sql querying. 
- They can be composed with OR/AND to make complex constraints.
- Ideally, this is the 'language' for specifying which sets of jobs you want to initialize, submit, analyze, or manipulate in any way
"""

class Constraint(object):
	def __init__(self,column,value,relation):
		self.col = column
		self.val = value #strings must have ' in the (python) string
		self.rel = relation
	def __str__(self):
		if   self.rel == '=':    		return self.col + '='+ str(self.val)
		elif self.rel == '!=':   		return self.col + '!='+ str(self.val)
		elif self.rel == '=list': 		return ' ('+self.col + ' = ' + (' OR '+self.col+' = ').join(self.val)  + ' ) '
		elif self.rel == '>':    		return self.col + '>'+ str(self.val)
		elif self.rel == '<':    		return self.col + '<'+ str(self.val)
		elif self.rel == 'btw':  		return self.col + str(self.val[0]) + ' AND ' + str(self.val[1])
		
		elif self.rel == 'like': 		return self.col + ' LIKE ' + self.val 												# https://www.techonthenet.com/sqlite/like.php
		elif self.rel == 'notlike': 	return self.col + 'NOT LIKE ' + self.val 											# https://www.techonthenet.com/sqlite/like.php
		elif self.rel == 'notlikelist': return self.col + ' NOT LIKE ' + (' AND '+self.col+' NOT LIKE ').join(self.val) 	# https://www.techonthenet.com/sqlite/like.php
		elif self.rel == 'likelist': 	return ' (' + self.col + ' LIKE ' + (' OR ' + self.col + ' LIKE ').join(self.val) + ')'	# https://www.techonthenet.com/sqlite/like.php
		elif self.rel == 'or':          return '((' + str(self.val[0]) + ' ) OR (' + str(self.val[1]) + '))'
		elif self.rel == 'and':         return '((' + str(self.val[0]) + ' ) AND (' + str(self.val[1]) + '))'
		elif self.rel == 'andlist':     return '((' + ' ) AND ('.join([str(x) for x in self.val])+ '))'
		elif self.rel == 'orlist':     	return '((' +  ' ) OR ('.join([str(x) for x in self.val])+ '))'
		else: raise NotImplementedError
		# change AND/OR to take in a list of arbitrary length?


#############################
### KELD DATABASE COMPOSITION
#############################
AG   = Constraint('name',"'%Ag-fcc%'",		'like')
AL   = Constraint('name',"'%Al-fcc%'",		'like')
AU   = Constraint('name',"'%Au-fcc%'",		'like')
BA   = Constraint('name',"'%Ba-bcc%'",		'like')
BE   = Constraint('name',"'%Be-hcp%'",		'like')
C    = Constraint('name',"'%C-diamond%'",	'like')
CA   = Constraint('name',"'%Ca-fcc%'",		'like')
CD   = Constraint('name',"'%Cd-hcp%'",		'like')
CO   = Constraint('name',"'%Co-hcp%'",		'like')
CU   = Constraint('name',"'%Cu-fcc%'",		'like')
FE   = Constraint('name',"'%Fe-bcc%'",		'like')
GE   = Constraint('name',"'%Ge-diamond%'",	'like')
IR   = Constraint('name',"'%Ir-fcc%'",		'like')
K    = Constraint('name',"'%K-bcc%'",		'like')
LI   = Constraint('name',"'%Li-bcc%'",		'like')
MG   = Constraint('name',"'%Mg-hcp%'",		'like')
MO   = Constraint('name',"'%Mo-bcc%'",		'like')
NA   = Constraint('name',"'%Na-bcc%'",		'like')
NB   = Constraint('name',"'%Nb-bcc%'",		'like')
NI   = Constraint('name',"'%Ni-fcc%'",		'like')
OS   = Constraint('name',"'%Os-hcp%'",		'like')
PB   = Constraint('name',"'%Pb-fcc%'",		'like')
PD   = Constraint('name',"'%Pd-fcc%'",		'like')
RB   = Constraint('name',"'%Rb-bcc%'",		'like')
RH   = Constraint('name',"'%Rh-fcc%'",		'like')
RU   = Constraint('name',"'%Ru-hcp%'",		'like')
SC   = Constraint('name',"'%Sc-hcp%'",		'like')
SI   = Constraint('name',"'%Si-diamond%'",	'like')
SN   = Constraint('name',"'%Sn-diamond%'",	'like')
SR   = Constraint('name',"'%Sr-fcc%'",		'like')
TI   = Constraint('name',"'%Ti-hcp%'",		'like')
ZN   = Constraint('name',"'%Zn-hcp%'",		'like')
ZR   = Constraint('name',"'%Zr-hcp%'",		'like')

ELEMENTS = [AG,AL,AU,BA,BE,C,CA,CD,CO,CU,FE,GE,IR,K,LI,MG,MO,NA,NB,NI
			,OS,PB,PD,RB,RH,RU,SC,SI,SN,SR,TI,ZN,ZR]

SINGLE_ELEMENTS 	= Constraint('name',[x.val for x in ELEMENTS], 			'likelist')
VARIOUS_ELEMENTS	= Constraint('name',[x.val for x in [LI,BE,AL,CU,C]],	'likelist')


####################
# SIMPLE CONSTRAINTS
####################


def STRUCTURE(x): return Constraint('structure',"'"+x+"'", '=')
if True:
	HCP 	= STRUCTURE('hexagonal')
	FCC 	= STRUCTURE('fcc')
	BCC 	= STRUCTURE('bcc')
	DIAMOND = STRUCTURE('diamond')





SYMMETRIC 	= Constraint('symmetric',"'True'",	'=')
ASYMMETRIC 	= Constraint('symmetric',"'False'",	'=')


def JOBID(x): 	return Constraint('jobid',	x,			'=')

def JOBKIND(x): return Constraint('jobkind',"'"+x+"'",	'=')	
if True:
	BULKRELAX 	= JOBKIND('bulkrelax')
	SURFRELAX 	= JOBKIND('surfrelax')
	VIB 		= JOBKIND('vib')
	DOS 		= JOBKIND('dos')
	NEB 		= JOBKIND('neb')

def PW(x): return Constraint('pw',x,'=')

def XC(x): return Constraint('xc',"'"+x+"'",'=')
if True:
	PBE 	= XC('PBE')
	RPBE 	= XC('RPBE')
	BEEF 	= XC('BEEF')
	MBEEF 	= XC('MBEEF')


def KPTDEN(x): return Constraint('kptden',x,'=')
if True:
	KPTDEN2 = KPTDEN(2)

def PSP(x): return Constraint('psp',"'"+x+"'",'=')
if True:
	SG15 		= PSP('sg15')
	GBRV15PBE 	= PSP('gbrv15pbe')


def XTOL(x):	return Constraint('xtol',	x,'=')
def STRAIN(x): 	return Constraint('strain',	x,'=')
def CONVID(x): 	return Constraint('convid',	x,'=')

def PRECALC(x): return Constraint('precalc',"'"+x+"'",'=')
if True:
	PBE_PRECALC = PRECALC('PBE')
	NO_PRECALC 	= PRECALC('None')

def DFTCODE(x): return Constraint('dftcode',"'"+x+"'",	'=')
if True:
	GPAW = DFTCODE('gpaw')
	QE 	 = DFTCODE('quantumespresso')

def STATUS(x): 	return Constraint('status',	"'"+x+"'",	'=')
if True:
	INITIALIZED = STATUS('initialized')
	COMPLETED 	= STATUS('complete')
	FIZZLED  	= STATUS('fizzled')
	QUEUED      = STATUS('queued')

def KIND(x): 	return Constraint('kind',	"'"+x+"'", 	'=')
if True:
	SURFACE 	= KIND('surface')
	BULK 		= KIND('bulk')
	MOLECULE 	= KIND('molecule')

def LAYERS(x): 	return Constraint('layers',	x,'=')
def SURFXY(x):	return Constraint('xy',		x,'=')
def VACUUM(x):	return Constraint('vacuum',	x,'=')

def FACET(x): return Constraint('facet',x,'=')
if True:
	F111 	= FACET(111)
	F100	= FACET(100)




######
# LIKE
#######

KOHNSHAM 	= Constraint('error',"'%KohnShamConvergenceError%'",'like')
TIMEOUT 	= Constraint('error',"'%TIMEOUT%'",'like')

##########
# < > > <
##########

BFIT_OK 		= Constraint('bfit', 0.5,	'>')
PW_OVER_1000 	= Constraint('pw', 	1000,	'>')
PW_UNDER_1000 	= Constraint('pw', 	1000,	'<')

###########
# Negations
###########
INCOMPLETED = Constraint('status',"'complete'",'!=')

NOT_MBEEF 	= Constraint('xc', 		"'mBEEF'", 			'!=')
NOT_BEEF    = Constraint('xc', 		"'BEEF'", 			'!=')	
NOT_GPAW 	= Constraint('dftcode',"'gpaw'",			'!=')
NOT_QE 		= Constraint('dftcode',"'quantumespresso'",	'!=')

NO_VACANCY 	 = Constraint('vacancies',"null",	'=') #how should this work?
NO_ADSORBATE = Constraint('adsorbates',"null",	'=')

###############
# COMPOSITE
###############

GPAW_SG15 	= Constraint('',(GPAW,SG15), 	'and')
QE_GBRV15PBE = Constraint('',(QE,GBRV15PBE),	'and')

PBE_MBEEF_CALCS 	= Constraint('',(MBEEF,PBE_PRECALC),		'and')
NO_PRECALC_NOT_MBEEF 	= Constraint('',(NOT_MBEEF,NO_PRECALC),	'and')

PSPGUARD 			= Constraint('',(GPAW_SG15,QE_GBRV15PBE),					'or')
PRECALC_WHEN_NEEDED = Constraint('',(PBE_MBEEF_CALCS,NO_PRECALC_NOT_MBEEF),	'or')
NO_BEEF_IN_GPAW 	= Constraint('',(NOT_GPAW,NOT_BEEF),	'or')
NO_MBEEF_IN_QE 		= Constraint('',(NOT_QE,NOT_MBEEF),	'or')

######################
# REFERENCES
######################
"""
Constraints for references should specify a linearly-independent combination of elements, 
given a set of parameters: 
	(xc,pw,kptden,psp,xtol,strain,convid,precalc,dftcode)
"""

refNames = [u'Al-fcc', u'Li-bcc', u'Be-hcp', u'C-diamond', u'Rh-fcc', u'Ir-fcc', u'Pb-fcc', u'Ni-fcc', u'Pd-fcc', u'Au-fcc', u'Cu-fcc', u'Ag-fcc', u'Ca-fcc', u'Sr-fcc', u'Na-bcc', u'K-bcc', u'Rb-bcc', u'Ba-bcc', u'Mo-bcc', u'Fe-bcc', u'Nb-bcc', u'Cd-hcp', u'Co-hcp', u'Os-hcp', u'Ru-hcp', u'Zn-hcp', u'Zr-hcp', u'Sc-hcp', u'Mg-hcp', u'Ti-hcp', u'Sn-diamond', u'Si-diamond', u'Ge-diamond']
REF_CONSTRAINTS = [Constraint('name',["'"+x+"'" for x in refNames],'=list')]

######################
# LISTS OF CONSTRAINTS
######################

SINGLE_ELEMENT_FILTERS 	= [[x] for x in ELEMENTS]
SINGLE_ELEMENT_LABELS	= [x.val.split('%')[1] for x in ELEMENTS]
TEST_FILTERS 			= [INITIALIZED,SINGLE_ELEMENTS,PW_UNDER_1000,KPTDEN(2)]

ESSENTIAL_BULKJOB_CONSTRAINTS = [PSPGUARD
								,PRECALC_WHEN_NEEDED
								,NO_BEEF_IN_GPAW
								,NO_MBEEF_IN_QE]   # change with caution

ESSENTIAL_SURFJOB_CONSTRAINTS = []

