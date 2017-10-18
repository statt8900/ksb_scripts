#Internal Modules
import misc


def s(x): 		return "'"+str(x)+"'"

def SIMPLEBINARY(a,b,c): return '%s %s %s'%(a,b,c)
def EQ(a,b): 		return SIMPLEBINARY(a,'=',b)
def NEQ(a,b):		return SIMPLEBINARY(a,'!=',b)
def GT(a,b):		return SIMPLEBINARY(a,'>',b)
def LT(a,b):		return SIMPLEBINARY(a,'<',b)
def LIKE(a,b):		return SIMPLEBINARY(a,'like',s(b))
def NOTLIKE(a,b):	return SIMPLEBINARY(a,'not like',s(b))
def NULL(a):		return SIMPLEBINARY(a,'is','null')

def IN(a,lst):
	assert isinstance(lst,list)
	strlst = '('+','.join([str(x) for x in lst]) + ')'
	if   len(lst)==1: strlst = strlst.replace(',','')
	elif len(lst)==0: return ' 0 '
	return SIMPLEBINARY(a,'in',strlst)

def OR(*lst): 	return '('+' OR '.join(lst)+')'
def AND(*lst): return '('+' AND '.join(lst)+')'

DEL = 'not deleted'

def FWID(x): 		return EQ('fwid',x)
def FWIDS(xs): 		return IN('fwid',xs)
def STRJOB(x):	 	return EQ('strjob',s(x))

def JOBKIND(x): 	return EQ('job_type_ksb',s(x))
LATTICEOPT 	= JOBKIND('latticeopt')
BULKMOD 	= JOBKIND('bulkmod')
XCCONTRIBS 	= JOBKIND('xc')
RELAX 		= JOBKIND('relax')
VCRELAX 	= JOBKIND('vcrelax')
VIB 		= JOBKIND('vib')
DOS 		= JOBKIND('dos')
NEB 		= JOBKIND('neb')

RELAXORLAT 	= OR(LATTICEOPT,RELAX,VCRELAX)
BULKOPT 	= OR(LATTICEOPT,VCRELAX)

def NAME(x): 		return EQ('job_name',	s(x))

def DFTCODE(x): return EQ('dftcode',s(x))
GPAW = DFTCODE('gpaw')
QE 	 = DFTCODE('quantumespresso')
def PW(x): return EQ('planewave_cutoff',x)
def XC(x): return EQ('xc',s(x))
PBE 	= XC('PBE')
RPBE 	= XC('RPBE')
BEEF 	= XC('BEEF')
MBEEF 	= XC('mBEEF')
def KPTDEN(x): return EQ('kptden_ksb',x)
def PSP(x): return EQ('psp_ksb',s(x))
SG15 		= PSP('sg15')
GBRV15PBE 	= PSP('gbrv15pbe')
PAW 		= PSP('paw')

def DWRAT(x): 	return EQ('dwrat_ksb',	x)
def ECONV(x): 	return EQ('econv_ksb',	x)
def MIXING(x): 	return EQ('mixing_ksb',	x)
def NMIX(x): 	return EQ('nmix_ksb',	x)
def MAXSTEP(x): return EQ('maxstep_ksb',x)
def NBANDS(x): 	return EQ('nbands_ksb',	x)
def SIGMA(x): 	return EQ('sigma_ksb',	x)
def FMAX(x): 	return EQ('fmax_ksb',	x)
def XTOL(x):	return EQ('xtol_ksb',	x)
def STRAIN(x): 	return EQ('strain_ksb',	x)

H2 	= EQ('numbers_str',"'[1, 1]'")
O2 	= EQ('numbers_str',"'[8, 8]'")
N2 	= EQ('numbers_str',"'[7, 7]'")
LI 	= EQ('numbers_str',"'[3]'")
BE 	= EQ('numbers_str',"'[4, 4]'")
NA 	= EQ('numbers_str',"'[11]'")
MG 	= EQ('numbers_str',"'[12, 12]'")
AL 	= EQ('numbers_str',"'[13]'")
K   = EQ('numbers_str',"'[19]'")
CA  = EQ('numbers_str',"'[20]'")
NI 	= EQ('numbers_str',"'[28]'")
CU 	= EQ('numbers_str',"'[29]'")
ZN 	= EQ('numbers_str',"'[30, 30]'")
RB  = EQ('numbers_str',"'[37]'")
SR  = EQ('numbers_str',"'[38]'")
RU  = EQ('numbers_str',"'[44, 44]'")
PD 	= EQ('numbers_str',"'[46]'")
CS  = EQ('numbers_str',"'[55]'")
BA  = EQ('numbers_str',"'[56]'")
IR  = EQ('numbers_str',"'[77]'")
PT 	= EQ('numbers_str',"'[78]'")
AU 	= EQ('numbers_str',"'[79]'")
H2O = EQ('name',"'H2O'")

PAWELEM = OR(LI,BE,NA,MG,K,CA,RB,SR,CS,BA,ZN)

PAWPSP = OR(AND(GPAW,PAW,'pawjob'),AND(GPAW,SG15,'not pawjob'),QE)
# need be, mg, 

##########################
# Kind-related constraints
##########################
def KIND(x): return EQ('system_type_ksb',s(x))
SURFACE 	= KIND('surface')
BULK 		= KIND('bulk')
MOLECULE 	= KIND('molecule')


def STRUCTURE(x): return EQ('structure_ksb',s(x))
HCP 	= STRUCTURE('hexagonal')
FCC 	= STRUCTURE('fcc')
BCC 	= STRUCTURE('bcc')
DIAMOND = STRUCTURE('diamond')


def SURFXY(x):		return EQ('xy_ksb',			x)
def LAYERS(x): 		return EQ('layers_ksb',		x)
def CONSTRAINED(x): return EQ('constrained_ksb',x)
def SYMMETRIC(x): 	return EQ('symmetric_ksb',	int(x))
def VACUUM(x): 		return EQ('vacuum_ksb',		x)
def ADSORBATES(x): 	return EQ('adsorbates_ksb',s(x)) #  https://stackoverflow.com/questions/603572/how-to-properly-escape-a-single-quote-for-a-sqlite-database

#######

def EXISTS(filename): return EQ(filename.replace('.','_'),1)

BFIT_OK 		= GT('bulk_modulus_quadfit', 0.5)
PW_OVER_1000 	= GT('planewave_cutoff', 	1000)
PW_UNDER_1000 	= LT('planewave_cutoff', 	1000)

NONMETAL 		= NEQ('metalspecies_str','species_str')
INTERSTITIAL 	= AND(LATTICEOPT,NONMETAL)
ZEROFORCE 		= EQ('emtsym',1)
DOF1 			= EQ('dof',1)

### REFERENCES
#Constraints for references should specify a linearly-independent combination of elements, 
#given a set of parameters:  (xc,planewave_cutoff,kptden_ksb,psp,xtol,strain,convid,precalc,dftcode)
REF_CONSTRAINT = AND(RELAXORLAT,IN('name',[s(x) for x in misc.refNames]))
NOREF = 'refeng is null'

"""
###################################################################
### For equivalence relations - used in conjunction with queryTuple
###################################################################
def EQUAL(col): return 'j0.%s=j1.%s'%(col,col)
EQUALCALC = AND(EQUAL('planewave_cutoff'),EQUAL('xc'),EQUAL('kptden_ksb'),EQUAL('dftcode') ,EQUAL('psp'),EQUAL('dwrat'),EQUAL('econv'),EQUAL('fmax'))
LATOPTDOF 	= AND("j0.jobkind='latticeopt'",EQUAL('jobkind'),EQUAL('dof'))
"""
