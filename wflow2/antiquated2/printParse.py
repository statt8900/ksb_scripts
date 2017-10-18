#External Modules
import re,sys,math
import numpy as np
#Internal Modules
import misc
"""
Auxillary functions that are occasionally needed to translate between strings and things
"""

def printSleep(n):
	for i in reversed(range(n)):
		print "sleeping .... %d ..."%i ; sleep(1)

def s(x): 		return "'"+str(x)+"'"

def ask(x): return raw_input(x+'\n(y/n)---> ').lower() in ['y','yes']

def parseStoich(x): 
	if not x: return {}
	return {alpha(y) : int(digit(y)) for y in cleanSplit(x,'-')}
def digit(x):    return re.sub("[^0-9]","",x)
def alpha(x):    return re.sub("[^a-zA-Z]", "", x)
def alphaNumSplit(x): return (alpha(x),digit(x))

def printTime(floatHours): 
	intHours = int(floatHours)
	return "%02d:%02d" % (intHours,(floatHours-intHours)*60)

def doubleTime(t):
	"""Doubles time in either HH:MM::SS or HH:MM format. Min time = 1 hr, max time = 40 hr"""
	times = [int(x) for x in t.split(':')]
	HHMMSS = len(times) == 3
	tot = times[0]+times[1]/60.0 + (times[2]/3600.0 if HHMMSS else 0)
	return printTime(min(40,math.ceil(2*tot)))+(':00' if HHMMSS else '')

def abbreviateDict(d): return '\n'.join([str(k)+':'+(str(v) if len(str(v))<100 else '...') for k,v in d.items()])


def parseChemicalFormula(x): return {ele:(1 if n is '' else int(n)) for ele,n in re.findall(r'([A-Z][a-z]*)(\d*)', x)}

def printAds(x): 
	return '__'.join([k+'_'+'-'.join(v) for k,v in x.items()])
def parseAds(x):
	items	= x.split('__')
	keys 	= [i.split('_')[0] for i in items]
	vals	= [i.split('_')[1].split('-') for i in items]
	return {k:v for k,v in zip(keys,vals)}

def parseAds(x): raise NotImplementedError


def printNp(x): return str(x.tolist())
def parseNp(x):
	if x is None: return None
	stripped = ''.join(str(x).split('\n'))
	strIn = stripped.replace('] [',' ; ').replace('[','').replace(']','').replace('(','').replace(')','')
	mOut= np.matrix(strIn,dtype=np.float_)	
	return mOut

def int2List(xnum): return [int(x) for x in str(xnum)]


###################
# STRING GENERATORS
###################
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

def AND(lst): 	
	assert isinstance(lst,list)
	return '('+' AND '.join(lst)+')'
def OR(lst): 	
	assert isinstance(lst,list)
	return '('+' OR '.join(lst)+')'

DEL = 'not deleted'

def INITTRAJ(x): 	return EQ('inittraj',s(x))
def FWID(x): 		return EQ('fwid',x)
def FWIDS(xs): 		return IN('fwid',xs)
def LAUNCHDIR(x): 	return EQ('launchdir',s(x))
def LAUNCHDIRS(xs): return IN('launchdir',map(s,xs))
def STRJOB(x):	 	return EQ('strjob',s(x))

def JOBKIND(x): 	return EQ('jobkind',s(x))
LATTICEOPT 	= JOBKIND('latticeopt')
BULKMOD 	= JOBKIND('bulkmod')
XCCONTRIBS 	= JOBKIND('xc')
RELAX 		= JOBKIND('relax')
VCRELAX 	= JOBKIND('vcrelax')
VIB 		= JOBKIND('vib')
DOS 		= JOBKIND('dos')
NEB 		= JOBKIND('neb')

RELAXORLAT 	= OR([LATTICEOPT,RELAX,VCRELAX])
BULKOPT 	= OR([LATTICEOPT,VCRELAX])

def NAME(x): 		return EQ('name',	s(x))

def DFTCODE(x): return EQ('dftcode',s(x))
GPAW = DFTCODE('gpaw')
QE 	 = DFTCODE('quantumespresso')
def PW(x): return EQ('pw',x)
def XC(x): return EQ('xc',s(x))
PBE 	= XC('PBE')
RPBE 	= XC('RPBE')
BEEF 	= XC('BEEF')
MBEEF 	= XC('mBEEF')
def KPTDEN(x): return EQ('kptden',x)
def PSP(x): return EQ('psp',s(x))
SG15 		= PSP('sg15')
GBRV15PBE 	= PSP('gbrv15pbe')

def DWRAT(x): 	return EQ('dwrat',	x)
def ECONV(x): 	return EQ('econv',	x)
def MIXING(x): 	return EQ('mixing',	x)
def NMIX(x): 	return EQ('nmix',	x)
def MAXSTEP(x): return EQ('maxstep',x)
def NBANDS(x): 	return EQ('nbands',	x)
def SIGMA(x): 	return EQ('sigma',	x)
def FMAX(x): 	return EQ('fmax',	x)
def XTOL(x):	return EQ('xtol',	x)
def STRAIN(x): 	return EQ('strain',	x)

def RELAXED(x): 	return LIKE('trajinfo',"'u'relaxed': True''")
def UNRELAXED(x): 	return LIKE('trajinfo',"'u'relaxed': False''")

def STATUS(x): return EQ('status',s(x))
READY 		= STATUS('ready')
QUEUED 		= STATUS('queued')
FIZZLED 	= STATUS('fizzled')
TIMEOUT 	= STATUS('timeout')
CANCELLED	= STATUS('cancelled')
FAILED  	= STATUS('failed')
COMPLETED 	= STATUS('completed')
RUNNING     = STATUS('running')
ARCHIVED  	= STATUS('archived')
NOTCOMPLETED = NEQ('status',s('complete'))

KOHNSHAM = LIKE('trace', '%KohnShamConvergenceError%')
NOTKOHN = OR([NULL('trace'), NOTLIKE('trace','%KohnShamConvergenceError%')])

FWARCHIVED = EQ('fwstatus',s('ARCHIVED'))

H2 	= EQ('numbers_str',"'[1, 1]'")
H2O = EQ('name',"'H2O'")
O2 	= EQ('numbers_str',"'[8, 8]'")
N2 	= EQ('numbers_str',"'[7, 7]'")
LI 	= EQ('numbers_str',"'[3]'")
CU 	= EQ('numbers_str',"'[29]'")
BE 	= EQ('numbers_str',"'[4, 4]'")
AL 	= EQ('numbers_str',"'[13]'")
PD 	= EQ('numbers_str',"'[46]'")
AU 	= EQ('numbers_str',"'[79]'")
PT 	= EQ('numbers_str',"'[78]'")
NI 	= EQ('numbers_str',"'[28]'")
ZN 	= EQ('numbers_str',"'[30, 30]'")

##########################
# Kind-related constraints
##########################
def KIND(x): return EQ('kind',s(x))
SURFACE 	= KIND('surface')
BULK 		= KIND('bulk')
MOLECULE 	= KIND('molecule')


def STRUCTURE(x): return EQ('structure',s(x))
HCP 	= STRUCTURE('hexagonal')
FCC 	= STRUCTURE('fcc')
BCC 	= STRUCTURE('bcc')
DIAMOND = STRUCTURE('diamond')


def SURFXY(x):		return EQ('xy',			x)
def LAYERS(x): 		return EQ('layers',		x)
def CONSTRAINED(x): return EQ('constrained',x)
def SYMMETRIC(x): 	return EQ('symmetric',	int(x))
def VACUUM(x): 		return EQ('vacuum',		x)
def ADSORBATES(x): 	return EQ('stradsorbates',s(x)) #  https://stackoverflow.com/questions/603572/how-to-properly-escape-a-single-quote-for-a-sqlite-database

#######

def EXISTS(filename): return EQ(filename.replace('.','_'),1)

BFIT_OK 		= GT('bfit', 0.5)
PW_OVER_1000 	= GT('pw', 	1000)
PW_UNDER_1000 	= LT('pw', 	1000)


NONMETAL 		= NEQ('metalspecies_str','species_str')#AND([,NOTLIKE('name',s('%C-diamond%'))])
INTERSTITIAL 	= AND([LATTICEOPT,NONMETAL])
ZEROFORCE 		= EQ('emtsym',1)
DOF1 			= EQ('dof',1)

QOVER40 		= OR(["walltime like '%s'"%x for x in [str(i)+'_:__:__' for i in range(4,10)]+[str(i)+'_:__' for i in range(4,10)]])
### REFERENCES
#Constraints for references should specify a linearly-independent combination of elements, 
#given a set of parameters:  (xc,pw,kptden,psp,xtol,strain,convid,precalc,dftcode)

REF_CONSTRAINT = AND([RELAXORLAT,IN('name',[s(x) for x in misc.refNames])])
NOREF = 'refeng is null'
###################################################################
### For equivalence relations - used in conjunction with queryTuple
###################################################################
def EQUAL(col): return 'j0.%s=j1.%s'%(col,col)

EQUALCALC = AND([EQUAL('pw'),EQUAL('xc'),EQUAL('kptden'),EQUAL('dftcode')
								,EQUAL('psp'),EQUAL('dwrat'),EQUAL('econv'),EQUAL('fmax')])

LATOPTDOF 	= AND(["j0.jobkind='latticeopt'",EQUAL('jobkind'),EQUAL('dof')])


#class Query(object):
#	def __init__(self,select,frm='job',where='1',groupby=None,having=None,orderby=None):
#		self.select=select;self.frm=frm;self.where=where;self.groupby=groupby;self.having=having;self.orderby=orderby

