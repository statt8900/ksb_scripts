"""
Designed as an enhancement to constraint.py, we now allow constraints that 
can be used to select multiple rows of a column, with constraints that relate the rows
"""
import itertools

from printParse import s
from misc 		import partition,refNames,makeEquivRelation
from dbase 		import queryTuple,queryCol

class Constraint(object):
	def __init__(self,n,col,val,rel):
		self.n=n;self.col=col;self.val=val;self.rel=rel
	def table(self): return 'detail%d'%(self.n - 1)
	def sql(self):

		def simpleBinary(symb): return '%s %s %s'%(self.col,self.rel,str(self.val))

		if self.rel in ['=','!=','like','not like','>','<']: 
			return simpleBinary(self.rel)
		elif self.rel == 'in': 
			strval = '('+','.join(self.val) + ')'
			if   len(self.val)==1: strval = strval.replace(',','')
			elif len(self.val)==0: return ' 0 '
			return '%s %s %s'%(self.col,self.rel,strval)
		elif self.rel == 'and':         
			assert isinstance(self.val,list), 							'"and" Constraint was not given a list in its "val" field'
			assert all([isinstance(x,Constraint) for x in self.val]), 	'"and" Constraint has something weird in the list'
			assert all([x.n == self.val[0].n for x in self.val]), 		'"and" Constraint combining constraints with different arities'
			assert self.col is None, 									'"and" Constraint has a col field defined'
			return '('+' AND '.join([x.sql() for x in self.val])+')'

		elif self.rel is None: return '1'

def NONE(n): 	return Constraint(n,None,None,None)
def AND(cList): return Constraint(cList[0].n,None,cList,'and')


##############################################################
def INITTRAJ(x): return Constraint(1,'inittraj',"'"+x+"'",'=')
def FWID(x): return Constraint(1,'fwid',x,'=')
def FWIDS(x): return Constraint(1,'fwid',x,'in')
def LAUNCHDIR(x): return Constraint(1,'launchdir',s(x) ,'=')
def LAUNCHDIRS(x): return Constraint(1,'launchdir',[s(y) for y in x],'in')
def JOBKIND(x): return Constraint(1,'jobkind',s(x),'=')	
LATTICEOPT 	= JOBKIND('latticeopt')
BULKMOD 	= JOBKIND('bulkmod')
XCCONTRIBS 	= JOBKIND('xc')
RELAX 		= JOBKIND('relax')
VIB 		= JOBKIND('vib')
DOS 		= JOBKIND('dos')
NEB 		= JOBKIND('neb')

RELAXORLAT 	= OR([LATTICEOPT,RELAX])
CALC 		= OR([LATTICEOPT,RELAX,BULKMOD,XCCONTRIBS,VIB,DOS,NEB]) #all jobs currently use a calc

def SYMBOLS(x): 	return Constraint(1,'strsymbols',	s(x),'=')
def NUMBERS(x): 	return Constraint(1,'strnumbers',	s(x),'=')
def INITPOS(x): 	return Constraint(1,'strinitpos',	s(x),'=')
def INITCELL(x): 	return Constraint(1,'strinitcell',	s(x),'=')
def MAGMOMINIT(x): 	return Constraint(1,'strmagmominit',s(x),'=')

def DFTCODE(x): return Constraint(1,'dftcode',s(x),'=')
GPAW = DFTCODE('gpaw')
QE 	 = DFTCODE('quantumespresso')

def PW(x): return Constraint(1,'pw',x,'=')

def XC(x): return Constraint(1,'xc',s(x),'=')
PBE 	= XC('PBE')
RPBE 	= XC('RPBE')
BEEF 	= XC('BEEF')
MBEEF 	= XC('MBEEF')

def KPTDEN(x): return Constraint(1,'kptden',x,'=')

def PSP(x): return Constraint(1,'psp',s(x),'=')
SG15 		= PSP('sg15')
GBRV15PBE 	= PSP('gbrv15pbe')

def DWRAT(x): 	return Constraint(1,'dwrat',	x,'=')
def ECONV(x): 	return Constraint(1,'econv',	x,'=')
def MIXING(x): 	return Constraint(1,'mixing',	x,'=')
def NMIX(x): 	return Constraint(1,'nmix',		x,'=')
def MAXSTEP(x): return Constraint(1,'maxstep',	x,'=')
def NBANDS(x): 	return Constraint(1,'nbands',	x,'=')
def SIGMA(x): 	return Constraint(1,'sigma',	x,'=')
def FMAX(x): 	return Constraint(1,'fmax',		x,'=')
def XTOL(x):	return Constraint(1,'xtol',		x,'=')
def STRAIN(x): 	return Constraint(1,'strain',	x,'=')

def RELAXED(x): 	return Constraint(1,'trajinfo',"'u'relaxed': True''",'like')
def UNRELAXED(x): 	return Constraint(1,'trajinfo',"'u'relaxed': False''",'like')

def ROW(arity,rowids): 	return Constraint(arity,'id%d'%(arity-1),rowids,'in') 

def STATUS(x): return Constraint(1,'status',s(x),'=')
READY 		= STATUS('ready')
QUEUED 		= STATUS('queued')
FIZZLED 	= STATUS('fizzled')
TIMEOUT 	= STATUS('timeout')
CANCELLED	= STATUS('cancelled')
FAILED  	= STATUS('failed')
COMPLETED 	= STATUS('completed')
RUNNING     = STATUS('running')

NOTCOMPLETED = Constraint(1,'status',s('completed'),'!=')


# 'unlocked' means that we presently INTEND to eventually update the row
# a job that has not finished should be unlocked
# if a job is finished, then it should be unlocked 

##########
# ELEMENTS
##########
AG   = Constraint(1,'name',s('%Ag-fcc%'),	'like')
AL   = Constraint(1,'name',s('%Al-fcc%'),	'like')
AU   = Constraint(1,'name',s('%Au-fcc%'),	'like')
BA   = Constraint(1,'name',s('%Ba-bcc%'),	'like')
BE   = Constraint(1,'name',s('%Be-hcp%'),	'like')
C    = Constraint(1,'name',s('%C-diamond%'),'like')
CA   = Constraint(1,'name',s('%Ca-fcc%'),	'like')
CD   = Constraint(1,'name',s('%Cd-hcp%'),	'like')
CO   = Constraint(1,'name',s('%Co-hcp%'),	'like')
CU   = Constraint(1,'name',s('%Cu-fcc%'),	'like')
FE   = Constraint(1,'name',s('%Fe-bcc%'),	'like')
GE   = Constraint(1,'name',s('%Ge-diamond%'),'like')
IR   = Constraint(1,'name',s('%Ir-fcc%'),	'like')
K    = Constraint(1,'name',s('%K-bcc%'),	'like')
LI   = Constraint(1,'name',s('%Li-bcc%'),	'like')
MG   = Constraint(1,'name',s('%Mg-hcp%'),	'like')
MO   = Constraint(1,'name',s('%Mo-bcc%'),	'like')
NA   = Constraint(1,'name',s('%Na-bcc%'),	'like')
NB   = Constraint(1,'name',s('%Nb-bcc%'),	'like')
NI   = Constraint(1,'name',s('%Ni-fcc%'),	'like')
OS   = Constraint(1,'name',s('%Os-hcp%'),	'like')
PB   = Constraint(1,'name',s('%Pb-fcc%'),	'like')
PD   = Constraint(1,'name',s('%Pd-fcc%'),	'like')
RB   = Constraint(1,'name',s('%Rb-bcc%'),	'like')
RH   = Constraint(1,'name',s('%Rh-fcc%'),	'like')
RU   = Constraint(1,'name',s('%Ru-hcp%'),	'like')
SC   = Constraint(1,'name',s('%Sc-hcp%'),	'like')
SI   = Constraint(1,'name',s('%Si-diamond%'),'like')
SN   = Constraint(1,'name',s('%Sn-diamond%'),'like')
SR   = Constraint(1,'name',s('%Sr-fcc%'),	'like')
TI   = Constraint(1,'name',s('%Ti-hcp%'),	'like')
ZN   = Constraint(1,'name',s('%Zn-hcp%'),	'like')
ZR   = Constraint(1,'name',s('%Zr-hcp%'),	'like')

ELEMENTS = [AG,AL,AU,BA,BE,C,CA,CD,CO,CU,FE,GE,IR,K,LI,MG,MO,NA,NB,NI
			,OS,PB,PD,RB,RH,RU,SC,SI,SN,SR,TI,ZN,ZR]

H2 = Constraint(1,'numbers','[1,1]','=')

##########################
# Kind-related constraints
##########################
def KIND(x): return Constraint(1,'kind',s(x), '=')
SURFACE 	= KIND('surface')
BULK 		= KIND('bulk')
MOLECULE 	= KIND('molecule')


def STRUCTURE(x): return Constraint(1,'structure',s(x), '=')
HCP 	= STRUCTURE('hexagonal')
FCC 	= STRUCTURE('fcc')
BCC 	= STRUCTURE('bcc')
DIAMOND = STRUCTURE('diamond')


def SURFXY(x):		return Constraint(1,'xy',			x,			'=')
def LAYERS(x): 		return Constraint(1,'layers',		x,			'=')
def CONSTRAINED(x): return Constraint(1,'constrained',	x,			'=')
def SYMMETRIC(x): 	return Constraint(1,'symmetric',	s(str(x)),	'=')
def VACUUM(x): 		return Constraint(1,'vacuum',		x,			'=')
def ADSORBATES(x): 	return Constraint(1,'stradsorbates',s(x),		'=') #  https://stackoverflow.com/questions/603572/how-to-properly-escape-a-single-quote-for-a-sqlite-database

#######

def EXISTS(filename): return Constraint(1,filename.replace('.','_'),s(True),'=')

BFIT_OK 		= Constraint(1,'bfit', 0.5,		'>')
PW_OVER_1000 	= Constraint(1,'pw', 	1000,	'>')
PW_UNDER_1000 	= Constraint(1,'pw', 	1000,	'<')

INTERSTITIAL 	= AND([LATTICEOPT,Constraint(1,'stoich','metalstoich','!='),Constraint(1,'name',s('%C-diamond%'),'not like')])
ZEROFORCE 		= Constraint(1,'structure',['hexagonal','hcp','fcc','bcc','diamond'], 'in')


### REFERENCES
#Constraints for references should specify a linearly-independent combination of elements, 
#given a set of parameters:  (xc,pw,kptden,psp,xtol,strain,convid,precalc,dftcode)

REF_CONSTRAINT = Constraint(1,'name',[s(x) for x in refNames],'in')

############################################################################################
############################################################################################
############################################################################################


class Equivalence(object):
	def __init__(self,detail2,domainConstraint=NONE(1)):
		assert domainConstraint.n == 1
		self.detail2 			= detail2 			# String given to sqlite that picks out pairs of jobs. All jobs that are linked by the relation will be grouped together
		self.domainConstraint 	= domainConstraint 		# CONSTRAINT of arity 1
	
	def output(self,extraConstraint=NONE(1)): 
		print 'entering OUTPUT'
		dom = queryCol('launchdir',self.domainConstraint)# AND([self.domainConstraint]))#,extraConstraint]))
		#print 'in equiv output: len dom = ',len(dom)
		#print "queryTuple(2,self.querystrs) ",queryTuple(2,self.querystrs)
		return [(x,y) for (x,y) in sqlexecute('SELECT ld0,ld1 from detail1 where '+self.detail2) if (x in dom) and (y in dom)]

	def classes(self,extraConstraint=NONE(1)): 
		print 'entering classes'
		
		output 				= self.output(extraConstraint)
		domain 				= allElems(output)
		outputEquivRelation = makeEquivRelation(output,domain)
		#print 'in classes: output = ',output
		#print 'in classes: domain = ',domain
		#print 'in classes: outputEquivRelation = ',outputEquivRelation
		return partition(domain,lambda x,y: (x,y) in outputEquivRelation)
	def getClass(self,launchdir): 
		print 'entering getClass'
		#print 'in GetClass: launchdir = ',launchdir
		#print 'computing self.classes()\n\t ',self.classes()
		return [x for x in self.classes() if launchdir in x][0]			

def allElems(listoflists):
	return list(set(itertools.chain(*listoflists))) 

# general equality of some field for 2 things
def EQUAL(col): return 'j0.%s=j1.%s'%(col,col)

EQUALCALC = Equivalence([EQUAL('pw'),EQUAL('xc'),EQUAL('kptden'),EQUAL('dftcode')
								,EQUAL('psp'),EQUAL('dwrat'),EQUAL('econv'),EQUAL('fmax')])

EQUALNAME = Equivalence([EQUAL('name')])

#TRIVIALEQ = Equivalence([EQUAL('id0')])
TRIVIALEQ 	= Equivalence('1',NONE(1))
JOBKIND 	= Equivalence([EQUAL('jobkind')],NONE(1))
LATOPTDOF 	= Equivalence([EQUAL('jobkind'),EQUAL('dof')],NONE(1))

#EQUALCALC.classes()
# equal calcs
