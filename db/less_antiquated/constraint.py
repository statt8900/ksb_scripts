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
		elif self.rel == '=list': 		return ' ('+self.col + ' = ' + (' OR '+self.col+' = ').join(self.val) + ')'
		elif self.rel == '>':    		return self.col + '>'+ str(self.val)
		elif self.rel == '<':    		return self.col + '<'+ str(self.val)
		elif self.rel == 'btw':  		return self.col + str(self.val[0]) + ' AND ' + str(self.val[1])
		
		elif self.rel == 'like': 		return self.col + ' LIKE ' + self.val 												# https://www.techonthenet.com/sqlite/like.php
		elif self.rel == 'notlike': 	return self.col + 'NOT LIKE ' + self.val 											# https://www.techonthenet.com/sqlite/like.php
		elif self.rel == 'notlikelist': return self.col + ' NOT LIKE ' + (' AND '+self.col+' NOT LIKE ').join(self.val) 	# https://www.techonthenet.com/sqlite/like.php
		elif self.rel == 'likelist': 	return ' (' + self.col + ' LIKE ' + (' OR ' + self.col + ' LIKE ').join(self.val) + ')'	# https://www.techonthenet.com/sqlite/like.php
		elif self.rel == 'or':          return '((' + str(self.val[0]) + ' ) OR (' + str(self.val[1]) + '))'
		elif self.rel == 'and':         return '((' + str(self.val[0]) + ' ) AND (' + str(self.val[1]) + '))'
		else: raise NotImplementedError
		# change AND/OR to take in a list of arbitrary length?


############################
### BULK COMPOSITION RELATED
############################
ag   = Constraint('name',"'%Ag-fcc%'",		'like')
al   = Constraint('name',"'%Al-fcc%'",		'like')
au   = Constraint('name',"'%Au-fcc%'",		'like')
ba   = Constraint('name',"'%Ba-bcc%'",		'like')
be   = Constraint('name',"'%Be-hcp%'",		'like')
c    = Constraint('name',"'%C-diamond%'",	'like')
ca   = Constraint('name',"'%Ca-fcc%'",		'like')
cd   = Constraint('name',"'%Cd-hcp%'",		'like')
co   = Constraint('name',"'%Co-hcp%'",		'like')
cu   = Constraint('name',"'%Cu-fcc%'",		'like')
fe   = Constraint('name',"'%Fe-bcc%'",		'like')
ge   = Constraint('name',"'%Ge-diamond%'",	'like')
ir   = Constraint('name',"'%Ir-fcc%'",		'like')
k    = Constraint('name',"'%K-bcc%'",		'like')
li   = Constraint('name',"'%Li-bcc%'",		'like')
mg   = Constraint('name',"'%Mg-hcp%'",		'like')
mo   = Constraint('name',"'%Mo-bcc%'",		'like')
na   = Constraint('name',"'%Na-bcc%'",		'like')
nb   = Constraint('name',"'%Nb-bcc%'",		'like')
ni   = Constraint('name',"'%Ni-fcc%'",		'like')
os   = Constraint('name',"'%Os-hcp%'",		'like')
pb   = Constraint('name',"'%Pb-fcc%'",		'like')
pd   = Constraint('name',"'%Pd-fcc%'",		'like')
rb   = Constraint('name',"'%Rb-bcc%'",		'like')
rh   = Constraint('name',"'%Rh-fcc%'",		'like')
ru   = Constraint('name',"'%Ru-hcp%'",		'like')
sc   = Constraint('name',"'%Sc-hcp%'",		'like')
si   = Constraint('name',"'%Si-diamond%'",	'like')
sn   = Constraint('name',"'%Sn-diamond%'",	'like')
sr   = Constraint('name',"'%Sr-fcc%'",		'like')
ti   = Constraint('name',"'%Ti-hcp%'",		'like')
zn   = Constraint('name',"'%Zn-hcp%'",		'like')
zr   = Constraint('name',"'%Zr-hcp%'",		'like')

elements = [ag,al,au,ba,be,c,ca,cd,co,cu,fe,ge,ir,k,li,mg,mo,na,nb,ni
			,os,pb,pd,rb,rh,ru,sc,si,sn,sr,ti,zn,zr]

singleElements  = Constraint('name',[x.val for x in elements], 			'likelist')
variousElements	= Constraint('name',[x.val for x in [li,be,al,cu,c]],	'likelist')
liAndBe			= Constraint('name',[x.val for x in [li,be]], 			'likelist')


badElements  = Constraint('name',[x.val for x in [fe,au,ni,ca,ag,os,c,ru,fe,mg,si,ba,zr,pd,co,sr,ti,rb,cd,ge,k,mo,sc]],'notlikelist')
goodElements = Constraint('name',[x.val for x in [fe,au,ni,ca,ag,os,c,ru,fe,mg,si,ba,zr,pd,co,sr,ti,rb,cd,ge,k,mo,sc]],'likelist')

####################
# SIMPLE CONSTRAINTS
####################

########
# Equals
########
hcp 	= Constraint('structure',"'hexagonal'", '=')
fcc 	= Constraint('structure',"'fcc'", 		'=')
bcc 	= Constraint('structure',"'bcc'", 		'=')
diamond = Constraint('structure',"'diamond'",	'=')

pbe  	= Constraint('xc',         "'PBE'",   	'=')
mbeef 	= Constraint('xc',         "'mBEEF'", 	'=')
rpbe 	= Constraint('xc',         "'RPBE'", 	'=')

gpaw 	= Constraint('dftcode',"'gpaw'",			'=')
qe 		= Constraint('dftcode',"'quantumespresso'",	'=')

pbePrecalc 	= Constraint('precalc',	"'PBE'",	'=')
noPrecalc 	= Constraint('precalc',	"'None'",	'=')

sg15 		= Constraint('psp',"'sg15'",		'=')
gbrv15pbe 	= Constraint('psp',"'gbrv15pbe'",	'=')

f111		= Constraint('facet',"'[1,1,1]'",'=')
f100 		= Constraint('facet',"'[1,0,0]'",'=')

symmetric 	= Constraint('symmetric',"'True'",	'=')
asymmetric 	= Constraint('symmetric',"'False'",	'=')
novacancy 	= Constraint('vacancies',"null",	'=') #how should this work?
noadsorbate = Constraint('adsorbates',"null",	'=')

def status(x): 	return Constraint('status',"'"+x+"'",'=')
def pw(x): 		return Constraint('pw',x,'=')
def kptden(x): 	return Constraint('kptden',x,'=')
def layers(x): 	return Constraint('layers',x,'=')
def surfXY(x,y):return Constraint('xy',"'[%d,%d]'"%(x,y),'=')
def vacuum(x):	return Constraint('vacuum',x,'=')

initialized = status('initialized')
completed 	= status('complete')
fizzled  	= status('fizzled')
queued      = status('queued')

######
# LIKE
#######
bn   = Constraint('name', "'%BN-zincblende%'",	'like')
licl = Constraint('name', "'%LiCl%'",			'like')
li 	 = Constraint('name',"'%Li-bcc%'", 			'like')

##########
# < > > <
##########

bfitok 		= Constraint('bfit', 0.5,	'>')
pwOver1000 	= Constraint('pw', 	1000,	'>')
pwUnder1000 = Constraint('pw', 	1000,	'<')

###########
# Negations
###########
incompleted = Constraint('status',"'complete'",'!=')

notMBeef 	= Constraint('xc', 		"'mBEEF'", 			'!=')
notBeef    	= Constraint('xc', 		"'BEEF'", 			'!=')	
notGPAW 	= Constraint('dftcode',"'gpaw'",			'!=')
notQE 		= Constraint('dftcode',"'quantumespresso'",	'!=')

###############
# COMPOSITE
###############

gpawsg15 	= Constraint('',(gpaw,sg15), 	'and')
qegbrv15pbe = Constraint('',(qe,gbrv15pbe),	'and')

pbe_mbeef_calcs 	= Constraint('',(mbeef,pbePrecalc),		'and')
noPrecalc_notMBEEF 	= Constraint('',(notMBeef,noPrecalc),	'and')

pspGuard 			= Constraint('',(gpawsg15,qegbrv15pbe),					'or')
preCalcWhenNeeded 	= Constraint('',(pbe_mbeef_calcs,noPrecalc_notMBEEF),	'or')
noBEEFinGPAW 		= Constraint('',(notGPAW,notBeef),	'or')
noMBEEFinQE 		= Constraint('',(notQE,notMBeef),	'or')

######################
# LISTS OF CONSTRAINTS
######################
singleElementFilters 	= [[x] for x in elements]
singleElementLabels	 	= [x.val.split('%')[1] for x in elements]
testFilters 			= [initialized,singleElements,pwUnder1000,kptden(2)]

essentialBulkJobConstraints = 	[pspGuard
								,preCalcWhenNeeded
								,noBEEFinGPAW
								,noMBEEFinQE]   # change with caution

essentialSurfJobConstraints = []

