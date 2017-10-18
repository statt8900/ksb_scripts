from bulk    import BulkTraj,BulkPure
from surface import SurfaceTraj,SurfacePure
from calc    import Calc
from traj    import TrajDataBulk,TrajDataSurf
from job     import BulkJob,SurfaceJob

from miniDatabase import * #really we just want all the domains

from itertools import product
from copy import deepcopy

"""
This module contains the Filter object class, functions for working filters, and lots of Filter instances.
"""

###############################################
# FILTER HELPER FUNCTIONS
###############################################
boolDomain = [True,False]

def filterBulk(filterList):
	"""
	[Filter] -> [BulkJob], applies filters to globally-defined domains, returning all jobs that remain
	"""	
	filters = andFilters([noFilterBulk]+[f.filters for f in filterList]) # noFilter guarantees that all dictionary keys present
	filteredBulk     = [b     for b     in bulkDomain         if filters['bulk'](b)]
	filteredDFT      = [dft   for dft   in dftDomain          if filters['dft'](dft)]
	filteredXC       = [xc    for xc    in xcDomain           if filters['xc'](xc)]
	filteredPW       = [pw    for pw    in pwDomain           if filters['pw'](pw)]
	filteredKpt      = [kpt   for kpt   in kptDomain          if filters['kpt'](kpt)]
	filteredEconv    = [eConv for eConv in eConvDomain        if filters['eConv'](eConv)]
	filteredMixing   = [mix   for mix   in mixDomain          if filters['mix'](mix)]
	filteredNmix     = [nmix  for nmix  in nMixDomain         if filters['nmix'](nmix)]
	filteredMaxSteps = [ms    for ms    in maxStepDomain      if filters['maxstep'](ms)]
	filteredNbands   = [nb    for nb    in nBandDomain        if filters['nbands'](nb)]
	filteredSigma    = [sig   for sig   in sigmaDomain        if filters['sigma'](sig)]
	filteredXtol     = [x     for x     in xTolDomain         if filters['xtol'](x)]
	filteredMag      = [m     for m     in magDomain          if filters['mag'](m)]
	filteredPSP      = [p     for p     in pspDomain          if filters['psp'](p)]

	assert len(filteredDFT)==1, 'Only one DFT code must be specified at a time'
	
	if not filteredBulk:     print "No valid bulk stoich/structures"
	if not filteredDFT:      print "No valid DFT"
	if not filteredXC:       print "No valid XC"
	if not filteredPW:       print "No valid PW"
	if not filteredKpt:      print "No valid KPT"
	if not filteredEconv:    print "No valid Econv"
	if not filteredMixing:   print "No valid Mixing "
	if not filteredNmix:     print "No valid Nmix "
	if not filteredMaxSteps: print "No valid Maxsteps"
	if not filteredNbands:   print "No valid Nbands"
	if not filteredSigma:    print "No valid Sigma"
	if not filteredXtol:     print "No valid Xtol"
	if not filteredMag:      print "No valid Mag"
	if not filteredPSP:      print "No valid PSP"

	combinations = product(filteredBulk,filteredDFT,filteredXC,filteredPW,filteredKpt,filteredEconv
							,filteredMixing,filteredNmix,filteredMaxSteps,filteredNbands,filteredSigma
							,filteredXtol,filteredMag,filteredPSP)
	def makeJob(job):
		if isinstance(job[0],TrajDataBulk): 
			bulk = BulkTraj(job[0].name)
		else: 
			bulk=BulkPure(job[0])
		dft = job[1]
		calc = Calc(*job[2:])
		return BulkJob(bulk,calc,dft)
	return map(makeJob,combinations)

def filterSurf(bulkFilters,surfFilters):
	"""
	[Filter] -> [SurfJob], applies filters to globally-defined domains, returning all jobs that remain
	"""
	
	bulkJobs =  filterBulk(bulkFilters) 
	for j in bulkJobs: 
		if not j.checkResult():
			print 'Bulk structure not yet completed; cannot create surface for '+str(j)
	bulkJobs = [j for j in bulkJobs if j.checkResult()]	
	
	filters       = andFilters([noFilterSurf]+surfFilters) # noFilter guarantees that all dictionary keys present
	
	filteredFacet = [f     for f     in facetDomain        if filters['facet'](f)]
	filteredScale = [s     for s     in scaleDomain        if filters['scale'](s)]
	filteredSym   = [sym   for sym   in boolDomain         if filters['sym'](sym)]
	filteredFixed = [fixed for fixed in fixedDomain        if filters['fixed'](fixed)]
	filteredVac   = [vac   for vac   in vacDomain          if filters['vac'](vac)]
	filteredAds   = [ads   for ads   in adsDomain          if filters['ads'](ads)]

	filteredDFT      = [dft   for dft   in dftDomain          if filters['dft'](dft)]
	filteredXC       = [xc    for xc    in xcDomain           if filters['xc'](xc)]
	filteredPW       = [pw    for pw    in pwDomain           if filters['pw'](pw)]
	filteredKpt      = [kpt   for kpt   in kptDomain          if filters['kpt'](kpt)]
	filteredEconv    = [eConv for eConv in eConvDomain        if filters['eConv'](eConv)]
	filteredMixing   = [mix   for mix   in mixDomain          if filters['mix'](mix)]
	filteredNmix     = [nmix  for nmix  in nMixDomain         if filters['nmix'](nmix)]
	filteredMaxSteps = [ms    for ms    in maxStepDomain      if filters['maxstep'](ms)]
	filteredNbands   = [nb    for nb    in nBandDomain        if filters['nbands'](nb)]
	filteredSigma    = [sig   for sig   in sigmaDomain        if filters['sigma'](sig)]
	filteredXtol     = [x     for x     in xTolDomain         if filters['xtol'](x)]
	filteredMag      = [m     for m     in magDomain          if filters['mag'](m)]
	filteredPSP      = [p     for p     in pspDomain          if filters['psp'](p)]
	
	filteredTraj     = [t     for t     in sTrajDomain        if filters['traj'](x)]
	
	trajCombinations = product(filteredTraj,filteredScale,filteredAds,
				filteredDFT,filteredXC,filteredPW,filteredKpt,filteredEconv,filteredMixing,filteredNmix,filteredMaxSteps,filteredNbands,filteredSigma,filteredXtol,filteredMag,filteredPSP)

	pureCombinations = product(bulkJobs,filteredFacet,filteredScale,filteredSym,filteredFixed,filteredVac,filteredAds,
				filteredDFT,filteredXC,filteredPW,filteredKpt,filteredEconv,filteredMixing,filteredNmix,filteredMaxSteps,filteredNbands,filteredSigma,filteredXtol,filteredMag,filteredPSP)
	
	def makePure(job):
		surf  = SurfacePure(*job[:7])
		dft   = job[7]
		calc  = Calc(*job[8:])
		return SurfaceJob(surf,calc,dft)

	def makeTraj(job):
		surf = SurfaceTraj(job[0].name,job[1],job[2])
		dft  = job[3]
		calc = Calc(*job[4:])
		return SurfaceJob(surf,calc,dft)
	return map(makeTraj,trajCombinations)+map(makePure,pureCombinations)
############################################################################################
############################################################################################
# BEGIN COLLECTING FILTERS / FILTERTOOLS TO IMPORT TO BATCHJOBS
############################################################################################
__all__ = ['filterBulk','filterSurf']


class Filter(object):
	def __init__(self,name,filters):
		assert isinstance(name,str)
		assert isinstance(filters,list)
		assert isinstance(filters[0],dict)
		self.name    = name                # String
		self.filters = andFilters(filters) # {keyword : Job.keyword -> Bool}



############################################################################################
# PRIMITIVE FILTER FUNCTIONS
############################################################################################
def true(y):   return True                   #trivially true function
def eq(y):     return lambda x: x == y
def neq(y):    return lambda x: x != y
def fIN(y):    return lambda x: x in y
def fNIN(y):   return lambda x: x not in y
def gt(y):     return lambda x: x > y
def lt(y):     return lambda x: x < y
def fstEQ(y):  return lambda x: x[0] == y
def fstIN(y):  return lambda x: x[0] in y
def fstGT(y):  return lambda x: x[0] > y
def fstLT(y):  return lambda x: x[0] < y
def fAND(x,y): return lambda z: x(z) and y(z)
def fOR(x,y):  return lambda z: x(z) or y(z)
def fINST(y):  return lambda x: isinstance(x,y)

__all__ += ['true','eq','neq','fIN','fNIN','gt','lt','fstEQ','fstIN','fstGT','fstLT','fAND','fOR','fINST']
############################################################################################
# FILTER COMBINATORS
############################################################################################
def merge(A, B, f):
	"""
	Combine two dictionaries, applying the rule 'f' when the two share a key
	"""
	merged = {k: A.get(k, B.get(k)) for k in A.viewkeys() ^ B.viewkeys()}  # Start with symmetric difference; keys either in A or B
	merged.update({k: f(A[k], B[k]) for k in A.viewkeys() & B.viewkeys()}) # Update with `f()` applied to the intersection
	return merged

def andFilters(fList): 
	"""
	Merge two filter dictionaries, combining with logical AND
	"""
	outFilter = fList[0]
	for f in fList[1:]:	outFilter = merge(f,outFilter,fAND)
	return outFilter

def orFilters(fList): 
	"""
	Merge two filter dictionaries, combining with logical OR
	"""
	outFilter = fList[0]
	for f in fList[1:]:	outFilter = merge(f,outFilter,fOR)
	return outFilter


def andFilterObjects(name,filterObjectList): 
	assert isinstance(name,str)
	assert isinstance(filterObjectList,list)
	assert isinstance(filterObjectList[0],Filter)
	return Filter(name,[deepcopy(f.filters) for f in filterObjectList])

def orFilterObjects(name,filterObjectList):   return Filter(name,[orFilters([deepcopy(f.filters) for f in filterObjectList])])

__all__ +=['andFilters','orFilters','andFilterObjects','orFilterObjects']

#############################################################################################
# TRIVIAL FILTERS
#############################################################################################
noFilterBulk =	{'bulk':    true
				,'dft':     true
				,'xc':	    true
				,'pw':	    true
				,'kpt':     true
				,'eConv':   true
				,'mix':     true
				,'nmix':    true
				,'maxstep': true
				,'nbands':  true
				,'sigma':   true
				,'xtol':    true
				,'mag':     true
				,'psp':     true}

noFilterSurf =	{'dft':     true
				,'xc':	    true
				,'pw':	    true
				,'kpt':     true
				,'eConv':   true
				,'mix':     true
				,'nmix':    true
				,'maxstep': true
				,'nbands':  true
				,'sigma':   true
				,'xtol':    true
				,'facet':   true
				,'scale':   true
				,'sym'  :   true
				,'fixed':   true
				,'vac'  :   true
				,'ads'  :   true
				,'traj'  :   true}

__all__ +=['noFilterBulk','noFilterSurf']

##################################################################################################
# Bulk
##################################################################################################
hcp         = Filter('hcp',        [{'bulk': fIN(hcps)}])
bcc         = Filter('bcc',        [{'bulk': fIN(bccs)}])
fcc         = Filter('fcc',        [{'bulk': fIN(fccs)}])
diamond     = Filter('diamond',    [{'bulk': fIN(diamonds)}])
pureMetal   = Filter('pureMetal',  [{'bulk': fIN(pureMetals)}])
cscl        = Filter('cscl',       [{'bulk': fIN(cscls)}])
rockSalt    = Filter('rockSalt',   [{'bulk': fIN(rocksalts)}])
zincBlende  = Filter('zincBlende', [{'bulk': fIN(zincblendes)}])
binaryAlloy = Filter('binaryAlloy',[{'bulk': fIN(binaryAlloys)}])


bulkPure   = Filter('FromPure',[{'bulk': fINST(tuple)}])
bulkTraj   = Filter('FromTraj',[{'bulk': fINST(TrajDataBulk)}])


gpaw     = Filter('gpaw',  [{'dft': eq('gpaw')}
							,{'xc':fNIN(['beef'])}
							,{'psp': neq('gbrv15pbe')}])

qe       = Filter('quantumEspresso',[{'dft': eq('quantumEspresso')}
									,{'psp': neq('sgs15')}])


pbe          = Filter('PBE',  [{'xc': eq('PBE')}])
beef         = Filter('BEEF', [{'xc': eq('BEEF')}])
mbeef        = Filter('mBEEF',[{'xc': eq('mBEEF')}])

pwDomain      = [500,700,900,1100,1300,1500,1700]

pw500        = Filter('pw500', [{'pw': eq(500)}])
pw700        = Filter('pw700', [{'pw': eq(700)}])
pw900        = Filter('pw900', [{'pw': eq(900)}])
pw1100       = Filter('pw1100',[{'pw': eq(1100)}])
pw1300       = Filter('pw1300',[{'pw': eq(1300)}])
pw1500       = Filter('pw1500',[{'pw': eq(1500)}])
pw1700       = Filter('pw1700',[{'pw': eq(1700)}])

kpt555       = Filter('kpt (5,5,5)',   [{'kpt': eq((5,5,5))}])
kpt111111    = Filter('kpt (11,11,11)',[{'kpt': eq((11,11,11))}])
kpt151515    = Filter('kpt (15,15,15)',[{'kpt': eq((15,15,15))}])
isoMetricKpt = andFilterObjects('kpt (x,x,x)',[kpt555,kpt111111,kpt151515])

kpt441       = Filter('kpt (4,4,1)', [{'kpt': eq((4,4,1))}])
surfKpt      = andFilterObjects('kpt (x,x,1)',[kpt441])

eConv5       = Filter('eConv = 5e-5 eV',[{'eConv':eq(5)}])
mix5         = Filter('mixing = 0.05',  [{'mix':eq(5)}])
mix10        = Filter('mixing = 0.1',   [{'mix':eq(10)}])
nMix5        = Filter('nMix = 5',       [{'nmix':eq(5)}])
sigma1       = Filter('sigma = 0.1',    [{'sigma':eq(1)}])
nband12      = Filter('12 extra bands', [{'nband':eq(12)}])
maxStep500   = Filter('maxstep = 500',  [{'maxstep':eq(500)}])

xtol5        = Filter('xTol 0.5 pm',[{'xtol': eq(5)}])
xtol20       = Filter('xTol 2 pm',  [{'xtol': eq(20)}])

mag          = Filter('Magnetic materials',   [{'mag': neq(0),'bulk': fIN(magList)}])
nonmag       = Filter('Nonmagnetic materials',[{'mag':  eq(0),'bulk': fNIN(magList)}])

sg15         = Filter('sg15 setup',[{'psp': eq('sgs15')}])
gbrv         = Filter('gbrv psp',  [{'psp': eq('gbrv15pbe')}])

__all__ +=['hcp','bcc','fcc','diamond','pureMetal','cscl','rockSalt','zincBlende','binaryAlloy'
			,'bulkPure','bulkTraj'
			,'gpaw','qe'
			,'pbe','beef','mbeef'
			,'pw500','pw700','pw900','pw1100','pw1300','pw1500','pw1700'
			,'kpt555','kpt111111','kpt151515','isoMetricKpt'
			,'kpt441','surfKpt'
			,'eConv5','mix5','mix10','nMix5','sigma1','nband12','maxStep500'
			,'xtol5','xtol20'
			,'mag','nonmag'
			,'sg15','gbrv']

#########
# Surface
#########
f001        = Filter('100 Facet',[{'facet': eq((0,0,1))}])
f111        = Filter('111 Facet',[{'facet': eq((1,1,1))}])
f100        = Filter('100 Facet',[{'facet': eq((1,0,0))}])

scale114    = Filter('1x1x4',[{'scale': eq((1,1,4))}])
scale224    = Filter('2x2x4',[{'scale': eq((2,2,4))}])

sym         = Filter('Symmetric',[{'sym':eq(True)}])
assym       = Filter('Asymmetric',[{'sym':eq(False)}])

fixed2      = Filter('2 layers fixed',[{'fixed':eq(2)}])

vac9        = Filter('vacuum = 9 A',[{'vac': eq(9)}])

bare        = Filter('Bare',[{'ads': eq({})}])
notBare     = Filter('not Bare',[{'ads': neq({})}])

__all__ +=  ['f001','f111','f100'
			,'scale114','scale224'
			,'sym','assym'
			,'fixed2','vac9'
			,'bare','notBare']

##############################
# COMMON COMBINATIONS
##############################
defaultCalc = andFilterObjects('Default Calc Misc',[eConv5,mix5,nMix5,maxStep500,xtol5,sigma1])

gpawsg15    = andFilterObjects('GPAW - sg15',[gpaw,sg15])
qegbrv      = andFilterObjects('Quantum Espresso - gbrv',[qe,gbrv])

highAcc     = andFilterObjects('pw = 1500 eV, kpt = (15,15,15)',[pw1500,kpt151515])
lowAcc      = andFilterObjects('pw = 500 eV, kpt (5,5,5)',[pw500,kpt555])
__all__ +=['defaultCalc','gpawsg15','qegbrv'
			,'highAcc','lowAcc']
###############################################
# Stoich specific FILTERS
###############################################
liBcc    = Filter('Li bcc',   [{'bulk': eq(('Li2','bcc'))}])
alFcc    = Filter('Al fcc',   [{'bulk': eq(('Al','fcc'))}])
cDiamond = Filter('C diamond',[{'bulk':eq(('C2','diamond'))}])
beHcp    = Filter('Be hcp',   [{'bulk':eq(('Be2','hcp'))}])

li2Generic = andFilterObjects('Li bcc',   [liBcc,nonmag,defaultCalc])
alGeneric  = andFilterObjects('Al fcc',   [alFcc,nonmag,defaultCalc])
cGeneric   = andFilterObjects('C diamond',[cDiamond,nonmag,defaultCalc])
beGeneric  = andFilterObjects('Be hcp',   [beHcp,nonmag,defaultCalc])

__all__ +=['liBcc','alFcc','cDiamond','beHcp'
			,'li2Generic','alGeneric','cGeneric','beGeneric']

###############################################
# Other FILTERS
###############################################
normalPW = orFilterObjects('',[pw900,pw1100,pw1300,pw1500])
normalKPT = orFilterObjects('',[kpt111111])
johannesBulkFilter = andFilterObjects('Johannes',[defaultCalc,gpawsg15,mbeef,normalPW,normalKPT,bulkPure])

variousMetals = orFilterObjects('Various metals',[li2Generic,alGeneric,cGeneric,beGeneric])

__all__ +=['normalPW','normalKPT'
			,'johannesBulkFilter','variousMetals']

##############

pureLowAcc     = andFilterObjects('All materials\nlow acc', [bulkPure,lowAcc])
pureHighAcc    = andFilterObjects('All materials\nhigh acc',[bulkPure,highAcc])
pureHCP        = andFilterObjects('hcp materials\nhigh acc',[bulkPure,hcp,highAcc])
pureFCC        = andFilterObjects('fcc materials\nhigh acc',[bulkPure,fcc,highAcc])
pureBCC        = andFilterObjects('bcc materials\nhigh acc',[bulkPure,bcc,highAcc])

__all__ +=['pureHighAcc','pureLowAcc'
			,'pureHCP','pureFCC','pureBCC']

##############
LiBccPW500     = andFilterObjects('Li bcc, pw 500',  [liBcc,pw500])
LiBccPW700     = andFilterObjects('Li bcc, pw 700',  [liBcc,pw700])
LiBccPW900     = andFilterObjects('Li bcc, pw 900',  [liBcc,pw900])
LiBccPW1100    = andFilterObjects('Li bcc, pw 1100', [liBcc,pw1100])
LiBccPW1300    = andFilterObjects('Li bcc, pw 1300', [liBcc,pw1300])
LiBccPW1500    = andFilterObjects('Li bcc, pw 1500', [liBcc,pw1500])
LiBccPW1700    = andFilterObjects('Li bcc, pw 1700', [liBcc,pw1700])
LiBccKPT555    = andFilterObjects('Li bcc, kpt=(5,5,5)',   [liBcc,kpt555])
LiBccKPT111111 = andFilterObjects('Li bcc, kpt=(11,11,11)',[liBcc,kpt111111])
LiBccKPT151515 = andFilterObjects('Li bcc, kpt=(15,15,15)',[liBcc,kpt151515])

li555PBE   =andFilterObjects('Li, Kpt 5x5x5',   [LiBccKPT555,pbe])
li111111PBE=andFilterObjects('Li, Kpt 11x11x11',[LiBccKPT111111,pbe])
li151515PBE=andFilterObjects('Li, Kpt 15x15x15',[LiBccKPT151515,pbe])
li500PBE   =andFilterObjects('Li, pw500', [LiBccPW500,pbe])
li700PBE   =andFilterObjects('Li, pw700', [LiBccPW700,pbe])
li900PBE   =andFilterObjects('Li, pw900', [LiBccPW900,pbe])
li1100PBE  =andFilterObjects('Li, pw1100', [LiBccPW1100,pbe])
li1300PBE  =andFilterObjects('Li, pw1300', [LiBccPW1300,pbe])
li1500PBE  =andFilterObjects('Li, pw1500',[LiBccPW1500,pbe])
li1700PBE  =andFilterObjects('Li, pw1700',[LiBccPW1700,pbe])

li555MBEEF   =andFilterObjects('Li, Kpt 5x5x5 mBEEF',   [LiBccKPT555,mbeef])
li111111MBEEF=andFilterObjects('Li, Kpt 11x11x11 mBEEF',[LiBccKPT111111,mbeef])
li151515MBEEF=andFilterObjects('Li, Kpt 15x15x15 mBEEF',[LiBccKPT151515,mbeef])
li500MBEEF   =andFilterObjects('Li, pw500 mBEEF', [LiBccPW500,mbeef])
li700MBEEF   =andFilterObjects('Li, pw700 mBEEF', [LiBccPW700,mbeef])
li900MBEEF   =andFilterObjects('Li, pw900 mBEEF', [LiBccPW900,mbeef])
li1100MBEEF  =andFilterObjects('Li, pw1100 mBEEF', [LiBccPW1100,mbeef])
li1300MBEEF  =andFilterObjects('Li, pw1300 mBEEF', [LiBccPW1300,mbeef])
li1500MBEEF  =andFilterObjects('Li, pw1500 mBEEF',[LiBccPW1500,mbeef])
li1700MBEEF  =andFilterObjects('Li, pw1700 mBEEF',[LiBccPW1700,pbe])


__all__ +=['LiBccPW500','LiBccPW700','LiBccPW900','LiBccPW1100','LiBccPW1300','LiBccPW1500','LiBccPW1700'
			,'LiBccKPT555','LiBccKPT111111','LiBccKPT151515'
			,'li555PBE','li111111PBE','li151515PBE','li500PBE','li700PBE','li900PBE'
			,'li1100PBE','li1300PBE','li1500PBE','li1700PBE'
			,'li555MBEEF','li111111MBEEF','li151515MBEEF','li500MBEEF','li700MBEEF','li900MBEEF'
			,'li1100MBEEF','li1300MBEEF','li1500MBEEF','li1700MBEEF']
