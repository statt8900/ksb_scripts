import numpy as np

from db 				import plotQuery,addCol,defaulttable,tentativetable,updateDB,allJobIDs,dropTable,sqlexecute,countDB,query1
from job 				import cell2param,db2object
from createDB 			import detailstable
from constraint 		import *
from data_solids_wPBE 	import data,getKey
from printParse 		import parseNp

from ase.db 			import connect
from ase.data 			import chemical_symbols

from sqlite3 			import OperationalError
from ast 				import literal_eval


"""
A Detail is some derived piece of data from a job that we are interested in. 
It has an associated name (sql column header) and function (jobID -> a)
Some Details rely on other Detail columns existing, so they must be ordered by some Rank (rank 1 = no dependence)
"""
asedb		= connect('/scratch/users/ksb/db/ase.db')
nELEM		= len(chemical_symbols)
toAtomicNum = {symb:num for symb,num in zip(chemical_symbols,range(nELEM))}
nonmetals 	= ['H','C','N','P','O','S','Se','F','Cl','Br','I','At','He','Ne','Ar','Kr','Xe','Rn']



#############################################################################
# Irreducible information from job LEFT JOIN convergence LEFT JOIN result:  #
#############################################################################
# 	jobkind,aseidinitial,vibids,nebids,xc 									#
# 		pw,kptden,psp,xtol,strain,convid,precalc,dftcode 					#
#		comments,error,status 												#
# 	dwrat,econv,mixing,nmix,maxstep,nbands,sigma,fmax,delta,climb 			#
#	launchdir,aseid,energy,forces,bulkmodulus,bfit 							#
#		xccoeffs,viblist,dos,barrier,time,niter 							#
#############################################################################


class Detail(object):
	def __init__(self,name,kind,rank,cols,func,axislabel=None,convtol=None,derivtol=None):
		self.name 		= name
		self.kind 		= kind 
		self.rank 		= rank
		self.cols 		= cols
		self.func 		= func
		self.axislabel	= axislabel
		self.convtol 	= convtol
		self.derivtol 	= derivtol

	def apply(self,table): return [self.func(x) for x in plotQuery(self.cols,[],table)]


#####################
# Auxillary Functions
#####################
def identity(x): 	return x
def unbox(x): 		return x[0]
def maybe(f,x): 	return None if x is None else f(x) 
def intersect(listx,listy): 
	return list(set(listx) & set(listy))

def keldData(k): 	
	"""Presumes only col requested is 'name'. Accesses key of data entry in Keld's dataset"""
	def getKeldData(name):
		try: 
			return data[getKey[str(name[0])]][k]
		except KeyError: return None
	return getKeldData

def getFromAse(key,f=identity): 	
	return lambda aseid: None if aseid[0] is None else f(asedb.get(id=aseid[0]).get(key))

def fromJobObject(f):
	global jobtable
	def getFromJobObject(jobid):
		jobObject = db2object(jobid[0],jobtable)
		return f(jobObject)
	return getFromJobObject

def diff(col1col2): return None if None in col1col2 else float(col1col2[0])-float(col1col2[1])

def getJobCalcConstraints(jobid):
	"""Constraints that isolate all jobs done with same calculator/convergence parameters"""
	xc,pw,kptden,psp,xtol,strain,convid,precalc,dftcode = plotQuery(['xc','pw','kptden','psp','xtol','strain','convid','precalc','dftcode'],[JOBID(jobid)])[0]
	return [XC(xc),PW(pw),KPTDEN(kptden),PSP(psp),XTOL(xtol),STRAIN(strain),CONVID(convid),PRECALC(precalc),DFTCODE(dftcode)]
def similarJobs(jobid): return [x[0] for x in plotQuery(['jobid'],getJobCalcConstraints(jobid[0]))]


################################################################################
########
# Rank 1
########

name 		= Detail('name'
					,'varchar',1,['aseidinitial']
					,getFromAse('name'))

relaxed 	= Detail('relaxed'
					,'varchar',1,['aseidinitial']
					,getFromAse('relaxed'))

symbols 	= Detail('symbols'
					,'varchar',1,['aseidinitial']
					,getFromAse('symbols',str))

emt 		= Detail('emt'
					,'varchar',1,['aseidinitial']
					,getFromAse('emt'))

comments	= Detail('comments'
					,'varchar',1,['aseidinitial']
					,getFromAse('comments'))

kind 		= Detail('kind'
					,'varchar',1,['aseidinitial']
					,getFromAse('kind'))

structure 	= Detail('structure'
					,'varchar',1,['aseidinitial']
					,getFromAse('structure'))

vacancies 	= Detail('vacancies'
					,'varchar',1,['aseidinitial']
					,getFromAse('vacancies'))

sites 		= Detail('sites'
					,'varchar',1,['aseidinitial']
					,getFromAse('sites'))

facet 		= Detail('facet'
					,'varchar',1,['aseidinitial']
					,getFromAse('facet'))

xy 			= Detail('xy'
					,'varchar',1,['aseidinitial']
					,getFromAse('xy'))

layers 		= Detail('layers'
					,'varchar',1,['aseidinitial']
					,getFromAse('layers')
					,'Number of layers')

constrained = Detail('constrained'
					,'varchar',1,['aseidinitial']
					,getFromAse('constrained')
					,'Number of constrained layers')

symmetric 	= Detail('symmetric'
					,'varchar',1,['aseidinitial']
					,getFromAse('symmetric'))

vacuum	 	= Detail('vacuum'
					,'varchar',1,['aseidinitial']
					,getFromAse('vacuum')
					,r'Vacuum, $\AA$')

adsorbates 	= Detail('adsorbates'
					,'varchar',1,['aseidinitial']
					,getFromAse('adsorbates'))


pw = Detail('pwcutoff'
			,'numeric',1,['pw']
			,unbox
			,'Planewave cutoff, eV')


kpt 	= Detail('kpt'
				,'varchar',1,['jobid']
				, fromJobObject(lambda x: str(x.kpt()))) 

kptden 	= Detail('kptdensity'
				,'numeric',1,['kptden']
				,unbox
				,r'K-point density, pts/$\AA^{-1}$, eV')

bulkmod = Detail('bulkmod_calc'
					,'numeric',1,['bulkmodulus']
					,unbox
					,'Expt Bulk Mod, GPa')
energy = Detail('energy_calc'
					,'numeric',1,['energy']
					,unbox
					,'Electronic Energy, eV')

timeperstep = Detail('timeperstep'
					,'numeric',1,['time']
					,unbox
					,'Time per ionic step, min')

vibs 		= Detail('vibs'
					,'varchar',1,['viblist']
					,unbox)


cell 		= Detail('cell'
					,'varchar',1,['aseid']
					,getFromAse('cell',str))


########
# Rank 2
########
calcLabel 	= Detail('calclabel'
					,'varchar',2,['dftcode','xc','pw','kpt']
					,lambda x: '%s_%s_%d_%s'%(x[0],x[1],x[2],x[3]))


params 	= Detail('params'
				,'numeric',2,['cell']
				,lambda x: maybe(lambda y: str(cell2param(parseNp(y))),x[0]))

exptA 	= Detail('a_expt'
				,'numeric',2,['name']
				,keldData('lattice parameter')
				,r'Expt Lattice A, $\AA$')

exptBM 	= Detail('bulkmod_expt'
				,'numeric',2,['name']
				,keldData('bulk modulus')
				,'Expt Bulk Mod, GPa')


metalStoich = Detail('metalstoich'
				,'varchar',2,['symbols']
				,lambda x: str([y for y in literal_eval(x[0]) if y not in nonmetals]))


"""
Default vib analysis:
	- metal atoms have zero vibration
	- energies,etc. get multiplied by (# vibrating atoms/(#number of nonmetal atoms)
	- print warning in cases like vibrating C,O when nonmetal atoms are C,O,C,O,H


"""

########
# Rank 3
########
errBM 		= Detail('errBM'
					,'numeric',3,['bulkmodulus','bulkmod_expt']
					,diff
					,'Error in BulkMod, GPa')

calcA 		= Detail('a'
					,'numeric',3,['params']
					,lambda x:maybe(lambda y: str(parseNp(y).item(0)),x[0])
					,r'Calculated Lattice A, $\AA$')

########
# Rank 4
########
def errAFunc(a_aExpt_name):
	a,aexpt,name = a_aExpt_name
	if a is None or aexpt is None: return None
	elif 'hcp' in name: multFactor = 1
	elif 'bcc' in name: multFactor = 3**(0.5)/2.
	else:				multFactor = 2**(-0.5)
	return a - multFactor*aexpt

errA 	= Detail('errA'
				,'numeric',4,['a','a_expt','name']
				,errAFunc
				,r'Error in Lattice A, $\AA$')

###################################################
# DEPENDENT ON REFERENCE (defined in constraint.py)
###################################################
def referenceIDs(): return  [x[0] for x in plotQuery(['jobid'],REF_CONSTRAINTS)] #THE REFERENCES MUST ALREADY HAVE RANK 1 RESULTS - need to run twice if a new references is added?

isRef 		= Detail('isref'
					,'integer',2,['name']
					,lambda jobid: int(jobid in [x[0] for x in plotQuery(['jobid'],REF_CONSTRAINTS)]))

def getStoichVec(jobid):
	stoich = literal_eval(plotQuery(['symbols'],[JOBID(jobid)])[0][0])
	stoichNum = [toAtomicNum[x] for x in stoich]
	output = [0] * nELEM
	for n in stoichNum: output[n-1]+=1
	return np.array(output)

def getReferenceMatrix(refIDs,keyword):
	mat = np.transpose(np.matrix([getStoichVec(j) for j in refIDs]))	
	vec = np.array([query1(keyword,'jobid',x) for x in refIDs])
	return mat,vec

def getRefEnergy(jobid):
	stoichVec 		= getStoichVec(jobid[0])
	refIDs 			= intersect(similarJobs(jobid),referenceIDs())
	refMat, refVec 	= getReferenceMatrix(refIDs,'energy')
	contribs	 	= np.linalg.lstsq(refMat,stoichVec)	#refMat * x = stoichVec
	def filterZeros(x,y): return [(x,y) for x,y in zip(x,y) if x!=0 and y!=0]
	def dot(x,y): 		return sum([x*y for x,y in zip(x,y)])
	def safeDot(x,y): 	return None if x.count(None)+y.count(None) > 0 else dot(x,y)
	filtered = filterZeros(contribs[0],refVec)
	output= safeDot(*zip(*filtered)) # WHY DOES LINALG RETURN TUPLE??? maybe should be try: except block
	return output

"""This Detail deceptively depends on 'energy' and calculator parameters"""
refEnergy  = Detail('refeng'
					,'numeric',2,['jobid'] 
					,getRefEnergy)

eFormation = Detail('eform'
					,'numeric',3,['energy','refeng']
					,diff,r'$\Delta H_{form},eV$')

###########################################################################
### KEEP ALLDETAILS UP TO DATE
# abandoned details: 'isref'
allDetails = [name,relaxed,symbols,emt,comments,kind,structure,vacancies,sites,facet
				,xy,layers,constrained,symmetric,vacuum,adsorbates,cell,params
				,exptA,exptBM,calcA,errBM,errA,pw,kptden,kpt,energy,calcLabel,bulkmod,timeperstep
				,isRef,refEnergy,eFormation,metalStoich]

activeDetails = allDetails # allDetails #or select a subset [exptA,exptBM,errBM,errA]

activeDetails.sort(key=lambda x: x.rank)

###############################################################################
def applyDetails(tentative=False,overwrite=False):	
	"""To-do: Time doing overwrite vs checking and overwriting if new """
	global jobtable #make global so that functions like fromJobObject() can access correct one
	jobtable 	= 'tentativejob' 		if tentative else 'job'
	querytable 	= tentativetable 		if tentative else defaulttable
		
	newjobs =  [i[0]for i in plotQuery(['jobid'],[],jobtable) if i[0] > countDB('details')]
	for i in newjobs:
		sqlexecute('INSERT INTO details (jobtableid) VALUES (%d)'%i)
	
	for i,d in enumerate(activeDetails):
		print d.name, " %d/%d"%(i+1,len(activeDetails))

		try:  addCol(d.name,d.kind,'details')
		except OperationalError: pass
		
		vals = d.apply(querytable)
		ids	 = allJobIDs()
		for ID,val in zip(ids,vals):
			if overwrite or (val is not None and query1(d.name,'jobtableid',ID,'details') is None):
				updateDB(d.name,'jobtableid',ID,val,None,'details')

if __name__=='__main__':
	applyDetails()
