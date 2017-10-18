import os,json,time
import numpy as np

from objectClass import Object
from bulk        import BulkTraj,db2Bulk
from calc        import db2Calc
from db          import query1all
from printParse  import multiLine,printTime
"""
"""

def safeF(f,x): return None if x is  None else f(x)

class BulkJob(Object):

	def __init__(self,ID,created,createdBy,lastModified,kind,comments,status,bulkpth,calcid,precalcxc,dftCode,timeLim):
		"""Creates a bulk optimization job."""
		self.ID           = safeF(int,ID)            # Int
		self.created      = safeF(float,created)       # Double
		self.createdBy    = safeF(str,createdBy)     # String
		self.lastModified = safeF(float,lastModified)  # Double
		self.kind         = safeF(str,kind)          # String
		self.comments     = safeF(str,comments)       # String
		self.bulkpth      = safeF(str,bulkpth)       # Bulk
		self.calcid       = safeF(int,calcid)        # Int
		self.precalcxc    = safeF(str,precalcxc)     # Int
		self.dftCode      = safeF(str,dftCode)       # String
		self.timeLim      = safeF(int,timeLim)       # Int
	
	def __str__(self):  return "Job#%d"%(self.ID)
	
	def path(self): return  '/scratch/users/ksb/db/jobs/'+str(self.ID)
	
	def submit(self): 
		os.chdir(self.path())
		print "Submitting "+str(self)+" from sherlock"
		os.system('./metarun.sh')

	def prepare(self):
		path    = self.path()
		try: os.makedirs(path)
		except OSError: pass
		mPath   = path + '/metarun.sh'
		rPath   = path + '/run.sh'
		oPath   = path + '/opt.py'
		if os.path.exists(mPath): os.remove(mPath)
		if os.path.exists(rPath): os.remove(rPath);
		if os.path.exists(oPath): os.remove(oPath);

		BulkTraj(*query1all('bulk','pth',self.bulkpth)).writeInit(path) # save initial structure to traj file
		with open (mPath,'w') as f: f.write(self.metarunFile())
		with open (oPath,'w') as f: f.write(self.optFile())
		with open (rPath,'w') as f: f.write(self.runFile())
		os.chmod(mPath,0o775); os.chmod(rPath,0o775)



	def metarunFile(self):
		#Contents of file that is run (for any cluster, any dftcode) in order to submit batch job
		return multiLine(	['#!/bin/bash'
							,'if [ ! -f result.json ]; then'	
							,('sbatch ' if self.dftCode == 'gpaw' else '')+'./run.sh'	# Don't repeat already completed jobs
							,'fi\n'])

	def runFile(self): 
		return (multiLine(	['#!/bin/bash'
							,sherlockHeader(self.timeLim) if self.dftCode == 'gpaw' else ''
							,"#Slurm parameters"
							,"NTASKS=`echo $SLURM_TASKS_PER_NODE|tr '(' ' '|awk '{print $1}'`"
							,"NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`"
							,'NCPU=`echo " $NTASKS * $NNODES " | bc`'
							,'echo "$NCPU"'
							,"#load gpaw-specific paths"
							,"source /scratch/users/ksb/gpaw/paths.bash"
							,'echo "$1"'
							,"#run parallel gpaw"
							,"mpirun -n $NCPU gpaw-python ./opt.py\n"])
				 if self.dftCode == 'gpaw' else 'sbatch ./opt.py')
 
	def optFile(self): #Contents of opt.py, which optmizes lattice parameters (and atomic coordinates if necessary)
		return multiLine(	['#!/scratch/PI/suncat/sw/bin/python'
							,sherlockHeader(self.timeLim) if self.dftCode=='quantumespresso' else ''
							,'import time,os,traceback'
							,'from auto           import getBulkEnergy,getBulkModulus,getXCcontribs,saveResults,rank'
							,'from db             import updateStatus,updateDB'
							,'from scipy.optimize import fmin'
							,'from ase.parallel   import parallel_function'
							,'try:'
							,'\n\t#Initialize'
							,'\t# time.sleep(60) # sleep 1 min before modifying DB in case other jobs are still being submitted'
							,'\tupdateStatus({0},"queued","running")'.format(self.ID)
							,'\tstart=time.time()'
							,"\ttry: os.remove('lattice_opt.log')"
							,"\texcept OSError: pass"
							,'\n\t#Energy Minimizing Function'
							,'\tlatticeParams=fmin(getBulkEnergy,{0},ftol=2,xtol={1},args=({2},))'.format(self.iGuess(),self.calc().xtol,self.ID)
							,'\n\t#Bulk Modulus Function'
							,'\tcell,nIter,pos,magmom,e0,b,bFit = getBulkModulus({0},latticeParams)'.format(self.ID)
							,'\n\t#XC contribution Function'
							,'\txcContribs = getXCcontribs({0}) #(reads inp.pw generated in BulkMod)'.format(self.ID)
							,'\n\t#Save results'
							,'\ttotalTime = (time.time()-start)/3600'
							,'\tsaveResults({0},cell.tolist(),pos.tolist(),magmom.tolist(),e0,b,bFit,xcContribs.tolist(),totalTime,nIter)'.format(self.ID)
							,'except Exception as e:'
							,'\tprint e'
							,'\ttb = traceback.format_exc()'
							,'\tif rank() == 0:'
							,'\t\tupdateStatus({0},"running","failed")'.format(self.ID)
							,"\t\tupdateDB('bulkjob','comments','id',{0},tb)".format(self.ID)])

	def iGuess(self): return db2Bulk(self.bulkpth).iGuess()
	def calc(self): return db2Calc(self.calcid)
	########################################################################
	########################################################################
	# SQL METHODS
	########################################################################
	def sqlTable(self):  return 'bulkjob'

	def sqlCols(self):   return ['created'
								,'createdby'
								,'lastmodified'
								,'kind'
								,'comments'
								,'status'
								,'bulk'
								,'calcid'
								,'precalcxc'
								,'dft'
								,'timelim']

	def sqlInsert(self): return [time.time()
								,'kris'
								,time.time()
								,'bulkrelax'
								,None
								,'initialized'
								,self.bulkpth
								,self.calcid
								,self.precalcxc
								,self.dftCode
								,self.timeLim]
		

	def sqlEq(self): return ['kind','bulk','calcid','precalcxc','dft','timelim']

##############################################################################
##############################################################################
# OTHER FUNCTIONS
##############################################################################
def db2BulkJob(i): return  BulkJob(*query1all('bulkjob','id',i))

def sherlockHeader(t): return multiLine(['#SBATCH -p iric'
										,'#SBATCH -x gpu-14-1,sh-20-35'
										,'#SBATCH --job-name=myjob'
										,'#SBATCH --output=myjob.out'
										,'#SBATCH --error=myjob.err'
										,'#SBATCH --time={0}:00'.format(printTime(t))
										,'#SBATCH --qos=iric'
										,'#SBATCH --nodes=1'
										,'#SBATCH --mem-per-cpu=4000'
										,'#SBATCH --mail-type=END,FAIL'
										,'#SBATCH  --mail-user=ksb@stanford.edu'
										,'#SBATCH --ntasks-per-node=16\n'])
















"""
	def nAtoms(self): return len(self.bulk.makeAtoms())
	def getStruct(self):   return self.bulk
	def magmom(self):      return self.calc.magmom > 0
	def checkResult(self): return self.getBulkResult() is not None

	def getBulkResult(self): 
		#		Searches for result in current cluster, returning a Result object if found. (Warning printed if not found)
		
		import json
		pth = self.jobPth(getCluster(2).name)+'/result.json'
		if os.path.exists(pth): 
			with open(pth,'r') as f: data = f.read()
			return parseBulkResult(data)
		else: 
			print 'No result yet for ',self
			return None

	def hasPreCalc(self): return self.preCalc is not None  # Bool

	def zeros(self): return np.zeros(self.nAtoms()) # for magmom of non-magnetic materials

	def set_initial_magnetic_moments(self,atoms):
		ms =  [(self.calc.magmom if (ele in magElems) else 0) for ele in self.getStruct().atomList()]
		if self.calc.magmom>0:
			atoms.set_initial_magnetic_moments(ms)


###########################
# Job length prediction
###########################
def getScaling(dft):
	#data=scrape() -- filter just for a particular dft code
	#fit data to  t = A*n_e^3 *pw^(3/2) *prod_kpt
	# maybe read summary file that can be generated on its own to prevent weird import problems
	if dft == 'quantumEspresso': return 5e-10
	else: return 5e-9

# #BANDS^2*#PW
# BANDS = nElectrons/2 + (nbandsInput) #### use cell volume!

def getNelectrons(struct,calc): return sum(map(lambda x: pspDict[calc.psp].nElec[x],struct.atomList()))

def bulkTimeFunc(bulk,calc,dft):
	guess =  getScaling(dft)*(getNelectrons(bulk,calc)**3) * ((calc.pw)**1.5) * (calc.kpt[0]*calc.kpt[1]*calc.kpt[2])
	return min(max(guess,0.2),20)
def surfTimeFunc(surf,calc,dft):
	guess = 0.1+ getScaling(dft)*(getNelectrons(surf,calc)**3) * ((calc.pw)**1.5) * (calc.kpt[0]*calc.kpt[1]*calc.kpt[2])
	return min(max(guess,0.5),20)

###########################
# Auxillary functions
###########################
def parseBulkJob(x):
	#print 'parsing ',x
	bulk,calc,dftCode = derscoreSep(2,x)
	return BulkJob(parseBulk(bulk),parseCalc(calc),dftCode)

def parseSurfaceJob(x):
	surf,calc,dftCode = derscoreSep(3,x)
	return SurfaceJob(parseSurface(surf),parseCalc(calc),dftCode)

def parseBulkResult(strJson):
	#for i, line in enumerate(strJson.split('\n')): print i,line #WHY IS THERE AN EXTRA BRACKET RANDOMLY AT EOF
	load = json.loads(strJson)
	return Result(dictToBulkJob(load['job']),load['cell']
				,load['pos'],load['magmom']
				,load['energy'],load['bulkModulus'],load['bFit']
				,load['xcCoeffs'],load['time'],load['nIter'])

def parseSurfResult(strJson):
	load = json.loads(strJson)
	return Result(dictToSurfJob(load['job']),load['cell']
				,load['pos'],load['magmom']
				,load['energy'],load['bulkModulus'],load['bFit']
				,load['xcCoeffs'],load['time'],load['nIter'])

def dictToCalc(d): return Calc(d['xc'],int(d['pw']),(d['kpt']),int(d['eConv']),int(d['mixing']),int(d['nmix'])
								,int(d['maxsteps']),int(d['nbands']),int(d['sigma']),int(d['xtol']),int(d['magmom']),d['psp'])

def dictToBulkJob(dictJob): 
	return BulkJob(dictToBulk(dictJob['bulk']),dictToCalc(dictJob['calc']),dictJob['dftCode'],dictJob['timeLim'])

def dictToPureSurf(d):
	 return SurfacePure(dictToBulkJob(d['bulkJob']),d['facet'],d['scale'],d['symmetric'],d['constrained'],d['vacuum'],d['adsorbates'])

def dictToSurfJob(dictJob): return SurfaceJob(dictToSurf(dictJob['surf']),dictToCalc(dictJob['calc']),dictJob['dftCode'])

def safeMakeDir(x): 
	try: os.makedirs(str(x))
	except OSError: pass

def jobKind(strJob): 
	return 'surface' if '___' in strJob else 'bulk' #surface jobs have '___' in their string representation
"""
