import os,json
import numpy as np

from objectClass import Object
from result      import Result

from cluster     import clusterDict,getCluster
from dft         import dftDict
from lattice     import structDict
from printParse  import underscoreSep,derscoreSep,multiLine
from calc        import Calc,parseCalc,pspDict
from bulk        import parseBulk,dictToBulk
from surface     import parseSurface,SurfacePure
from miniDatabase import magElems

"""
This module defines the Job object
A Job object should contain sufficient information to uniquely determine a DFT optimization
SurfaceJobs and BulkJobs are subclasses of Job.
A variety of auxillary functions are found at the bottom of this module.
Many of these interface with JSON
"""

class Job(Object):
	def jobPth(self,clusterName): return clusterDict[clusterName].root+'/'+str(self)

	def initialize(self,writeCluster,runCluster):
		"""
		Creates directory (if it doesn't exist), writes files such that only metarun.sh must be entered to submit the job
		"""
		#try: 
		path    = self.jobPth(writeCluster)
		mPath   = path+'/'+runCluster+'_metarun.sh'
		rPath   = path+'/'+runCluster+'_run.sh'
		oPath   = path+'/'+runCluster+'_opt.py'
		if os.path.exists(mPath): os.remove(mPath)
		if os.path.exists(rPath): os.remove(rPath);
		if os.path.exists(oPath): os.remove(oPath);
		cluster = clusterDict[runCluster]
		safeMakeDir(path)
		try: self.surf.makeAdsorptionPlot(path)
		except AttributeError: pass
		self.getStruct().guessTraj(path) # save initial structure to traj file
		with open (mPath,'w') as f: f.write(self.metarunFile(runCluster))
		with open (oPath,'w') as f: f.write(self.genericScript(cluster))
		with open (rPath,'w') as f: f.write(self.runFile(cluster))
		os.chmod(mPath,0o775); os.chmod(rPath,0o775)

	def metarunFile(self,runCluster):
		"""
		Contents of file that is run (for any cluster, any dftcode) in order to submit batch job
		"""
		return multiLine(	['#!/bin/bash'
							,'if [ ! -f result.json ]; then'	
							,dftDict[self.dftCode].metarunLine(runCluster)	# Don't repeat already completed jobs
							,'fi\n'])

	def runFile(self,cluster): return cluster.runFileHeader(self) #Intermediate file, calls opt.py

	def makeCalc(self,preCalc=False): return dftDict[self.dftCode].calcFunc(self,preCalc) 

	def getBulkResult(self): 
		"""
		Searches for result in current cluster, returning a Result object if found. (Warning printed if not found)
		"""
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



class SurfaceJob(Job):
	def __init__(self,surf,calc,dftCode,timeLim=None,preCalc=None):
		self.surf    = surf       # Surface object
		self.calc    = calc       # Calculator object
		self.preCalc = preCalc    # Calculator object
		self.dftCode = dftCode    # String
		self.timeLim = surfTimeFunc(surf,calc,dftCode) if timeLim is None else timeLim # Float, hours

	def __str__(self):  return underscoreSep(3,[str(self.surf),str(self.calc),self.dftCode])

	def nAtoms(self): 	return len(self.surf.makeAtoms()) # number of atoms in final structure

	def getStruct(self): return self.surf

	def genericScript(self,cluster):
		"""
		Contents of opt.py, a simple relaxation of a slab 
		"""
		return multiLine(	[cluster.python
							,cluster.header(self.timeLim)+'\n'
							,'import auto,job,time,os,copy'
							,'\n# Initialize'
							,'jobObject = job.parseSurfaceJob("{0}")'.format(self)
							,'start=time.time()'
							,'\n#Energy Minimizing Function'
							,'atoms=auto.surfRelax(copy.deepcopy(jobObject)) '
							,'cell = atoms.get_cell()'
							,'pos    = atoms.get_positions() # NOT scaled'  
							,'magmom = atoms.get_magnetic_moments() if jobObject.calc.magmom() else jobObject.zeros()'
							,'e0     = atoms.get_potential_energy()'
							,'b      = None'
							,'totalTime = (time.time()-start)/3600'
							,'\n#Save intermediate results'
							,'auto.saveResults(jobObject,cell,pos.tolist(),magmom.tolist(),e0,None,[[]],totalTime)'
							,'\n#If time remaining, calculate XC contributions'
							,'#xcContribs = auto.getXCcontribs(jobObject) #reads inp.pw generated in surfRelax'
							,'\n#Save results'
							,'#auto.saveResults(jobObject,cell,pos.tolist(),magmom.tolist(),e0,None,xcContribs.tolist(),totalTime)'])


class BulkJob(Job):
	def __init__(self,bulk,calc,dftCode,timeLim=None,preCalc=Calc('PBE',500,(5,5,5),5,1,5,500,12,1,1,None,None)):
		"""Creates a bulk optimization job."""
		self.bulk    = bulk       # Bulk
		self.calc    = calc       # Calculator object
		self.preCalc = preCalc    # Calculator object
		self.dftCode = dftCode    # String
		self.timeLim = bulkTimeFunc(bulk,calc,dftCode) if timeLim is None else timeLim # Float, hours
		setattr(self.preCalc, 'magmom', self.calc.magmom) # make precalc use same pseudopotential
		setattr(self.preCalc, 'psp', self.calc.psp)       # make precalc use same pseudopotential

	def __str__(self):  return underscoreSep(2,[str(self.bulk),str(self.calc),self.dftCode])

	def nAtoms(self): return len(self.bulk.makeAtoms())

	def getStruct(self):   return self.bulk
	def magmom(self):      return self.calc.magmom > 0
	def checkResult(self): return self.getBulkResult() is not None
	def genericScript(self,cluster):
		"""
		Contents of opt.py, which optmizes lattice parameters (and atomic coordinates if necessary)
		"""
		return multiLine(	[cluster.python
							,cluster.header(self.timeLim)+'\n'
							,'import time,os'
							,'from auto import getEnergy,getBulkModulus,getXCcontribs,saveResults'
							,'from job import parseBulkJob'
							,'from scipy.optimize import fmin'
							,'\n# Initialize'
							,'jobObject = parseBulkJob("{0}")'.format(self)
							,'start=time.time()'
							,'\n#Energy Minimizing Function'
							,'latticeParams=fmin(getEnergy,{0},ftol=2,xtol=jobObject.calc.xtol*0.001,args=({0},))'
							,'\n#Bulk Modulus Function'
							,'cell,nIter,pos,magmom,e0,b,bFit = getBulkModulus(jobObject,latticeParams)'
							,'\n#XC contribution Function'
							,'xcContribs = getXCcontribs(jobObject) #(reads inp.pw generated in BulkMod)'
							,'\n#Save results'
							,'totalTime = (time.time()-start)/3600'
							,'saveResults(jobObject,cell.tolist(),pos.tolist(),magmom.tolist(),e0,b,bFit,xcContribs.tolist(),totalTime,nIter)'])


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
