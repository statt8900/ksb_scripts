import os,sys,time
from db         import getCluster,insertObject,updateDB
from bulkjob    import db2BulkJob
from bulkresult import BulkResult
from surfacejob import db2SurfaceJob,SurfaceResult

from fireworks.core.firework          import FiretaskBase,FWAction
from fireworks.utilities.fw_utilities import explicit_serialize

from ase.parallel   import paropen,parprint
from ase.eos        import EquationOfState
from ase.units      import kJ
from ase            import io
from copy 			import deepcopy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from scipy.optimize import fmin,leastsq
import numpy as np

"""
Defines Firetasks to be called in a submit method of a Job object

If having problems, maybe get rid of setupEnv function?
"""

########################
# FIRETASKS
######################
@explicit_serialize
class OptimizeLattice(FiretaskBase):
	def run_task(self,fw_spec):
		global nIter
		nIter = 0
		
		jobID = fw_spec['jobID']
		job   = db2BulkJob(jobID)

		setupEnv(job)

		job.printFields() # confirm job details

		existsPrecalc = job.precalcxc != 'None'
		t0 = time.time()

		def getBulkEnergy(latticeParams):
			"""For a given set of bravais lattice parameters, optimize atomic coordinates and return minimum energy"""
			global nIter
			nIter += 1
			atoms = job.fromParams(latticeParams)
			calc  = job.calc(restart=True) if existsPrecalc else job.calc()

			if existsPrecalc: 
				preCalc = job.calc(precalc=True)
				job.optimizePos(atoms,preCalc,saveWF = True)
				if job.dftcode =='gpaw': 
					atoms,calc = job.gpawRestart()

			job.optimizePos(atoms,calc)

			energy = atoms.get_potential_energy()
			
			with paropen('lattice_opt.log','a') as logfile: 
				logfile.write('%s\t%s\n' %(energy,latticeParams))
				parprint('%s\t%s\n' %(energy,latticeParams))

			if job.dftcode=='gpaw': 
				atoms.calc.write('inp.gpw', mode='all') #for use in getXCContribs

			io.write('out.traj',atoms)
			return energy

		optimizedLatticeParams = fmin(getBulkEnergy,job.iGuess(),ftol=1,xtol=job.xtol)
		avgtime                = (time.time()-t0)/float(nIter)/60.0 #min

		### Get info from final state
		optAtoms = io.read('out.traj')
		optCell  = optAtoms.get_cell()
		pos      = optAtoms.get_scaled_positions()
		e0       = optAtoms.get_potential_energy()
		f0       = optAtoms.get_forces()
		magmom   = optAtoms.get_magnetic_moments() if any([x>0 for x in job.magmoms]) else np.zeros(len(optAtoms))
		resultDict = {'latticeParams': optimizedLatticeParams,'avgtime':avgtime,'optCell':optCell,'pos':pos,'e0':e0,'f0':f0,'magmom':magmom,'niter':nIter}

		return FWAction( stored_data=resultDict,mod_spec=[{'_push': resultDict}])

@explicit_serialize
class GetBulkModulus(FiretaskBase):
	def run_task(self,fw_spec): 

		jobID,latticeParams = fw_spec['jobID'],fw_spec['latticeParams'] # WHY IS LATTICEPARAMS [[FLOAT]] INSTEAD OF [FLOAT]?
		job                 = db2BulkJob(jobID)

		setupEnv(job)

		existsPrecalc = job.precalcxc != 'None'
		optAtoms  = job.fromParams(latticeParams[0]) 
		optCell   = optAtoms.get_cell()
		strains   = np.linspace(1-job.strain,1+job.strain,9)
		calc      = job.calc(restart=True) if existsPrecalc else job.calc()

		vol,eng = [],[]

		for i, strain in enumerate(strains):
			atoms = deepcopy(optAtoms)
			atoms.set_cell(optCell*strain,scale_atoms=True)

			if existsPrecalc: 
				preCalc = job.calc(precalc=True)
				job.optimizePos(atoms,preCalc,saveWF = True)
				if job.dftcode =='gpaw': 
					atoms,calc = job.gpawRestart()

			job.optimizePos(atoms,calc)

			energy = atoms.get_potential_energy()
			volume = atoms.get_volume()
			vol.append(deepcopy(volume))
			eng.append(deepcopy(energy))
			parprint('%f %f'%(strain,energy))

		aHat,quadR2 = quadFit(np.array(deepcopy(vol)),np.array(deepcopy(eng)))

		try:		
			eos = EquationOfState(vol,eng)
			v0, e0, b = eos.fit()
			eos.plot(filename='bulk-eos.png',show=False)
			b0= b/kJ*1e24 #GPa use this value if EOS doesn't fail
		except ValueError:  # too bad of a fit for ASE to handle
			b0 = aHat*2*vol[4]*160.2 # units: eV/A^6 * A^3 * 1, where 1 === 160.2 GPa*A^3/eV


		return FWAction( stored_data=        {'b0':b0,'quadR2':quadR2}
						,mod_spec=[{'_push': {'b0':b0,'quadR2':quadR2}}])
@explicit_serialize
class GetXCcontribs(FiretaskBase):
	def run_task(self,fw_spec): 

		jobID = fw_spec['jobID']
		job   = db2BulkJob(jobID)
		setupEnv(job)
	
		if job.dftcode == 'gpaw':
			from gpaw        import restart
			from gpaw.xc.bee import BEEFEnsemble
			atoms,calc = restart('inp.gpw', setups='sg15', xc='mBEEF', convergence={'energy': 5e-4}, txt='mbeef.txt')
			atoms.get_potential_energy()
			beef = BEEFEnsemble(calc)
			xcContribs = beef.mbeef_exchange_energy_contribs()
		else: 
			xcContribs = np.zeros((8,8))

		return FWAction( stored_data={'xcContribs': xcContribs}
						,mod_spec=[{'_push': {'xcContribs': xcContribs}}])


@explicit_serialize
class SaveResults(FiretaskBase):
	def run_task(self,fw_spec): 
		jobID	= fw_spec['jobID'] 
		job 	= db2BulkJob(jobID)

		resultNames = ['optCell','pos','magmom','e0','f0','b0','quadR2','xcContribs','avgtime','niter']
		results     = [fw_spec[x][0] for x in resultNames]
		resultDict  = {n:fw_spec[n] for n in resultNames}
	
		result = BulkResult(*([jobID]+results))
		
		result.addToASEdb(job)

		insertObject(result) # add result to bulkresult table
		
		updateDB('bulkjob','status','id',jobID,'complete')

		return FWAction(stored_data= resultDict)
#################################################################
#################################################################
#################################################################
@explicit_serialize
class RelaxSurface(FiretaskBase):
	def run_task(self,fw_spec): 
		jobID	= fw_spec['jobID'] 
		job 	= db2SurfaceJob(jobID)
		
		t0 				= time.time()
		# no precalc???
		atoms = job.atoms
		job.optimizePos(atoms,job.calc())
		
		t 		= (time.time()-t0)/60.0 #min
		niter 	= 1 # ??? how to calculate this (from log file? different for gpaw and qe)
 		avgtime = t/niter

		pos 	= atoms.get_scaled_positions()
		e0 		= atoms.get_potential_energy()
		f0      = atoms.get_forces()
		magmom  = atoms.get_magnetic_moments() if any([x>0 for x in job.magmoms]) else np.zeros(len(atoms))
		
		resultDict={'pos':pos,'magmom':magmom,'e0':e0,'f0':f0,'avgtime':avgtime,'niter':niter}
		
		if job.dftcode=='gpaw': 
			atoms.calc.write('inp.gpw', mode='all') #for use in getXCContribs
		
		io.write('out.traj',atoms)
		
		return FWAction(stored_data= resultDict,mod_spec=[{'_push': resultDict}])

@explicit_serialize
class GetSurfaceXCcontribs(FiretaskBase):
	def run_task(self,fw_spec): 

		jobID = fw_spec['jobID']
		job   = db2SurfaceJob(jobID)
		setupEnv(job)
	
		if job.dftcode == 'gpaw':
			from gpaw        import restart
			from gpaw.xc.bee import BEEFEnsemble
			atoms,calc = restart('inp.gpw', setups='sg15', xc='mBEEF', convergence={'energy': 5e-4}, txt='mbeef.txt')
			atoms.get_potential_energy()
			beef = BEEFEnsemble(calc)
			xcContribs = beef.mbeef_exchange_energy_contribs()
		else: 
			xcContribs = np.zeros((8,8))

		return FWAction( stored_data={'xcContribs': xcContribs}
						,mod_spec=[{'_push': {'xcContribs': xcContribs}}])


@explicit_serialize
class SaveSurfResults(FiretaskBase):
	def run_task(self,fw_spec): 
		jobID	= fw_spec['jobID'] 
		job 	= db2SurfaceJob(jobID)
		
		resultNames = ['pos','magmom','e0','f0','xcContribs','avgtime','niter']
		results     = [fw_spec[x][0] for x in resultNames]
		resultDict  = {n:fw_spec[n] for n in resultNames}
		
		result = SurfaceResult(*([jobID]+results))
		
		result.addToASEdb(job)
		
		insertObject(result) # add result to surface result table
		
		updateDB('surfacejob','status','id',jobID,'complete')
		
		return FWAction(stored_data= resultDict)

#############################################
#############################################

def quadFit(xIn,yIn):
	ax = plt.gca()
	
	# Center data around 4th data point
	x = xIn-xIn[4]; y = yIn-yIn[4]
	def model(a):  return a*np.square(x)
	def errVec(a): return a*np.square(x) - y   #create fitting function of form mx+b
	aHat, success = leastsq(errVec, [0.1])
	yhat 	= model(aHat)
	ybar 	= np.sum(y)/len(y)          # or sum(y)/len(y)
	ssresid = np.sum((yhat-y)**2)  
	sstotal = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
	r2 		= 1-(ssresid / float(sstotal))
	ax.plot(x,y,'ro')
	ax.set_xlabel('Volume, A^3')
	ax.set_ylabel('Energy, eV')
	ax.plot(x,yhat)
	ax.annotate('r2 = %f'%(r2),xy=(0,0),fontsize=12)
	plt.savefig('quadFit.png')
	plt.figure()
	return float(aHat[0]),float(r2)

def setupEnv(j):
	"""
	Not sure if this does anything, maybe remove?
	"""
	d = j.dftcode
	c = getCluster()

	if d == 'gpaw': 
		sys.path.append('/scratch/users/ksb/gpaw/gpaw_sg15/lib/python2.7/site-packages')
		os.environ["PATH"]+='/scratch/users/ksb/gpaw/gpaw_sg15/bin'
	elif d == 'quantumespresso': 
		pass
	else: 
		raise NotImplementedError, 'dftcode %s has no system environment in standardScripts.setupEnv()'%j.dftcode



