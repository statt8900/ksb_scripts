import os,sys,time

from db         import getCluster,insertObject,updateDB
from job    	import Job,db2object
from result 	import Result

from fireworks.core.firework          import FiretaskBase,FWAction
from fireworks.utilities.fw_utilities import explicit_serialize

from ase.parallel   import paropen,parprint
from ase.db 		import connect
from ase.eos        import EquationOfState
from ase.units      import kJ
from ase            import io,Atoms
from copy 			import deepcopy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from scipy.optimize import fmin,leastsq
import numpy as np

#############
# FIRETASKS
#############

@explicit_serialize
class OptimizeLattice(FiretaskBase):
	def run_task(self,fw_spec):
		global nIter
		nIter = 0
		
		jobID = fw_spec['jobID']
		job   = db2object(jobID)
		
		existsPrecalc = job.precalc != 'None'
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
		magmom   = optAtoms.get_magnetic_moments() if any([x>0 for x in job.magmomsinit()]) else np.zeros(len(optAtoms))
		resultDict = {'latticeParams': optimizedLatticeParams,'avgtime':avgtime,'cell':optCell,'pos':pos,'e0':e0,'f0':f0,'magmom':magmom,'niter':nIter}
		
		
		return FWAction( stored_data=resultDict,mod_spec=[{'_push': resultDict}])
		
@explicit_serialize
class GetBulkModulus(FiretaskBase):
	def run_task(self,fw_spec): 
		
		jobID,latticeParams = fw_spec['jobID'],fw_spec['latticeParams'] # WHY IS LATTICEPARAMS [[FLOAT]] INSTEAD OF [FLOAT]?
		job   = db2object(jobID)

		existsPrecalc 	= job.precalc != 'None'
		optAtoms  		= job.fromParams(latticeParams[0]) 
		optCell   		= optAtoms.get_cell()
		strains   		= np.linspace(1-job.strain,1+job.strain,9)
		calc      		= job.calc(restart=True) if existsPrecalc else job.calc()
		
		vol,eng = [],[]
		
		for i, strain in enumerate(strains):
			atoms = optAtoms.copy()
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
			with paropen('bulkmod.log','a') as logfile: 
				logfile.write('%s\t%s\n' %(strain,energy))
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
		job   = db2object(jobID)
	
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
		job   	= db2object(jobID)

		if 'relax' in job.jobkind:
			pos  	= fw_spec['pos'][0]
			cell 	= fw_spec['cell'][0] if job.jobkind == 'bulkrelax' else job.cellinit()
			magmom 	= fw_spec['magmom'][0]

			db    = connect('/scratch/users/ksb/db/ase.db')	

			info  = {'name': 		job.name()
					,'emt': 		0 				# EMT for relaxed structures useless, only relevant for deciding when to relax something
					,'relaxed': 	True 			# Could always doing this be a problem?
					,'comments':	'Generated from #%d'%(job.aseidinitial)
					,'parent': 		job.aseidinitial
					,'jobid': 		job.jobid 		# Int

					###Unchanged Things###
					,'kind': 		job.kind()} 		# String
				
			if job.structure() 	is not None: info['structure'] = 	job.structure()
			if job.vacancies() 	is not None: info['vacancies'] = 	job.vacancies()
			if job.sites() 		is not None: info['sites'] = 		job.sites()
			if job.facet() 		is not None: info['facet'] = 		job.facet()
			if job.xy() 		is not None: info['xy'] = 			job.xy()
			if job.layers()		is not None: info['layers'] = 		job.layers()
			if job.constrained() is not None: info['constrained'] = job.constrained()
			if job.symmetric() 	is not None: info['symmetric'] = 	job.symmetric()
			if job.vacuum() 	is not None: info['vacuum'] = 		job.vacuum()
			if job.adsorbates() is not None: info['adsorbates'] = 	job.adsorbates()

			optAtoms = Atoms(symbols 		 	= job.symbols()
							,scaled_positions 	= pos
							,cell 				= cell
							,magmoms 			= magmom
							,constraint 		= job.constraints()
							,tags 				= job.tags())

			aseid = db.write(optAtoms,key_value_pairs=info)
		else: 
			aseid = None

		launchdir 	= fw_spec['_launch_dir']
		energy 		= fw_spec['e0'][0] 		if 'relax' in job.jobkind else None
		forces 		= fw_spec['f0'][0] 		if 'relax' in job.jobkind else None

		bulkmodulus = fw_spec['b0'][0] 		if job.jobkind == 'bulkrelax' else None 
		bfit 		= fw_spec['quadR2'][0] 	if job.jobkind == 'bulkrelax' else None

		xccoeffs 	= fw_spec['xcContribs'][0] if 'relax' in job.jobkind else None

		viblist 	= fw_spec['vibs'][0]	if job.jobkind == 'vib' else None
		dos 		= fw_spec['dos'][0] 	if job.jobkind == 'dos' else None	
		barrier 	= fw_spec['barrier'][0] if job.jobkind == 'neb' else None

		time 		= fw_spec['avgtime'][0] 
		niter 		= fw_spec['niter'][0]  	if 'relax' in job.jobkind else None

		result = Result(jobID,launchdir,aseid=aseid,energy=energy,forces=forces
				,bulkmodulus=bulkmodulus,bfit=bfit,xccoeffs=xccoeffs,viblist=viblist
				,dos=dos,barrier=barrier,time=time,niter=niter)
		

		resultDict  = {'_launch_dir':launchdir,'aseid':aseid,'energy':energy
						,'forces':str(forces),'bulkmodulus':bulkmodulus
						,'bfit':bfit,'xccoeffs':xccoeffs,'time':time
						,'niter':niter,'viblist':viblist,'barrier':barrier}

		insertObject(result) # add result to bulkresult table
		
		updateDB('status','jobid',jobID,'complete')

		return FWAction(stored_data= resultDict)
#################################################################
#################################################################
#################################################################
@explicit_serialize
class Relax(FiretaskBase):
	def run_task(self,fw_spec): 
		jobID	= fw_spec['jobID'] 
		job 	= db2object(jobID)
		
		t0 				= time.time()

		atoms = job.atoms()
		job.optimizePos(atoms,job.calc())
		
		t 		= (time.time()-t0)/60.0 #min
		niter 	= 1 # ??? how to calculate this (from log file? different for gpaw and qe)
 		avgtime = t/niter
		
		pos 	= atoms.get_scaled_positions()
		e0 		= atoms.get_potential_energy()
		f0      = atoms.get_forces()
		magmom  = atoms.get_magnetic_moments() if any([x>0 for x in job.magmomsinit()]) else np.zeros(len(atoms))
		
		resultDict={'pos':pos,'magmom':magmom,'e0':e0,'f0':f0,'avgtime':avgtime,'niter':niter}
		
		if job.dftcode=='gpaw': 
			atoms.calc.write('inp.gpw', mode='all') #for use in getXCContribs
		
		io.write('out.traj',atoms)
		
		return FWAction(stored_data= resultDict,mod_spec=[{'_push': resultDict}])

@explicit_serialize
class Vibrations(FiretaskBase):
	def run_task(self,fw_spec): 
		"""ONLY WORKS FOR QUANTUM ESPRESSO"""

		from ase.vibrations       import Vibrations

		jobID	= fw_spec['jobID']
		job 	= db2object(jobID)
		t0 		= time.time()
		atoms 	= job.atoms()
		
		atoms.set_calculator(job.vibCalc())

		vib = Vibrations(atoms,delta=job.delta(),indicies=job.vibids)

		vib.run()
		vib.summary(log='vibrations.txt')
		vib.write_jmol()

		vibs 	= vib.get_frequencies()
		t 		= (time.time()-t0)/60.0 #min
	
		return FWAction(stored_data= {'vibs':vibs,'avgtime':t},mod_spec=[{'_push': {'vibs':vibs,'time':t}}])

class DOS(FiretaskBase):
	def run_task(self,fw_spec):
		jobID	= fw_spec['jobID'] 
		job 	= db2object(jobID)
		t0 		= time.time()
		atoms 	= job.atoms()
		
		t 		= (time.time()-t0)/60.0 #min

		dos = None

		return FWAction(stored_data= {'dos':dos,'avgtime':t},mod_spec=[{'_push': {'vibs':vibs,'time':t}}])


class NEB(FiretaskBase):
	def run_task(self,fw_spec):
		return NotImplementedError
#################################################################
#################################################################
#################################################################
		

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

