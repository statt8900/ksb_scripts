
#External Modules
import os,json,traceback,base64,pickle,fireworks,ase,matplotlib,ast
import scipy.optimize 	as opt
import ase.optimize  	as aseopt
import ase.eos 			as aseeos
import ase.thermochemistry as asethermo
import copy 	as c		
import glob 	as g			
import numpy	as np

matplotlib.use('Agg')
from matplotlib import pyplot as plt

# Internal Modules
import jobs, misc

#############
# FIRETASKS
#############

@fireworks.utilities.fw_utilities.explicit_serialize
class OptimizeLattice(fireworks.core.firework.FiretaskBase):
	def run_task(self,fw_spec):
		global nIter
		job,params,initatoms = initialize(fw_spec)
		nIter,energies,lparams = 0,[],[]
		print 'initialized'
		if os.path.exists('lastguess.log'):
			with open('lastguess.log','r') as f: iGuess = pickle.loads(f.read())
			print 'read lastguess ',iGuess
		else: 
			iGuess = job.iGuess()
			print 'using initial guess ',iGuess
		
		for d in ['qn.traj','qn.log']:
			if os.path.exists(d): 
				os.remove(d)
				print 'removed existing file ',d
		
		def getBulkEnergy(latticeParams):
			"""For a given set of bravais lattice parameters, optimize atomic coordinates and return minimum energy"""
			global nIter
			nIter += 1
			with open('lastguess.log','w') as f: f.write(pickle.dumps(latticeParams))
			atoms = job.fromParams(latticeParams)
			optimizePos(atoms,job.calc(),job.fmax())
			energy = atoms.get_potential_energy()
			energies.append(energy);lparams.append(latticeParams)
			ase.io.write('out.traj',atoms)
			return energy
		
		optimizedLatticeParams 	= opt.fmin(getBulkEnergy,iGuess,ftol=1,xtol=job.xtol())
		print 'optimized lattice params',optimizedLatticeParams
		optAtoms 				= ase.io.read('out.traj')
		
		resultDict 	= misc.mergeDicts([params,trajDetails(optAtoms)
									,{'raw_energy': 		optAtoms.get_potential_energy()
										,'forces_pckl':		pickle.dumps(optAtoms.get_forces())
										,'niter': 			nIter
										,'latticeopt_pckl':	pickle.dumps(zip(energies,lparams))}])	
		
		with open('result.json', 'w') as outfile: json.dumps(resultDict, outfile)
		return fireworks.core.firework.FWAction(stored_data=resultDict)


@fireworks.utilities.fw_utilities.explicit_serialize
class GetBulkModulus(fireworks.core.firework.FiretaskBase):
	def run_task(self,fw_spec): 
		job,params,optAtoms = initialize(fw_spec)
		optCell   		= optAtoms.get_cell()
		strains   		= np.linspace(1-job.strain(),1+job.strain(),9)		
		vol,eng 		= [],[]
		
		for i, strain in enumerate(strains):
			atoms = optAtoms.copy()
			atoms.set_cell(optCell*strain,scale_atoms=True)
			optimizePos(atoms,job.calc(),job.fmax())
			volume,energy = atoms.get_volume(),atoms.get_potential_energy()
			vol.append(c.deepcopy(volume));eng.append(c.deepcopy(energy))
		aHat,quadR2 = quadFit(np.array(c.deepcopy(vol)),np.array(c.deepcopy(eng)))

		try:		
			eos = aseeos.EquationOfState(vol,eng)
			v0, e0, b = eos.fit()
			eos.plot(filename='bulk-eos.png',show=False)
			b0= b/ase.units.kJ*1e24 									#GPa: use this value if EOS doesn't fail
			with open('bulk-eos.png', 'rb') as f: img = base64.b64encode(f.read())

		except ValueError:  								# too bad of a fit for ASE to handle
			b0 = aHat*2*vol[4]*160.2						# units: eV/A^6 * A^3 * 1, where 1 === 160.2 GPa*A^3/eV
			img = None


		resultDict 	= misc.mergeDicts([params,trajDetails(optAtoms)
									,{'bulkmod':b0,'bfit':quadR2,'bulkmodimg_base64':img,'voleng_json':json.dumps(zip(vol,eng))}])

		with open('result.json', 'w') as outfile:   json.dumps(resultDict, outfile)
		return fireworks.core.firework.FWAction( stored_data=resultDict)

@fireworks.utilities.fw_utilities.explicit_serialize
class GetXCcontribs(fireworks.core.firework.FiretaskBase):
	def run_task(self,fw_spec): 
		from gpaw        import restart
		from gpaw.xc.bee import BEEFEnsemble
		job,params,atoms = initialize(fw_spec)
		atoms.set_calculator(job.PBEcalc())
		atoms.get_potential_energy()
		atoms.calc.write('inp.gpw', mode='all') 

		atoms,calc = restart('inp.gpw', setups='sg15', xc='mBEEF', convergence={'energy': 5e-4}, txt='mbeef.txt')
		atoms.get_potential_energy()
		beef 		= BEEFEnsemble(calc)
		xcContribs 	= beef.mbeef_exchange_energy_contribs()
		ase.io.write('final.traj',atoms)
		resultDict 	= misc.mergeDicts([params,trajDetails(ase.io.read('final.traj'))
									,{'xcContribs': pickle.dumps(xcContribs)}])

		with open('result.json', 'w') as outfile:   json.dumps(resultDict, outfile)
		return fireworks.core.firework.FWAction( stored_data=resultDict)

#################################################################
@fireworks.utilities.fw_utilities.explicit_serialize
class Relax(fireworks.core.firework.FiretaskBase):
	def run_task(self,fw_spec): 
		job,params,initatoms = initialize(fw_spec)
		print 'initialized'
		if os.path.exists('qn.traj'): 
			if os.stat('qn.traj').st_size > 100: 	
				initatoms = ase.io.read('qn.traj')
				print 'reading from qn.traj'
			else: 									
				os.remove('qn.traj')
				print 'removing empty qn.traj'
		if os.path.exists('qn.log') and os.stat('qn.log').st_size < 100: 
			os.remove('qn.log')
			print 'removing empty qn.log'

		optimizePos(initatoms,job.calc(),job.fmax())
		print 'optimized positions'
		ase.io.write('final.traj',initatoms)
		
		niter = getRelaxIter(job.dftcode())
		print 'got nIter'
		resultDict 	= misc.mergeDicts([params,trajDetails(ase.io.read('final.traj'))
									,{'raw_energy':initatoms.get_potential_energy()
									,'forces_pckl':pickle.dumps(initatoms.get_forces())
									,'niter':niter} ])
		print 'made result dict'
		with open('result.json', 'w') as outfile:   json.dumps(resultDict, outfile)
		return fireworks.core.firework.FWAction(stored_data=resultDict)

@fireworks.utilities.fw_utilities.explicit_serialize
class Vibrations(fireworks.core.firework.FiretaskBase):
	def run_task(self,fw_spec): 
		"""ONLY WORKS FOR QUANTUM ESPRESSO"""
		from ase.vibrations       import Vibrations
		job,params,atoms = initialize(fw_spec)
		
		prev = g.glob('*.pckl') #delete incomplete pckls - facilitating restarted jobs
		for p in prev:
			if os.stat(p).st_size < 100: os.remove(p)
		
		atoms.set_calculator(job.vibCalc())
		vib = Vibrations(atoms,delta=job.delta(),indices=job.vibids())


		vib.run(); vib.write_jmol()
		vib.summary(log='vibrations.txt')
		with open('vibrations.txt','r') as f: vibsummary = f.read()

		vib_energies,vib_frequencies = vib.get_energies(),vib.get_frequencies()
		resultDict 	= misc.mergeDicts([params,	{'vibfreqs_pckl': pickle.dumps(vib_frequencies)
												,'vibsummary':vibsummary
												,'vibengs_pckl':pickle.dumps(vib_energies)}])

		with open('result.json', 'w') as outfile: json.dumps(resultDict, outfile)
		return fireworks.core.firework.FWAction(stored_data=resultDict)

@fireworks.utilities.fw_utilities.explicit_serialize
class DOS(fireworks.core.firework.FiretaskBase):
	def run_task(self,fw_spec):
		job,params,atoms = initialize(fw_spec)
		dos = None
		raise NotImplementedError
		resultDict 	= misc.mergeDicts([params,	{}])
		with open('result.json', 'w') as outfile:   json.dumps(resultDict, outfile)
		return fireworks.core.firework.FWAction(stored_data=resultDict)

@fireworks.utilities.fw_utilities.explicit_serialize
class NEB(fireworks.core.firework.FiretaskBase):
	def run_task(self,fw_spec):
		job,params,atoms = initialize(fw_spec)
		neb = None
		raise NotImplementedError
		resultDict 	= misc.mergeDicts([params,	{}])
		with open('result.json', 'w') as outfile:   json.dumps(resultDict, outfile)
		return fireworks.core.firework.FWAction(stored_data=resultDict)

@fireworks.utilities.fw_utilities.explicit_serialize
class VCRelax(fireworks.core.firework.FiretaskBase):
	def run_task(self,fw_spec):

		job,params,atoms = initialize(fw_spec)
		
		if not os.path.exists('intermediate.traj'):
			atoms.set_calculator(job.vccalc())
			energy = atoms.get_potential_energy() #trigger espresso to be launched
			ase.io.write('intermediate.traj',atoms.calc.get_final_structure())

		atoms = ase.io.read('intermediate.traj')
		atoms.set_calculator(job.vccalc())
		energy = atoms.get_potential_energy() #trigger espresso to be launched
		ase.io.write('final.traj',atoms.calc.get_final_structure())
		e0,f0 	= atoms.get_potential_energy(),atoms.get_forces()
		atoms 	= ase.io.read('final.traj')

		resultDict 	= misc.mergeDicts([params,trajDetails(atoms),
									{'raw_energy': 	e0
									,'forces_pckl': pickle.dumps(f0)}])

		with open('result.json', 'w') as outfile:   json.dumps(resultDict, outfile)

		return fireworks.core.firework.FWAction(stored_data=resultDict)


#################################################################
#################################################################
#################################################################
		

def quadFit(xIn,yIn):
	""" Input x vector: units A^3, Input y vector: units eV """
	ax = plt.gca()
	# Center data around 4th data point
	x = xIn-xIn[4]; y = yIn-yIn[4]
	def model(a):  return a*np.square(x)
	def errVec(a): return a*np.square(x) - y   #create fitting function of form mx+b
	aHat, success = opt.leastsq(errVec, [0.1])
	yhat 	= model(aHat)
	ybar 	= np.sum(y)/len(y)          # or sum(y)/len(y)
	ssresid = np.sum((yhat-y)**2)  
	sstotal = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
	r2 		= 1 - (ssresid / float(sstotal)) # 
	ax.plot(x,y,'ro');	ax.plot(x,yhat)
	ax.set_xlabel('Volume, A^3'); ax.set_ylabel('Energy, eV')
	ax.annotate('r2 = %f'%(r2),xy=(0,0),fontsize=12)
	plt.savefig('quadFit.png')
	plt.figure()
	return float(aHat[0]),float(r2)

def optimizePos(atoms,calc,fmax):
	atoms.set_calculator(calc)
	dyn = aseopt.BFGS(atoms=atoms, logfile='qn.log', trajectory='qn.traj',restart='qn.pckl')
	dyn.run(fmax=fmax)

def clearOldFiles():
	def keepNewest(listx):	
		for x in sorted(listx)[:-1]: os.remove(x)
	try:
		keepNewest(g.glob('*.error'))
		keepNewest(g.glob('*.out'))
	except OSError: pass

def initialize(fwspec):
	clearOldFiles()
	prms 	= fwspec['params']
	jb 		= jobs.assignJob(prms)
	atms 	= jb.atoms().copy()
	ase.io.write('init.traj',atms)
	return jb,prms,atms

def trajDetails(atoms):
	try: mag = atoms.get_magnetic_moments()
	except: mag = np.array([0]*len(atoms))
	return {'finaltraj_pckl':pickle.dumps(atoms)
			,'finalpos_pckl':pickle.dumps(atoms.get_positions())
			,'finalcell_pckl':pickle.dumps(atoms.get_cell())
			,'finalmagmom_pckl':pickle.dumps( mag )}

def getRelaxIter(dftcode):
	if dftcode == 'gpaw':	
		calcs = g.glob('*.txt') #gpaw calculator OUGHT be only .txt file?
		assert len(calcs) == 1
		with open(calcs[0],'r') as f: out = f.read()
		return out.count('Converged')

	elif dftcode == 'quantumespresso':
		with open('calcdir/log','r') as f: out =  f.read()
		return out.count('convergence')
	else: 
		raise NotImplementedError

