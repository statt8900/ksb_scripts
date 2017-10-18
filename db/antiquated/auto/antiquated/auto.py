import sys, warnings
warnings.filterwarnings("ignore", message="Moved to ase.neighborlist")

import numpy as np 
from copy           import deepcopy
from scipy.optimize import leastsq

from ase            import optimize,io
from ase.parallel   import paropen,parprint,parallel_function
from ase.eos        import EquationOfState
from ase.units      import kJ

from bulkjob        import db2BulkJob
from calc           import db2Calc
from bulk           import db2Bulk
from result         import BulkResult

from db      import updateStatus,insertObject
"""
This module contains the main functions that are called by optimization scripts
For bulk relaxations, the function getEnergy is minimized by fmin.
Given an optimized cell, getBulkModulus returns a lot of data about it.
surfRelax is called for surface jobs
"""
####################################################################################
####################################################################################
# FUNCTIONS FOR SCRIPTS
####################################################################################
####################################################################################
def initialize(jobID,latticeParams):
	jobObject        = db2BulkJob(jobID)
	bulkObject       = db2Bulk(jobObject.bulkpth)
	calcObject       = db2Calc(jobObject.calcid) 
	dft              = jobObject.dftCode
	preCalcObject    = deepcopy(calcObject)
	preCalcObject.xc = jobObject.precalcxc
	existsPrecalc    = str(preCalcObject.xc) != 'None'
	atoms            = bulkObject.fromParams(latticeParams)
	magmoms          = bulkObject.magInit(calcObject.magmom)
	return (jobObject,bulkObject,calcObject,dft,preCalcObject,preCalcObject,existsPrecalc,atoms,magmoms)

def getBulkEnergy(latticeParams,jobID):
	"""
	For a given set of bravais lattice parameters, optimize atomic coordinates and return minimum energy
	"""
	### INITIALIZE
	jobObject,bulkObject,calcObject,dft,preCalcObject,preCalcObject,existsPrecalc,atoms,magmoms = initialize(jobID,latticeParams)
	### PRECALCULATE
	if existsPrecalc: optimizePos(atoms,dft,preCalcObject,magmoms,restart=False,saveWF = True)
	### CALCULATE	
	optimizePos(atoms,dft,calcObject,magmoms,restart=existsPrecalc,saveWF=False)
	### GET ENERGY
	energy = atoms.get_potential_energy()
	with paropen('lattice_opt.log','a') as logfile: logfile.write('%s\t%s\n' %(energy,latticeParams))
	parprint('%s\t%s\n' %(energy,latticeParams))
	return energy

def optimizePos(atoms,dft,calcObj,magInit,restart,saveWF):
	### CREATE CALC
	calcMaker   = gpawCalc if dft == 'gpaw' else qeCalc
	calc        = calcMaker(calcObj,restart=restart)
	### INITIALIZE ATOMS OBJECT
	atoms.set_calculator(calc)
	atoms.set_initial_magnetic_moments(magInit)

	### OPTIMIZE POSITIONS IF NECESSARY
	maxForce = np.amax(abs(atoms.get_forces()))
	if maxForce > calcObj.fMax:
		parprint("max force = %f, optimizing positions"%(maxForce))
		dyn = optimize.BFGS(atoms=atoms, logfile='qn.log', trajectory='qn.traj',restart='qn.pckl')
		dyn.run(fmax=calcObj.fMax)

	if saveWF and dft=='gpaw': 
		atoms.calc.write('preCalc_inp.gpw', mode='all') #for use in getXCContribs
		atoms,calc = gpawRestart(calcObj)

################################################
################################################
# Bulk Modulus
################################################
################################################
def getBulkModulus(jobID,latticeParams): 
	"""
	Use optimized lattice parameters from fmin(getEnergy). Calculate energies around the minimum. 
	Returns optimized atomic positions, energy of minimum, calculated bulk modulus, and writes wavefunctions to file inp.gpw
	"""			
	### INITIALIZE
	vol,eng = [],[]
	jobObject,bulkObject,calcObject,dft,preCalcObject,preCalcObject,existsPrecalc,optAtoms,magmoms = initialize(jobID,latticeParams)
	with open('lattice_opt.log','r') as f: 
		nIter = len(f.readlines())+9 #count number of ionic steps taken
	optCell = optAtoms.get_cell()
	strains = np.linspace(0.99,1.01,9)

	### GENERATE dE/dV plot
	for i, strain in enumerate(strains):
		atoms = deepcopy(optAtoms)
		atoms.set_cell(optCell*strain,scale_atoms=True)
		### PRECALCULATE
		if existsPrecalc: 
			optimizePos(atoms,dft,preCalcObject,magmoms,restart=False,saveWF = True)
		### CALCULATE	
		optimizePos(atoms,dft,calcObject,magmoms,restart=existsPrecalc,saveWF=False)
		### COLLECT RESULTS
		energy = atoms.get_potential_energy()
		volume = atoms.get_volume()
		vol.append(deepcopy(volume))
		eng.append(deepcopy(energy))
		parprint('%f %f'%(strain,energy))
		### COLLECT DETAILED DATA ABOUT MINIMUM ENERGY STATE
		if i==4:
			pos = atoms.get_scaled_positions()
			vMin  = deepcopy(volume)
			io.write('optimized.traj',atoms)
			magmom = atoms.get_magnetic_moments() if calcObject.magmom > 0 else np.zeros(len(atoms))
			if dft=='gpaw': 
				atoms.calc.write('inp.gpw', mode='all') #for use in getXCContribs

	aHat,quadR2 = quadFit(np.array(deepcopy(vol)),np.array(deepcopy(eng)))
	b0 = aHat*2*vMin*160.2 # units: eV/A^6 * A^3 * 1, where 1 === 160.2 GPa*A^3/eV
	parprint('quadR2 = '+str(quadR2))
	try:		
		eos = EquationOfState(vol,eng)
		v0, e0, b = eos.fit()
		eos.plot(filename='bulk-eos.png',show=False)
		b0= b/kJ*1e24 #GPa use this value if EOS doesn't fail
	except ValueError: e0=eng[4]
	return (optCell,nIter,pos,magmom,e0,b0,quadR2) #GPa


def quadFit(xIn,yIn):
	import matplotlib
	matplotlib.use('Agg')
	from matplotlib import pyplot as plt
	ax = plt.gca()

	# Center data around 4th data point
	x = xIn-xIn[4]; y = yIn-yIn[4]
	def model(a):  return a*np.square(x)
	def errVec(a): return a*np.square(x)-y   #create fitting function of form mx+b
	aHat, success = leastsq(errVec, [0.1])
	yhat     = model(aHat)
	ybar     = np.sum(y)/len(y)          # or sum(y)/len(y)
	ssresid  = np.sum((yhat-y)**2)  
	sstotal  = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
	r2       = 1-(ssresid / float(sstotal))
	ax.plot(x,y,'ro')
	ax.set_xlabel('Volume, A^3')
	ax.set_ylabel('Energy, eV')
	ax.plot(x,yhat)
	ax.annotate('r2 = %f'%(r2),xy=(0,0),fontsize=12)
	plt.savefig('quadFit.png')
	plt.figure()
	return float(aHat[0]),float(r2)

################################################
# Other Bulk Functions
################################################
def getXCcontribs(jobID): 
	jobObject = db2BulkJob(jobID)
	if jobObject.dftCode == 'gpaw':
		from gpaw        import restart
		from gpaw.xc.bee import BEEFEnsemble
		atoms,calc = restart('inp.gpw', setups='sg15', xc='mBEEF', convergence={'energy': 5e-4}, txt='mbeef.txt')
		atoms.get_potential_energy()
		beef=BEEFEnsemble(calc)
		return  beef.mbeef_exchange_energy_contribs()
	else: return np.zeros((8,8))

def saveResults(jobID,cell,pos,magmom,e0,b,bFit,xc,t,nIter):	
	result = BulkResult(jobID,cell,pos,magmom,e0,b,bFit,xc,t,nIter) #lattice constants determined downstream, depending on structure'
	if rank() == 0:
		insertObject(result)
		with open("result.json","w") as f: print >> f, result.toJSON()
		updateStatus(jobID,'running','completed')
################################################
################################################
# Calculator Functions
################################################
################################################
def gpawCalc(calc,restart=False,surf=False):
	"""
	Accepts a job (either Surface or Bulk) and returns a Calculator
	"""
	from gpaw import GPAW,PW,Davidson,Mixer,MixerSum,FermiDirac
		
	return GPAW(mode         = PW(calc.pw)                         #eV
				,xc          = calc.xc
				,kpts        = {'density': calc.kpt, 'gamma': True} if isinstance(calc.kpt,int) else calc.kpt
				,spinpol     = calc.magmom>0
				,convergence = {'energy':calc.eConv} #eV/electron
				,mixer       = ((MixerSum(beta=calc.mixing,nmaxold=calc.nmix,weight=100)) 
								if calc.magmom>0 else (Mixer(beta=calc.mixing,nmaxold=calc.nmix,weight=100)))
				,maxiter       = calc.maxsteps
				,nbands        =  -1*calc.nbands
				,occupations   = FermiDirac(calc.sigma)
				,setups        = 'sg15' 
				,eigensolver   = Davidson(5)
				,poissonsolver = None # {'dipolelayer': 'xy'} if isinstance(self,SurfaceJob) else None
				,txt=str(calc)+'.txt'
				,symmetry={'do_not_symmetrize_the_density': True} #ERROR IN LI bcc 111 surface
				)

def gpawRestart(calc):
	from gpaw import restart,PW,Davidson,Mixer,MixerSum,FermiDirac
	return restart('preCalc_inp.gpw',mode    = PW(calc.pw)
									,xc      = calc.xc
									,kpts    =  {'density': calc.kpt, 'gamma': True} if isinstance(calc.kpt,int) else calc.kpt
									,spinpol = calc.magmom>0
									,convergence = {'energy':calc.eConv} #eV/electron
									,mixer =                  ((MixerSum(beta=calc.mixing,nmaxold=calc.nmix,weight=100)) 
											if calc.magmom>0 else (Mixer(beta=calc.mixing,nmaxold=calc.nmix,weight=100)))
									,maxiter       = calc.maxsteps
									,nbands        =  -1*calc.nbands
									,occupations   = FermiDirac(calc.sigma)
									,setups        = 'sg15' #(pspDict[calc.psp]).pthFunc[cluster]
									,eigensolver   = Davidson(5)
									,poissonsolver = None # {'dipolelayer': 'xy'} if isinstance(self,SurfaceJob) else None
									,txt=str(calc)+'.txt'
									,symmetry={'do_not_symmetrize_the_density': True} #ERROR IN LI bcc 111 surface
									)

def qeCalc(calc,restart=False,surf=False):
	from espresso import espresso	
	return espresso( pw      = calc.pw
					,dw      = calc.pw*10
					,xc      = calc.xc
					,kpts    =  {'density': calc.kpt, 'gamma': True} if isinstance(calc.kpt,int) else calc.kpt
					,spinpol = calc.magmom > 0
					,convergence = {'energy':calc.eConv
									,'mixing':calc.mixing
									,'nmix':calc.nmix
									,'maxsteps':calc.maxsteps
									,'diag':'david'}
					,nbands = -1*calc.nbands
					,sigma  = calc.sigma
					,dipole = {'status':True if surf else False}
					,outdir = 'calcdir'	 # output directory
					,startingwfc='file' if restart else 'atomic+random'
					,psppath= 0 
					,output = {'removesave':True})

def rank():
	# Check for special MPI-enabled Python interpreters:
	if '_gpaw' in sys.builtin_module_names:
	    import _gpaw 	    # http://wiki.fysik.dtu.dk/gpaw
	    world = _gpaw.Communicator()
	elif '_asap' in sys.builtin_module_names:
	    import _asap # http://wiki.fysik.dtu.dk/asap, can't import asap3.mpi here (import deadlock)
	    world = _asap.Communicator()
	elif 'asapparallel3' in sys.modules: # Older version of Asap
	    import asapparallel3
	    world = asapparallel3.Communicator()
	elif 'Scientific_mpi' in sys.modules:
	    from Scientific.MPI import world
	elif 'mpi4py' in sys.modules:
	    world = MPI4PY()
	else:
	    world = DummyMPI()# This is a standard Python interpreter:
	rank = world.rank
	size = world.size
	return rank

################################################
################################################
# Surface functions
################################################
################################################
def surfRelax(jobID): raise NotImplementedError
