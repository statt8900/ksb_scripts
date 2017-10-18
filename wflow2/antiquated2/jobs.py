#External Modules
import warnings
warnings.filterwarnings("ignore", message="Moved to ase.neighborlist")
warnings.filterwarnings("ignore", message="matplotlib.pyplot")
import os,time,ast,random,json,pickle,re,math
import numpy 	as np
import math 	as m
import datetime as d		
# Internal Modules
import misc,dbase,cluster,emt
import printParse as pp


"""
Description:


Notes:
"Espresso executable doesn't seem to have been started"
	- caused by including /home/vossj/suncat/bin in PATH
"""

################################################
# Job base - must be inherited by all jobs
##############################################

class JobBase(object):
	"""Methods valid for all jobs"""

	def __init__(self,params):
		self.params = params
	
	### All jobs [should] have these attributes
	def jobkind(self): 		return self.params['jobkind']	 		# STRING
	def inittraj_pckl(self): return self.params['inittraj_pckl'] 	# PICKLE
	def comments(self): 	return self.params['comments']			# STRING
	def trajcomments(self): return self.params['trajcomments'] 		# STRING
	def name(self): 		return self.params['name'] 				# STRING
	def relaxed(self): 		return self.params['relaxed']			# INT # 0 means not relaxed, otherwise is FWID referring to relax job
	def kind(self): 		return self.params['kind'] 				# STRING (molecule/bulk/surface)

	def pw(self): 			return self.params['pw'] 				# INT
	def xc(self): 			return self.params['xc'] 				# STRING
	def kptden(self): 		return self.params['kptden'] 			# INT
	def psp(self): 			return self.params['psp'] 				# STRING
	def dwrat(self): 		return self.params['dwrat'] 			# FLOAT
	def econv(self): 		return self.params['econv'] 			# FLOAT
	def mixing(self): 		return self.params['mixing'] 			# FLOAT
	def nmix(self): 		return self.params['nmix'] 				# INT
	def maxstep(self):		return self.params['maxstep'] 			# INT
	def nbands(self): 		return self.params['nbands'] 			# INT
	def sigma(self): 		return self.params['sigma'] 			# FLOAT
	def fmax(self): 		return self.params['fmax'] 				# FLOAT
	def dftcode(self): 		return self.params['dftcode'] 			# STRING

	def generalCheck(self):
		essential = set(['jobkind','inittraj_pckl','comments','trajcomments','name','relaxed','kind'
						,'pw','xc','kptden','psp','dwrat','econv','mixing','nmix','maxstep','nbands','sigma','fmax','dftcode']) & self.kindEssentials() & self.jobEssentials()
		actual 	= set(self.params.keys())
		assert essential <= actual, 'Missing the following general keys: '+str(essential - actual)

	# Derived static info (known at the time of job creation), potentially add more 
	def natoms(self): 			return len(self.atoms())								# INT
	def symbols_pckl(self): 	return pickle.dumps(self.symbols())						# NP.ARRAY
	def symbols_str(self): 		return str(self.symbols())								# STRING
	def comp_pckl(self): 		return pickle.dumps(self.comp())						# NP.ARRAY
	def comp_str(self): 		return str(self.comp())								# STRING
	def species_pckl(self): 	return pickle.dumps(self.species())						# NP.ARRAY
	def species_str(self): 		return str(self.species())								# STRING
	def metalstoich_pckl(self): return pickle.dumps(self.metalstoich())						# NP.ARRAY
	def metalstoich_str(self): 	return str(self.metalstoich())								# STRING
	def metalspecies_pckl(self):return pickle.dumps(self.metalspecies())						# NP.ARRAY
	def metalspecies_str(self): return str(self.metalspecies())								# STRING
	def metalcomp_pckl(self): 	return pickle.dumps(self.metalcomp())						# NP.ARRAY
	def metalcomp_str(self): 	return str(self.metalcomp())								# STRING
	def numbers_pckl(self): 	return pickle.dumps(self.numbers())						# NP.ARRAY
	def numbers_str(self): 		return pp.printNp(self.numbers())						# STRING
	def cellinit_pckl(self): 	return pickle.dumps(self.cellinit()) 					# NP.MATRIX
	def cellinit_str(self): 	return pp.printNp(self.cellinit()) 						# STRING
	def posinit_pckl(self):	 	return pickle.dumps(self.posinit())						# NP.MATRIX
	def posinit_str(self):	 	return pp.printNp(self.posinit()) 						# STRING
	def magmomsinit_pckl(self):	return pickle.dumps(self.magmomsinit()) 				# NP.ARRAY
	def magmomsinit_str(self):	return pp.printNp(self.magmomsinit()) 					# STRING
	def tags_pckl(self): 		return pickle.dumps(self.tags())						# NP.ARRAY
	def tags_str(self): 		return pp.printNp(self.tags()) 							# STRING
	def constraints_pckl(self): return pickle.dumps(self.constraints())					# NP.ARRAY
	def constraints_str(self): 	return str(self.constraints()) 							# STRING
	def paramsinit_pckl(self): 	return pickle.dumps(self.paramsinit()) 					# NP.ARRAY ANGLES IN RADIANS
	def paramsinit_str(self): 	return str(self.paramsinit()) 					# STRING
	def spinpol(self):			return int(any([x>0 for x in self.magmomsinit()]))  	# INT (BOOL)
	def dipole(self): 			return int(True if self.kind() == 'surface' else False)	# INT (BOOL)
	def kpt_pckl(self): 		return pickle.dumps(self.kpt()) 						# NP.ARRAY
	def kpt_str(self): 			return pp.printNp(self.kpt()) 							# STRING

	def emt(self):				
		a = self.atoms().copy()
		a.set_calculator(emt.EMT())
		return a.get_potential_energy()/len(a) 				# eV/Atom

	def emtsym(self):				
		a = self.atoms().copy()
		a.set_calculator(emt.EMT())
		return int(np.amax(abs(a.get_forces())) < 0.0001) 	# Int (Bool)

	def calcLabel(self): return "%s-%s-%s-%s"%(self.xc(),self.pw(),self.kptden(),self.dftcode())

	def generalSummary(self): 
		"""summarize some of the invariant methods of jobs in a dictionary"""
		general = 	{'params_json':		json.dumps(self.params)
					,'natoms': 			self.natoms()
					,'symbols_pckl':	self.symbols_pckl()
					,'symbols_str':		self.symbols_str()
					,'comp_pckl':		self.comp_pckl()
					,'comp_str':		self.comp_str()
					,'species_pckl':	self.species_pckl()
					,'species_str':		self.species_str()
					,'numbers_pckl':	self.numbers_pckl()
					,'numbers_str':		self.numbers_str()
					,'metalstoich_pckl':self.metalstoich_pckl()
					,'metalstoich_str':	self.metalstoich_str()
					,'metalspecies_pckl':self.metalspecies_pckl()
					,'metalspecies_str':self.metalspecies_str()
					,'metalcomp_pckl':	self.metalcomp_pckl()
					,'metalcomp_str':	self.metalcomp_str()
					#,'cellinit_pckl':	self.cellinit_pckl()
					#,'cellinit_str':	self.cellinit_str()
					#,'posinit_pckl':	self.posinit_pckl()
					#,'posinit_str':		self.posinit_str()
					#,'magmomsinitpckl':	self.magmomsinit_pckl()
					#,'magmomsinitstr':	self.magmomsinit_str()
					#,'tags_pckl':		self.tags_pckl()
					#,'tags_str':		self.tags_str()
					#,'constraints_pckl':self.constraints_pckl()
					#,'constraints_str':	self.constraints_str()
					,'paramsinit_pckl':	self.paramsinit_pckl()
					#,'paramsinit_str':	self.paramsinit_str()
					,'kpt_pckl':		self.kpt_pckl()
					,'kpt_str':			self.kpt_str()
					,'emt':				self.emt()
					,'emtsym': 			self.emtsym()
					,'spinpol':			self.spinpol()
					,'dipole':			self.dipole()
					,'calclabel': 		self.calcLabel()
					,'strjob':	 		str(self)
				} 
		return misc.mergeDicts([self.params, general, self.jobSummary(), self.kindSummary()])

	def __str__(self): return '_'.join(	[ self.numbers_str() 
										+ self.posinit_str() 
										+ self.cellinit_str() 
										+ self.magmomsinit_str() 
										+ self.jobkind() 
										+ self.strcalc() 
										+ self.strkind() 
										+ self.strjob()])

	def strcalc(self): return '_'.join([str(x) for x in [self.params[y] for y in ['pw','xc','kptden','psp','dwrat','econv','mixing','nmix','maxstep','nbands','sigma','fmax','dftcode']]])
	##############################################################################
	# Dynamic info/methods, could change, or otherwise not stored in FW database
	##############################################################################
	def psppath(self):
		pspDict = 	{'sherlock':
						{'gbrv15pbe':'/home/vossj/suncat/psp/gbrv1.5pbe'}
					,'suncat':
						{'gbrv15pbe':'/nfs/slac/g/suncatfs/sw/external/esp-psp/gbrv1.5pbe'}}
		return pspDict[misc.getCluster()][self.psp()]

	def atoms(self): 		return pickle.loads(self.inittraj_pckl()) 				# BINARY
	def numbers(self): 		return self.atoms().get_atomic_numbers()				# NP.ARRAY
	def symbols(self): 		return self.atoms().get_chemical_symbols()				# [STRING]
	def species(self): 		return list(set(self.symbols()))						# [STRING] - set of elements in traj
	def comp(self): 		return misc.normalizeList(self.symbols())				# [STRING] - elements in traj but divided by GCD
	def metalstoich(self): 	return [x for x in self.symbols() if x not in misc.nonmetalSymbs] # [STRING]
	def metalspecies(self): return list(set(self.metalstoich())) 					# [STRING]
	def metalcomp(self): 	return misc.normalizeList(self.metalstoich())			# [STRING] - metal elements in traj but divided by GCD
	def posinit(self):	 	return self.atoms().get_scaled_positions()				# NP.MATRIX
	def cellinit(self): 	return self.atoms().get_cell()							# NP.MATRIX
	def magmomsinit(self):	return self.atoms().get_initial_magnetic_moments() 		# NP.ARRAY
	def paramsinit(self): 	return misc.cell2param(self.cellinit()) 				# NP.ARRAY ANGLES IN RADIANS
	def constraints(self):	return self.atoms().constraints 						# NP.ARRAY
	def tags(self): 		return self.atoms().get_tags()							# NP.ARRAY
		
	def kpt(self):  																# NP.ARRAY
		def kptdensity2monkhorstpack(atoms, kptdensity, even=True):
		    """Convert k-point density to Monkhorst-Pack grid size. """
		    recipcell,kpts = atoms.get_reciprocal_cell(),[]
		    for i in range(3):
		        k = 2 * 3.14159 * m.sqrt((recipcell[i]**2).sum()) * kptdensity
		        if even: 	kpts.append(2 * int(np.ceil(k / 2)))
		        else: 		kpts.append(int(np.ceil(k)))
		    return kpts

		kpt = kptdensity2monkhorstpack(self.atoms(),self.kptden())
		if   self.kind()=='surface': 	return np.array(kpt[:2]+[1])
		elif self.kind()=='molecule': 	return np.array([1,1,1])
		else: 							return np.array(kpt)

	def new(self): 
		return not dbase.anyQuery(pp.STRJOB(str(self)))#re.sub('[\s+]', '', str(self))))

	def submit(self): 
		from fireworks  import LaunchPad
		if self.new():
			launchpad 	= LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')
			wflow 		= self.wflow()
			launchpad.add_wf(wflow)
			time.sleep(2) # folder names are unique due to timestamp to the nearest second
			return 1
		else: 
			print 'repeat!'
			return 0

####################################################################
# Calculator selection: choose GPAWCalcJob, or QECalcJob
####################################################################

class GPAWCalcJob(object):
	def calc(self):
		from gpaw import GPAW,PW,Davidson,Mixer,MixerSum,FermiDirac
		return GPAW(mode         = PW(self.pw())                         #eV
					,xc          = self.xc()
					,kpts        = self.kpt()
					,spinpol     = self.spinpol()
					,convergence = {'energy':self.econv()} #eV/electron
					,mixer       = ((MixerSum(beta=self.mixing(),nmaxold=self.nmix(),weight=100)) 
									if self.spinpol() else (Mixer(beta=self.mixing(),nmaxold=self.nmix(),weight=100)))
					,maxiter       = self.maxstep()
					,nbands        = self.nbands()
					,occupations   = FermiDirac(self.sigma())
					,setups        = self.psp()
					,eigensolver   = Davidson(5)
					,poissonsolver = None # {'dipolelayer': 'xy'} if isinstance(self,SurfaceJob) else None
					,txt='%d_%s.txt'%(self.pw(),self.xc())
					,symmetry={'do_not_symmetrize_the_density': True}) #ERROR IN LI bcc 111 surface

	def PBEcalc(self):
		from gpaw import GPAW,PW,Davidson,Mixer,MixerSum,FermiDirac
		return GPAW(mode         = PW(self.pw())                         #eV
					,xc          = 'PBE'
					,kpts        = self.kpt()
					,spinpol     = self.spinpol()
					,convergence = {'energy':self.econv()} #eV/electron
					,mixer       = ((MixerSum(beta=self.mixing(),nmaxold=self.nmix(),weight=100)) 
									if self.spinpol() else (Mixer(beta=self.mixing(),nmaxold=self.nmix(),weight=100)))
					,maxiter       = self.maxstep()
					,nbands        = self.nbands()
					,occupations   = FermiDirac(self.sigma())
					,setups        = self.psp()
					,eigensolver   = Davidson(5)
					,poissonsolver = None # {'dipolelayer': 'xy'} if isinstance(self,SurfaceJob) else None
					,txt='%d_%s.txt'%(self.pw(),self.xc())
					,symmetry={'do_not_symmetrize_the_density': True}) #ERROR IN LI bcc 111 surface



	def allocate(self,n):
		return [cluster.sherlockgpaw]*n  

class QECalcJob(object):		
	def calc(self):
		from espresso import espresso	

		return espresso( pw      	= self.pw()
						,dw      	= self.pw()*self.dwrat()
						,xc      	= self.xc()
						,kpts    	= self.kpt()
						,spinpol 	= self.spinpol()
						,convergence= 	{'energy':	self.econv()
										,'mixing':	self.mixing()
										,'nmix': 	self.nmix()
										,'maxsteps':self.maxstep()
										,'diag':	'david'}
						,nbands 	= self.nbands()
						,sigma  	= self.sigma()
						,dipole 	= {'status': self.dipole()}
						,outdir 	= 'calcdir'	 # output directory
						,startingwfc= 'atomic+random' #alternatively, 'file' for restart
						,psppath 	= self.psppath()
						,output 	= {'removesave':True})

	def allocate(self,n):
		rand=[]
		def odds(randomx):
			if randomx < 0.6: 	return cluster.suncat
			elif randomx < 1: 		return cluster.suncat2
		for i in range(n): rand.append(random.random())
		return [odds(x) for x in rand]

############################################
# Methods depending on kind of ase.db object: select either BulkJob, SurfaceJob, or MoleculeJob
############################################
class BulkJob(object): 
	##########
	# Required
	##########
	def structure(self): 			return self.params['structure'] 	# String
	def bulkvacancy_json(self): 	return self.params['bulkvacancy_json'] # 
	def bulkscale_json(self): 		return self.params['bulkscale_json']   #

	#########
	# Derived
	##########
	def bulkvacancy_str(self): 		return str(json.loads(self.bulkvacancy_json())) # 
	def bulkscale_str(self): 		return str(json.loads(self.bulkscale_json()))   #
	########################################
	# Necessary methods for any traj 'kind'
	########################################

	def kindEssentials(self): return set(['structure','bulkvacancy_json','bulkscale_json'])
	def kindSummary(self): 	return {'bulkvacancy_str':self.bulkvacancy_str()
									,'bulkscale_str':self.bulkscale_str()}
	def strkind(self): 			return '_'.join([self.structure()+self.bulkvacancy_str(),self.bulkscale_str()])

class SurfaceJob(object):
	##################
	# Necessary params
	##################
	# From Bulk
	def bulkparent(self): 		return self.params['bulkparent'] 					# String
	def structure(self): 		return self.params['structure'] 					# String
	def bulkvacancy_json(self): return self.params['bulkvacancy_json'] 				# JSON [Int]
	def bulkscale_json(self): 	return self.params['bulkscale_json']   				# JSON [Int[
	# New for surface
	def sites_base64(self): 	return self.params['sites_base64'] 					# STRING (base64)
	def facet_json(self): 		return self.params['facet_json']					# JSON [Int]
	def xy_json(self): 			return self.params['xy_json']						# JSON [Int]
	def layers(self): 			return self.params['layers']						# Int
	def constrained(self): 		return self.params['constrained']					# Int
	def symmetric(self): 		return self.params['symmetric']						# Int
	def vacuum(self): 			return self.params['vacuum'] 						# Int
	def adsorbates_json(self): 	return self.params['adsorbates_json']	 			# JSON Map {String:[String]}

	#############################
	# Derived - stored in summary
	#############################
	def bulkvacancy_str(self): 		return str(json.loads(self.bulkvacancy_json())) 	# STRING
	def bulkscale_str(self): 		return str(json.loads(self.bulkscale_json()))   	# STRING
	def facet_str(self): 			return str(self.facet())							# 
	def xy_str(self): 				return str(self.xy()) 								# 
	def surfaceArea(self): 			return np.linalg.norm(np.cross(self.cellinit()[0],self.cellinit()[1])) * (2 if self.symmetric() else 1)
	def adsorbates_str(self): 		return pp.printAds(self.adsorbates())

	#############################
	# Derived, not stored
	#############################

	def adsorbates(self): 	return json.loads(self.adsorbates_json())	# MAP STRING [STRING]
	def facet(self):		return json.loads(self.facet_json()) 		# [INT]
	def xy(self): 			return json.loads(self.xy_json()) 			# [INT]

	# Necessary methods for any kind of traj
	def strkind(self): 		return '_'.join([self.structure()+self.bulkvacancy_str(),self.bulkscale_str()
									,self.facet_str(),self.xy_str(),str(self.layers()),str(self.constrained())
									,str(self.symmetric()),str(self.vacuum()),self.adsorbates_str()])

	def kindSummary(self): 	return {'bulkvacancy_str':self.bulkvacancy_str()
									,'bulkscale_str':self.bulkscale_str()
									,'facet_str':self.facet_str()
									,'xy_str':self.xy_str()
									,'surfaceArea': self.surfaceArea()
									,'adsorbates_str':self.adsorbates_str()}

	def kindEssentials(self): return  set(['bulkparent','structure','bulkvacancy_json','bulkscale_json','sites_base64'
											,'facet_json','xy_json','layers','constrained','symmetric'
											,'vacuum','adsorbates_json']) 

class MoleculeJob(object): 
	def strkind(self): 			return ''
	def kindSummary(self): 		return {}
	def kindEssentials(self): 	return set([])



##################################################################################################
# Methods depending on type of relaxation: select from latopt,bulkmod,xc,relax,dos,neb,vib,vcrelax
##################################################################################################

class LatOptMethods(object):
	#################
	# Required Params
	#################

	def xtol(self): return self.params['xtol']		
	######################
	# Required job methods
	#######################
	def strjob(self): return str(self.xtol())
	def jobEssentials(self): return set([])
	def jobSummary(self): return {'dof':self.dof()}

	#################
	# Other methods
	#################

	def fromParams(self,params): 
		""" 
		Params is a list of between 1 and 6 numbers (a,b,c,alpha,beta,gamma). 
		ANGLES OF INPUT ARE IN RADIANS
		ASE CELLS WANT ANGLES IN DEGREES
		Depending on structure, we can construct the cell from a subset of these parameters
		"""
		#UNIT CONVERSION
		
		atoms = self.atoms().copy() # returns scaled version of initial structure
		a = params[0]
		if self.structure() in ['cubic','cscl']:     
			cell = [a,a,a,90,90,90]
		elif self.structure() in ['fcc','diamond','zincblende','rocksalt']:  
			cell = [a,a,a,60,60,60]
		elif self.structure() in ['bcc']:
			alpha = 109.47122
			cell = [a,a,a,alpha,alpha,alpha]
		elif self.structure() in ['hcp','hexagonal']: 
			cell = [a,a,params[1],90,90,120]
		elif self.structure() =='triclinic': 
			cell = [params[0],params[1],params[2],m.degrees(params[3]),m.degrees(params[4]),m.degrees(params[5])]
		else: raise NotImplementedError, 'job.fromParams() cannot handle unknown structure = '+self.structure()

		atoms.set_cell(cell,scale_atoms=True)
		return atoms

	def iGuess(self):
		a,b,c,alpha,beta,gamma = self.paramsinit() #ANGLES IN RADIANS
		if self.structure() in ['fcc','bcc','rocksalt','diamond','cscl','zincblende']: return [a]
		elif self.structure() in ['hexagonal']: return [a,c]
		elif self.structure() in ['triclinic']: return self.paramsinit()
		else: raise ValueError, 'Bad entry in "structure" field for Atoms object info dictionary: '+self.structure()
	def dof(self): return len(self.iGuess())

	def guessTime(self): return min(40, 3 * len(self.iGuess()) * self.pw()/500.0 * self.kptden()/2*(4 if self.xc()=='mBEEF' else 1))
	def jobCheck(self): assert 'xtol' in self.params.keys()
	#def repeat(self):	return AND(self.calcEQ()+self.generalEQ()+self.kindEQ()+[XTOL(self.xtol())])
	def wflow(self): 
		from fireworks   		import Firework,Workflow
		from standardScripts 	import OptimizeLattice
		cluster = self.allocate(1)[0]
		timestamp = '_'+d.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
		name = 'OptimizeLattice'
		fw1 = Firework([OptimizeLattice()]
						,name = name
						,spec = {'params': 			self.params
								,'_fworker': 		cluster.fworker
								,'_add_launchpad_and_fw_id': 	True
								,"_queueadapter": 	cluster.qfunc(self.guessTime())
								,"_launch_dir": 	cluster.launchdir+name+timestamp})
		return Workflow([fw1],name=name)

class VCRelaxMethods(object):
	######################
	# Required job methods
	#######################

	def jobEssentials(self): return set([])
	def jobSummary(self): return {}

	def vccalc(self): 
		from espresso import espresso
		return  espresso( pw            	= self.pw()
				                ,dw           	= self.pw()*self.dwrat()
				                ,xc            	= self.xc()
				                ,kpts          	= self.kpt()
				                ,nbands        	= self.nbands()
								,dipole 	= {'status': self.dipole()}
				                ,sigma         	= self.sigma()
				                ,mode          	= 'vc-relax'
				                ,cell_dynamics	= 'bfgs'
				                ,opt_algorithm 	= 'bfgs'
				                ,cell_factor   	= 2.
				                ,spinpol       	= self.spinpol()
								,outdir 		= 'calcdir'	 
								,output 		= {'removesave':True}
				                ,psppath 	  	= self.psppath()
								,convergence= 	{'energy':	self.econv()
												,'mixing':	self.mixing()
												,'nmix': 	self.nmix()
												,'maxsteps':self.maxstep()
												,'diag':	'david'})


	def guessTime(self): return min(40,  3 * self.pw()/500.0 * self.kptden()/2.0*(4 if self.xc()=='mBEEF' else 1))
	def jobCheck(self): pass
	def strjob(self): return ''
	def wflow(self): 
		from fireworks   		import Firework,Workflow
		from standardScripts 	import VCRelax
		cluster = self.allocate(1)[0]
		timestamp = '_'+d.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
		name = 'VCRelax'
		fw1 = Firework([VCRelax()]
						,name = name
						,spec = {'params': 			self.params
								,'_fworker': 		cluster.fworker
								,'_add_launchpad_and_fw_id': 	True
								,"_queueadapter": 	cluster.qfunc(self.guessTime())
								,"_launch_dir": 	cluster.launchdir+name+timestamp})
		return Workflow([fw1],name=name)

class BulkModMethods(object):
	#################
	# Required Params
	#################
	def strain(self): return self.params['strain']
	######################
	# Required job methods
	#######################

	def jobEssentials(self): return set(['strain'])
	def jobSummary(self): return {}
	def strjob(self): return str(self.strain())

	def guessTime(self): return min(40,  3 * self.pw()/500.0*(4 if self.xc()=='mBEEF' else 1))

	def wflow(self): 
		from fireworks   		import Firework,Workflow
		from standardScripts 	import GetBulkModulus
		cluster = self.allocate(1)[0]
		timestamp = '_'+d.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
		name = 'GetBulkModulus'
		fw1 = Firework([GetBulkModulus()]
						,name = name
						,spec = {'params': self.params
								,'_fworker': 		cluster.fworker
								,'_add_launchpad_and_fw_id': 			True
								,"_queueadapter": 	cluster.qfunc(self.guessTime())
								,"_launch_dir": 	cluster.launchdir+name+timestamp})
		return Workflow([fw1],name=name)

class XCMethods(object):
	######################
	# Required job methods
	#######################
	def strjob(self): return ''
	def jobEssentials(self): return set([])
	def jobSummary(self): return {}

	def guessTime(self): return 10

	def wflow(self):
		from fireworks   		import Firework,Workflow
		from standardScripts 	import GetXCcontribs
		cluster = self.allocate(1)[0]
		timestamp = '_'+d.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
		name = 'GetXCcontribs'
		fw1 = Firework([GetXCcontribs()]
						,name = name
						,spec = {'params': self.params
								,'_fworker': 		cluster.fworker
								,'_add_launchpad_and_fw_id': 			True
								,"_queueadapter": 	cluster.qfunc(self.guessTime())
								,"_launch_dir": 	cluster.launchdir+name+timestamp})
		return Workflow([fw1],name=name)

class RelaxJobMethods(object):
	######################
	# Required job methods
	#######################
	def strjob(self): return ''

	def jobEssentials(self): return set([])
	def jobSummary(self): return {}

	def guessTime(self): return min(40,  5*self.pw()/500.0 * self.kptden()*(4 if self.xc()=='mBEEF' else 1))

	def repeat(self):	return AND(self.calcEQ()+self.generalEQ()+self.kindEQ())
	def wflow(self): 
		from fireworks   		import Firework,Workflow
		from standardScripts 	import Relax
		cluster = self.allocate(1)[0]
		timestamp = '_'+d.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
		name = 'Relaxation'
		fw1 = Firework([Relax()]
						,name = name
						,spec = {'params': 			self.params
								,'_fworker': 		cluster.fworker
								,'_add_launchpad_and_fw_id': 			True
								,"_queueadapter": 	cluster.qfunc(self.guessTime())
								,"_launch_dir": 	cluster.launchdir+name+timestamp})

		return Workflow([fw1],name=name)


class VibJobMethods(object):
	#################
	# Required Params
	#################

	def vibids_json(self): 	return self.params['vibids_json'] 	# JSON [Int]
	def delta(self): 		return self.params['delta'] 	# Float

	# Derived - stored in memory
	def vibids_str(self): return str(self.vibids())
	
	######################
	# Required job methods
	#######################
	def strjob(self): return str(self.vibids())+'_'+str(self.delta())
	def jobEssentials(self): return set(['vibids_json','delta'])
	def jobSummary(self): return {'vibids_str':self.vibids_str()}

	def guessTime(self): 
		print self.vibids(),self.pw(),self.kptden()
		return 5 * len(self.vibids()) * self.pw()/500.0 * self.kptden()

	def wflow(self): 
		from fireworks   		import Firework,Workflow
		from standardScripts 	import Vibrations
		cluster = self.allocate(1)[0]
		timestamp = '_'+d.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
		name = 'Vibration'
		fw1 = Firework([Vibrations()]
						,name = name
						,spec = {'params': 			self.params
								,'_fworker': 		cluster.fworker
								,'_add_launchpad_and_fw_id': 			True
								,"_pass_job_info": 	True
								,"_queueadapter": 	cluster.qfunc(self.guessTime())
								,"_launch_dir": 	cluster.launchdir+name+timestamp})

		return Workflow([fw1],name=name)

	###############
	# Other methods
	###############
	def vibids(self): return json.loads(self.vibids_json())

	def vibCalc(self):
		from espresso.vibespresso import vibespresso
		return vibespresso( pw      	= self.pw()
							,dw      	= self.pw()*self.dwrat()
							,xc      	= self.xc()
							,kpts    	= self.kpt()
							,spinpol 	= self.spinpol()
							,convergence= 	{'energy':	self.econv()
											,'mixing':	self.mixing()
											,'nmix': 	self.nmix()
											,'maxsteps':self.maxstep()
											,'diag':	'david'}
							,nbands 	= self.nbands()
							,sigma  	= self.sigma()
							,dipole 	= {'status': self.dipole()}
							,outdir 	= 'calcdir'	 # output directory
							,startingwfc= 'atomic+random' #alternatively, 'file' for restart
							,psppath 	= self.psppath()
							,output 	= {'removesave':True})


class NebJobMethods(object):
	#################
	# Required Params
	#################
	def nebfinal_pckl(self): return self.params['nebfinal_pckl']	

	######################
	# Required job methods
	#######################
	def strjob(self): raise NotImplementedError
	def jobEssentials(self): return set(['nebfinal_pckl'])
	def jobSummary(self): return {}

	def guessTime(self): return 20
	def wflow(self): raise NotImplementedError

	##############
	# Other methods
	##############
	def multiCalc(): 
		from espresso import multiespresso
		pass

class DosJobMethods(object):
	def strjob(self): return ''
	def jobEssentials(self): return set([])
	def jobSummary(self): return {}
	def guessTime(self): return 2
	def wflow(self): raise NotImplementedError


################
# Usable classes
################
###
class GPAWLatOptJob(JobBase,GPAWCalcJob,BulkJob,LatOptMethods): 
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'latticeopt'
		assert self.kind() 		== 'bulk'
		assert self.dftcode() 	== 'gpaw'

class QELatOptJob(	JobBase,QECalcJob,BulkJob,LatOptMethods): 
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'latticeopt'
		assert self.kind() 		== 'bulk'
		assert self.dftcode() 	== 'quantumespresso'

###
class GPAWBulkModJob(JobBase,GPAWCalcJob,BulkJob,BulkModMethods): 	
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'bulkmod'
		assert self.kind()		== 'bulk'
		assert self.dftcode() 	== 'gpaw'

class QEBulkModJob(	JobBase,QECalcJob,BulkJob,BulkModMethods): 	
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'bulkmod'
		assert self.kind()		== 'bulk'
		assert self.dftcode() 	== 'quantumespresso'
###
class BulkXCJob(JobBase,GPAWCalcJob,BulkJob,XCMethods): 	
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'xc'
		assert self.kind() 		== 'bulk'
		assert self.dftcode() 	== 'gpaw'

class SurfaceXCJob(JobBase,GPAWCalcJob,SurfaceJob,XCMethods): 	
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'xc'
		assert self.kind()		== 'surface'
		assert self.dftcode() 	== 'quantumespresso'

###
class QEBulkNebJob(JobBase,QECalcJob,BulkJob,NebJobMethods): 		
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'neb'
		assert self.kind()		== 'bulk'
		assert self.dftcode() 	== 'quantumespresso'


class QESurfaceNebJob(JobBase,QECalcJob,SurfaceJob,NebJobMethods): 		
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'neb'
		assert self.kind()		== 'surface'
		assert self.dftcode() 	== 'quantumespresso'

###
class GPAWBulkRelax(JobBase,GPAWCalcJob,BulkJob,RelaxJobMethods): 			
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'relax'
		assert self.kind() 		== 'bulk'
		assert self.dftcode() 	== 'gpaw'
class QEBulkRelax(JobBase,QECalcJob,BulkJob,RelaxJobMethods): 			
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'relax'
		assert self.kind() 		== 'bulk'
		assert self.dftcode() 	== 'quantumespresso'

class QEBulkVCRelax(JobBase,QECalcJob,BulkJob,VCRelaxMethods): 			
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'vcrelax'
		assert self.kind() 		== 'bulk'
		assert self.dftcode() 	== 'quantumespresso'
###
class GPAWSurfaceRelax(JobBase,GPAWCalcJob,SurfaceJob,RelaxJobMethods): 	
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'relax'
		assert self.kind() 		== 'surface'
		assert self.dftcode() 	== 'gpaw'
class QESurfaceRelax(JobBase,QECalcJob,SurfaceJob,RelaxJobMethods): 	
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'relax'
		assert self.kind() 		== 'surface'
		assert self.dftcode() 	== 'quantumespresso'
###
class GPAWMoleculeRelax(JobBase,GPAWCalcJob,MoleculeJob,RelaxJobMethods): 	
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'relax'
		assert self.kind() 		== 'molecule'
		assert self.dftcode() 	== 'gpaw'

class QEMoleculeRelax(JobBase,QECalcJob,MoleculeJob,RelaxJobMethods): 	
	def check(self):
		self.generalCheck()
		assert self.jobkind() 	== 'relax'
		assert self.kind()		== 'molecule'
		assert self.dftcode() 	== 'quantumespresso'

###
class QEDosJob(JobBase,QECalcJob,DosJobMethods):
	def check(self): 
		self.generalCheck()
		assert self.jobkind() 	== 'dos'
		assert self.dftcode() 	== 'quantumespresso'

class QEBulkVibJob(JobBase,QECalcJob,BulkJob,VibJobMethods): 	
	def check(self): 
		self.generalCheck()
		assert self.jobkind() 	== 'vib'
		assert self.kind() 		== 'bulk'
		assert self.dftcode() 	== 'quantumespresso'

class QESurfaceVibJob(JobBase,QECalcJob,SurfaceJob,VibJobMethods):
	def check(self): 
		self.generalCheck()
		assert self.jobkind() 	== 'vib'
		assert self.kind() 		== 'surface'
		assert self.dftcode() 	== 'quantumespresso'

class QEMoleculeVibJob(JobBase,QECalcJob,MoleculeJob,VibJobMethods): 
	def check(self): 
		self.generalCheck()
		assert self.jobkind() == 'vib';assert self.kind() == 'molecule';assert self.dftcode() == 'quantumespresso'

class GPAWBulkVibJob(JobBase,QECalcJob,BulkJob,VibJobMethods): 
	def check(self): 
		self.generalCheck()
		assert self.jobkind() 	== 'vib'
		assert self.kind() 		== 'bulk'
		assert self.dftcode() 	== 'gpaw'

class GPAWSurfaceVibJob(JobBase,QECalcJob,SurfaceJob,VibJobMethods): 
	def check(self): 
		self.generalCheck()
		assert self.jobkind() 	== 'vib'
		assert self.kind() 		== 'surface'
		assert self.dftcode() 	== 'gpaw'

class GPAWMoleculeVibJob(JobBase,QECalcJob,MoleculeJob,VibJobMethods): 
	def check(self): 
		self.generalCheck()
		assert self.jobkind() 	== 'vib'
		assert self.kind() 		== 'molecule'
		assert self.dftcode() 	== 'gpaw'

def assignJob(params):
	dftcode = params['dftcode']
	kind 	= params['kind']
	jobkind = params['jobkind']

	jobDict = 	{'quantumespresso':
					{'bulk':
						{'latticeopt':	QELatOptJob
						,'bulkmod':		QEBulkModJob
						,'dos':			QEDosJob
						,'neb':			QEBulkNebJob
						,'vib':			QEBulkVibJob
						,'relax': 		QEBulkRelax
						,'vcrelax': 	QEBulkVCRelax}
					,'surface':	
						{'vib': 		QESurfaceVibJob
						,'relax': 		QESurfaceRelax
						,'dos': 		QEDosJob
						,'neb': 		QESurfaceNebJob}
					,'molecule':
						{'vib': 		QEMoleculeVibJob
						,'relax': 		QEMoleculeRelax}}
				,'gpaw':
					{'bulk':
						{'latticeopt':	GPAWLatOptJob
						,'bulkmod':		GPAWBulkModJob
						,'xc':			BulkXCJob}
					,'surface':	
						{'xc':  		SurfaceXCJob
						,'relax': 		GPAWSurfaceRelax}
					,'molecule':
						{'vib':	 		GPAWMoleculeVibJob
						,'relax': 		GPAWMoleculeRelax
						,'xc':  		SurfaceXCJob}}}

	return jobDict[dftcode][kind][jobkind](params)

