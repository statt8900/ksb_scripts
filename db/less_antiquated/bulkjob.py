import os,time,ast,random
from datetime import datetime
import numpy as np
from copy import deepcopy
from math import degrees

from db  import getASEdb,query1all,getCluster,updateStatus

from ase  			import optimize
from ase.io    		import read
from ase.db     	import connect
from ase.parallel 	import parprint
from cluster 		import sherlock,suncat,suncat2
from job 			import kptdensity2monkhorstpack
"""
Defines the object of a bulkjob and some associated FUNCTIONS

TO DO:  Add DW_RATIO field (all jobs run so far, and in foreseeable future, are run at 10x

"""

class BulkJob(object):

	def __init__(self
				,bulkid
				,xc,psp,pw,kptden,kpt,econv,mixing,nmix,maxstep,nbands,sigma
				,fmax,xtol,strain,precalcxc,dftcode
				,ID=None,created=time.time(),createdby='Kris',lastmodified=time.time(),kind='relax',comments=None,status='initialized'):
		
		"""Initialize a bulk optimization job."""

		db                = connect(getASEdb())
		row               = db.get(bulkid)

		self.bulkid       = bulkid    # Bulk
		
		atoms             = db.get_atoms(id=bulkid)
		self.atoms        = atoms
		
		self.symbols      = atoms.get_chemical_symbols()
		cell              = atoms.get_cell()
		self.cell         = cell
		params            = cell2param(cell)
		self.a            = params[0]
		self.b            = params[1]
		self.c            = params[2]
		self.alpha        = params[3]
		self.beta         = params[4]
		self.gamma        = params[5]
		self.positions    = atoms.get_positions()
		self.magmoms      = atoms.get_initial_magnetic_moments()
		self.tags         = atoms.get_tags()
		self.constraints  = atoms.constraints
		self.bulkname     = row.get('name')
		self.structure    = row.get('structure')
		self.emt          = row.get('emt')
		self.bulkcomments = row.get('comments')
		self.kind 			= row.get('kind')

		self.xc           = xc            # String
		self.psp          = psp           # String
		self.pw           = pw            # Int
		self.kptden       = kptden

		if kpt is None: self.kpt = kptdensity2monkhorstpack(atoms,kptden)
		else:           self.kpt = ast.literal_eval(kpt)         

		self.econv        = econv         # Float
		self.mixing       = mixing        # Float
		self.nmix         = nmix          # Int
		self.maxstep      = maxstep       # Float
		self.nbands       = nbands        # Int
		self.sigma        = sigma         # Float
		self.fmax         = fmax          # Float
		self.xtol         = xtol          # Float
		self.strain       = strain        # Float 
		self.precalcxc    = precalcxc     # String
		self.dftcode      = dftcode       # String
		self.ID           = ID            # Int
		self.created      = created       # Double
		self.createdby    = createdby     # String
		self.lastmodified = lastmodified  # Double
		self.kind         = kind          # String
		self.comments     = comments      # String
		self.status       = status        # String

	def __str__(self):  return "Job#%d"%(self.ID)
			
	def printFields(self):	print """\t\tbulkid=%d,bulkname=%s,symbols=%s,xc=%s
		,psp=%s,pw=%d,kptden=%d,kpt=%s,econv=%f,mixing=%f,nmix=%f,maxstep=%f,nbands=%d,sigma=%f
		,fMax=%f,xTol=%f,strain=%f,precalcxc=%s,dftCode=%s
		,ID=%s,created=%s,createdby=%s,lastmodified=%s,kind=%s,comments=%s,status=%s):"""%(self.bulkid,self.bulkname,str(self.symbols)
		,self.xc,self.psp,self.pw,self.kptden,self.kpt,self.econv,self.mixing,self.nmix,self.maxstep,self.nbands,self.sigma
		,self.fmax,self.xtol,self.strain,self.precalcxc,self.dftcode
		,self.ID,self.created,self.createdby,self.lastmodified,self.kind,self.comments,self.status)

	def allocate(self):
		"""Philosophies:
			- short jobs should go to suncat/suncat2 
				- more open spots on queue for suncat generally, but jobs run slower

http://www.nersc.gov/users/computational-systems/edison/running-jobs/httpwww-nersc-govuserscomputational-systemsedisonrunning-jobsqueues-and-policies/
http://www.nersc.gov/users/computational-systems/cori/running-jobs/queues-and-policies/
			"""

		def subAllocate(randomx):
			if   randomx < 0.2: 	return sherlock
			elif randomx < 0.6: 	return suncat
			else: 					return suncat2
		
		rand = [subAllocate(rx) for rx in [random.random(),random.random(),random.random()]]+[sherlock]
		
		return  [sherlock]*4 if self.dftcode=='gpaw' else rand 


	def guessTime(self,jobkind): 
		mbeef = self.xc == 'mBEEF'
		if jobkind == 'optimizelattice':
			if mbeef: 	return 15
			else: 		return 10	
		elif jobkind == 'getbulkmodulus':
			if mbeef: 	return 15
			else: 		return 10	
		elif jobkind == 'getxccontribs': 	return 10
		elif jobkind == 'saveresults': 		return 1
		else: raise ValueError, 'unknown Job Kind '+jobkind
	
	def submit(self): 
		from fireworks   		import LaunchPad,Firework,Workflow,PyTask
		from standardScripts 	import OptimizeLattice,GetBulkModulus,GetXCcontribs,SaveResults

		launchpad=LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')

		names = [x+'_%d'%self.ID for x in ['OptimizeLattice','GetBulkModulus','GetXCcontribs','SaveResults']]
	
		clusters = self.allocate() # 
		timestamp = '_'+datetime.now().strftime('%Y_%m_%d_%H_%M')

		fw1 = Firework(	[OptimizeLattice()]
						,name = names[0]
						,spec = {"jobID": 			self.ID
								,'_fworker': 		clusters[0].fworker
								,"_pass_job_info": 	True
								,"_files_out": 		{"fw1":"inp.gpw"}
								,"_queueadapter": 	clusters[0].qfunc(self.guessTime('optimizelattice'))
								,"_launch_dir": 		clusters[0].launchdir+names[0]+timestamp})

		fw2 = Firework(	[GetBulkModulus()]
						,name 		= names[1]
						,parents	= [fw1]
						,spec		= 	{"jobID": 			self.ID
										,'_fworker': 		clusters[1].fworker
										,"_pass_job_info": 	True
										,"_queueadapter": 	clusters[1].qfunc(self.guessTime('getbulkmodulus'))
										,"_launch_dir": 		clusters[1].launchdir+names[1]+timestamp})


		fw3 = Firework(	[GetXCcontribs()]
						,name 		= names[2]
						,parents 	= [fw1]
						,spec 		= 	{"jobID": 			self.ID
										,'_fworker': 		clusters[2].fworker
										,"_pass_job_info": 	True
										,"_files_in": 		{"fw1":"inp.gpw"}
										,"_queueadapter": 	clusters[2].qfunc(self.guessTime('getxccontribs'))
										,"_launch_dir": 		clusters[2].launchdir+names[2]+timestamp})

		fw4 = Firework(	[SaveResults()]
						,name 		= names[3]
						,parents 	= [fw1,fw2,fw3]
						,spec 		= 	{"jobID":self.ID
										,'_fworker':clusters[3].fworker #MUST be sherlock
										,"_queueadapter": 	clusters[3].qfunc(self.guessTime('saveresults'))
										,"_launch_dir": 		clusters[3].launchdir+names[3]+timestamp})


		wflow=Workflow([fw1,fw2,fw3,fw4],name='BulkRelaxation_%d'%self.ID)

		print "Submitting job with ID = %d"%self.ID

		updateStatus('bulkjob',self.ID,'initialized','queued')

		launchpad.add_wf(wflow)
		

	##############
	# BULK RELATED
	##############

	def fromParams(self,params): 
		atoms = deepcopy(self.atoms)
		a = params[0]
		if   self.structure in ['cubic','cscl']:     
			cell = [a,a,a,90,90,90]
		elif self.structure in ['fcc','diamond','zincblende','rocksalt']:  
			cell = [a,a,a,60,60,60]
		elif self.structure in ['bcc']:
			alpha = 109.47122
			cell = [a,a,a,alpha,alpha,alpha]
		elif self.structure in ['hcp','hexagonal']: 
			cell = [a,a,params[1],90,90,120]
		elif self.structure =='triclinic': cell = params
		else: raise NotImplementedError

		atoms.set_cell(cell,scale_atoms=True)
		return atoms

	def iGuess(self):
		if self.structure in ['fcc','bcc','rocksalt','diamond','cscl','zincblende']:
			return [self.a]
		elif self.structure in ['hexagonal']:      
			return [self.a,self.c]
		elif self.structure in ['triclinic']:      
			return [self.a,self.b,self.c,self.alpha,self.beta,self.gamma]
		else: raise ValueError, 'Bad entry in "structure" field for Atoms object info dictionary: '+self.structure()

	####################
	# CALCULATOR RELATED
	####################
	def gpawRestart(self):
		from gpaw import restart,PW,Davidson,Mixer,MixerSum,FermiDirac

		spinpol = any([x>0 for x in self.magmoms])

		return restart('preCalc_inp.gpw',mode    = PW(self.pw)
										,xc      = self.xc
										,kpts    = self.kpt
										,spinpol = spinpol
										,convergence = {'energy':self.econv} #eV/electron
										,mixer =  ((MixerSum(beta=self.mixing,nmaxold=self.nmix,weight=100)) 
												if spinpol else (Mixer(beta=self.mixing,nmaxold=self.nmix,weight=100)))
										,maxiter       = self.maxstep
										,nbands        =  -1*self.nbands
										,occupations   = FermiDirac(self.sigma)
										,setups        = self.psp #(pspDict[calc.psp]).pthFunc[cluster]
										,eigensolver   = Davidson(5)
										,poissonsolver = None # {'dipolelayer': 'xy'} if isinstance(self,SurfaceJob) else None
										,txt='%d_%s.txt'%(self.pw,self.xc)
										,symmetry={'do_not_symmetrize_the_density': True}) #ERROR IN LI bcc 111 surface
										

	def gpawCalc(self,xc,spinpol):
		from gpaw import GPAW,PW,Davidson,Mixer,MixerSum,FermiDirac
		return GPAW(mode         = PW(self.pw)                         #eV
					,xc          = xc
					,kpts        = self.kpt
					,spinpol     = spinpol
					,convergence = {'energy':self.econv} #eV/electron
					,mixer       = ((MixerSum(beta=self.mixing,nmaxold=self.nmix,weight=100)) 
									if spinpol else (Mixer(beta=self.mixing,nmaxold=self.nmix,weight=100)))
					,maxiter       = self.maxstep
					,nbands        =  -1*self.nbands
					,occupations   = FermiDirac(self.sigma)
					,setups        = 'sg15' 
					,eigensolver   = Davidson(5)
					,poissonsolver = None # {'dipolelayer': 'xy'} if isinstance(self,SurfaceJob) else None
					,txt='%d_%s.txt'%(self.pw,xc)
					,symmetry={'do_not_symmetrize_the_density': True}) #ERROR IN LI bcc 111 surface
					
	def qeCalc(self,xc,spinpol,restart):
		from espresso import espresso	
				
		pspDict = 	{'sherlock':
						{'gbrv15pbe':'/home/vossj/suncat/psp/gbrv1.5pbe'}
					,'suncat':
						{'gbrv15pbe':'/nfs/slac/g/suncatfs/sw/external/esp-psp/gbrv1.5pbe'}}

		cluster = getCluster()
		psppath = pspDict[cluster][self.psp]

		return espresso( pw      = self.pw
						,dw      = self.pw*10
						,xc      = xc
						,kpts    = self.kpt
						,spinpol = spinpol
						,convergence = {'energy':self.econv
										,'mixing':self.mixing
										,'nmix':self.nmix
										,'maxsteps':self.maxstep
										,'diag':'david'}
						,nbands = -1*self.nbands
						,sigma  = self.sigma
						,dipole = {'status': False}
						,outdir = 'calcdir'	 # output directory
						,startingwfc='file' if restart else 'atomic+random'
						,psppath= psppath
						,output = {'removesave':True})

	def calc(self,precalc=False,restart=False):
		""" Creates an ASE calculator object """

		xc = self.precalcxc if precalc else self.xc
		spinpol = any([x>0 for x in self.magmoms])

		if   self.dftcode == 'gpaw':            return self.gpawCalc(xc,spinpol)
		elif self.dftcode == 'quantumespresso': return self.qeCalc(xc,spinpol,restart)
		else: raise ValueError, 'Unknown dftcode '+self.dftcode

	def optimizePos(self,atoms,calc,saveWF=False):
		atoms.set_calculator(calc)

		maxForce = np.amax(abs(atoms.get_forces()))
		if maxForce > self.fmax:
			parprint("max force = %f, optimizing positions"%(maxForce))
			dyn = optimize.BFGS(atoms=atoms, logfile='qn.log', trajectory='qn.traj',restart='qn.pckl')
			dyn.run(fmax=self.fMax)

		if saveWF and self.dftcode=='gpaw': 
			atoms.calc.write('preCalc_inp.gpw', mode='all') #for use in getXCContribs
			atoms,calc = self.gpawRestart()



	########################################################################
	########################################################################
	# SQL METHODS
	########################################################################
	def sqlTable(self,tentative=False):  return 'tentative' if tentative else 'bulkjob'

	def sqlCols(self):   return ['bulkid'
								,'symbols'
								,'cell'
								,'a','b','c','alpha','beta','gamma'
								,'positions','magmoms','tags','constraints'
								,'bulkname','structure','emt','bulkcomments'
								,'xc','psp','pw','kptden','kpt','econv','mixing','nmix','maxstep','nbands','sigma'
								,'fmax','xtol','strain'
								,'precalcxc','dftcode'
								,'created','createdby','lastmodified','kind','comments'
								,'status']

	def sqlInsert(self): return [str(self.bulkid)
								,str(self.symbols)
								,str(self.cell)
								,str(self.a)
								,str(self.b)
								,str(self.c)
								,str(self.alpha)
								,str(self.beta)
								,str(self.gamma)
								,str(self.positions)
								,str(self.magmoms)
								,str(self.tags)
								,str(self.constraints)
								,self.bulkname
								,self.structure
								,str(self.emt)
								,self.bulkcomments
								,self.xc
								,self.psp
								,str(self.pw)
								,str(self.kptden)
								,str(self.kpt)
								,str(self.econv)
								,str(self.mixing)
								,str(self.nmix)
								,str(self.maxstep)
								,str(self.nbands)
								,str(self.sigma)
								,str(self.fmax)
								,str(self.xtol)
								,str(self.strain)
								,self.precalcxc
								,self.dftcode
								,self.created
								,self.createdby
								,self.lastmodified
								,self.kind
								,self.comments
								,self.status]

	def sqlEq(self): return ['bulkid'
							,'xc'
							,'psp'
							,'pw'
							,'kptden'
							,'kpt'
							,'econv'
							,'mixing'
							,'nmix'
							,'maxstep'
							,'nbands'
							,'sigma'
							,'fmax'
							,'xtol'
							,'strain'
							,'precalcxc'
							,'dftcode']

##############################################################################
##############################################################################
# OTHER FUNCTIONS
##############################################################################
def db2BulkJob(i,tentative=False): 
	table = 'tentative' if tentative else 'bulkjob'
	row             = query1all(table,'id',i)
	immutableThings = (row[1],)+row[18:34]
	mutableThings   = (row[0],)+row[35:]
	return BulkJob(*(immutableThings+mutableThings))



def cell2param(cell):
	def angle(v1,v2): 
		return degrees(np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))))
	a = np.linalg.norm(cell[0])
	b = np.linalg.norm(cell[1])
	c = np.linalg.norm(cell[2])
	alpha = angle(cell[1],cell[2])
	beta  = angle(cell[0],cell[2])
	gamma = angle(cell[0],cell[1])
	return (a,b,c,alpha,beta,gamma)




