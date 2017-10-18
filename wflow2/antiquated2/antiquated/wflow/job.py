import os,time,ast,random
from datetime import datetime
import numpy as np
from copy import deepcopy
from math import degrees,sqrt

from db  import getASEdb,query1all,query1,getCluster,updateStatus,updateDB

from ase  			import optimize
from ase.io    		import read
from ase.db     	import connect
from ase.parallel 	import parprint
from cluster 		import sherlock,sherlockgpaw,suncat,suncat2
"""
Defines the object of a bulkjob and some associated FUNCTIONS

TO DO:  Add DW_RATIO field (all jobs run so far, and in foreseeable future, are run at 10x

"Espresso executable doesn't seem to have been started"
	- caused by including /home/vossj/suncat/bin in PATH
"""

db = connect(getASEdb())


class Job(object):

	def __init__(self,jobid,jobkind,aseidinitial,vibids,nebids,xc,pw,kptden,psp,xtol,strain,convid,precalc,dftcode,comments,error,status):
		self.jobid 			= jobid 		# 0 	--- Int
		self.jobkind 		= jobkind 		# 1 	--- String: 'bulkrelax','relax','dos','vib','neb'
		self.aseidinitial 	= aseidinitial  # 2 	--- Int, initial structure for relaxation/dos/vib --- 
		self.vibids 		= vibids 		# 3 	--- [Int] only used for vib
		self.nebids 		= nebids 		# 4 	--- [Int] only used for neb
		self.xc 			= xc 			# 5 	--- String, 
		self.pw 			= pw			# 6 	--- Int
		self.kptden 		= kptden		# 7 	--- Double
		self.psp 			= psp			# 8 	--- String
		self.xtol 			= xtol			# 9 	--- Double
		self.strain 		= strain		# 10 	--- Double
		self.convid 		= convid		# 11 	--- Int
		self.precalc 		= precalc		# 12 	--- String
		self.dftcode 		= dftcode		# 13 	--- String
		self.comments 		= comments		# 14 	--- String
		self.error 			= error			# 15 	--- String
		self.status 		= status		# 16 	--- String

	def row(self): 			return db.get(self.aseidinitial)
	def atoms(self): 		return db.get_atoms(id=self.aseidinitial)
	def symbols(self): 		return self.atoms().get_chemical_symbols()		
	def cellinit(self): 	return self.atoms().get_cell()
	def paramsinit(self): 	return cell2param(self.cellinit())
	def posinit(self):	 	return self.atoms().get_scaled_positions()
	def magmomsinit(self):	return self.atoms().get_initial_magnetic_moments()
	def tags(self): 		return self.atoms().get_tags()
	def constraints(self): 	return self.atoms().constraints
	### Defined for all
	def name(self): 		return self.row().get('name')
	def relaxed(self): 		return self.row().get('relaxed')
	def emt(self): 			return self.row().get('emt')
	def comments(self): 	return self.row().get('comments')
	def kind(self): 		return self.row().get('kind') #molecule/bulk/surface
	### Defined for bulk only
	def structure(self): 	return self.row().get('structure')
	def vacancies(self): 	return self.row().get('vacancies')
	### Defined for surface only
	def sites(self): 		return self.row().get('sites')
	def facet(self): 		return self.row().get('facet')
	def xy(self): 			return self.row().get('xy')
	def layers(self): 		return self.row().get('layers')
	def constrained(self): 	return self.row().get('constrained')
	def symmetric(self): 	return self.row().get('symmetric')
	def vacuum(self): 		return self.row().get('vacuum')
	def adsorbates(self): 	return self.row().get('adsorbates')

	def kpt(self): 
		kpt = kptdensity2monkhorstpack(self.atoms(),self.kptden)
		if   self.kind()=='surface': 	return kpt[:2]+[1]
		elif self.kind()=='molecule': 	return [1,1,1]
		else: 							return kpt

	def convRow(self): 	return query1all('convergenceid',self.convid,'convergence')
	def dwratio(self): 	return self.convRow()[1]
	def econv(self): 	return self.convRow()[2]
	def mixing(self): 	return self.convRow()[3]
	def nmix(self): 	return self.convRow()[4]
	def maxstep(self): 	return self.convRow()[5]
	def nbands(self): 	return self.convRow()[6]
	def sigma(self): 	return self.convRow()[7]
	def fmax(self): 	return self.convRow()[8]
	def delta(self): 	return self.convRow()[9] 	# for use in vibration analysis
	def climb(self): 	return self.convRow()[10] 	# for use in NEB analysis
	
	######################################################
	# PRINT INFORMATION
	#######################################################
	def __str__(self):  return "Job_%d"%(self.jobid)
	def summary(self):  return 	("\nJob ID: %d"%self.jobid
								+"\nKind: %s"%self.jobkind
								+"\nname: %s"%('' if self.name() is None else self.name())
								+"\nfacet: %d"%(0 if self.facet() is None else self.facet())
								+"\nxy: %d"%(0 if self.xy() is None else self.xy())
								+"\nlayers: %d"%(0 if self.layers() is None else self.layers())
								+"\nxc: %s"%self.xc
								+"\npw: %d"%self.pw
								+"\nkptden: %f"%self.kptden
								+"\ndftcode: %s"%self.dftcode
								+"\nerror: %s"%self.error
								+"\nstatus: %s"%self.status)

	def showSites(self):
		if self.kind()=='surface':
			with open('/scratch/users/ksb/img/sites.png','wb') as f:
				f.write(self.sites().decode('base64'))
			os.system('display /scratch/users/ksb/img/sites.png &')
		else: print "no sites to show. not a molecule! (%s)"%self.kind()
	#########################################################
	# JOB SUBMISSION
	#######################################################
	def allocate(self,n):
		"""Philosophies:
			- short jobs should go to suncat/suncat2 
				- more open spots on queue for suncat generally, but jobs run slower
				http://www.nersc.gov/users/computational-systems/edison/running-jobs/httpwww-nersc-govuserscomputational-systemsedisonrunning-jobsqueues-and-policies/
				http://www.nersc.gov/users/computational-systems/cori/running-jobs/queues-and-policies/
		"""
		
		def quantumespressoAllocate(randomx):
			if   randomx < 0: 		return sherlock
			elif randomx < 0.5: 	return suncat
			elif randomx < 1: 		return suncat2

		if self.dftcode=='gpaw': 	
			return [sherlockgpaw]*n  
		else:		
			rand=[]
			for i in range(n-1): rand.append(random.random())
			return [quantumespressoAllocate(x) for x in rand]+[sherlock] 

	def guessTime(self,jobkind): 
		baseGuess = {'OptimizeLattice': 	4 * len(self.iGuess()) 	if self.jobkind =='bulkrelax' 	else 0
					,'GetBulkModulus':		4 						if self.jobkind =='bulkrelax' 	else 0
					,'Relax': 				4 						if self.jobkind =='relax' 		else 0
					,'Vibrations': 			2 * len(self.vibids) 	if self.jobkind == 'vib' 		else 0
					,'GetXCcontribs': 		3 						if self.dftcode == 'gpaw'		else 0.05
					,'SaveResults': 		0.1}
		xcFactor = 3 if self.xc == 'mBEEF' else 1
		
		return baseGuess[jobkind]*xcFactor
	
	def submit(self): 
		from fireworks  import LaunchPad

		launchpad = LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')

		if 		self.jobkind == 'bulkrelax': 	wflow = self.submitBulkRelax()
		elif 	self.jobkind == 'relax': 		wflow = self.submitRelax()
		elif 	self.jobkind == 'vib': 			wflow = self.submitVib()
		elif 	self.jobkind == 'neb': 			wflow = self.submitNEB()
		elif 	self.jobkind == 'dos': 			wflow = self.submitDOS()

		print "Submitting job with ID = %d"%self.jobid

		updateDB('status','jobid',self.jobid,'queued',None,'job')
		launchpad.add_wf(wflow)

		"""if query1('status','jobid',self.jobid,'job')=='initialized': 	updateStatus(self.jobid,'initialized','queued')		"""

	def submitRelax(self): 
		from fireworks   		import Firework,Workflow
		from standardScripts 	import Relax,GetXCcontribs,SaveResults
		clusters = self.allocate(3) ; print 'submitting to ',clusters
		timestamp = '_'+datetime.now().strftime('%Y_%m_%d_%H_%M')
		names = [x+'_%d'%self.jobid for x in ['Relax','GetXCcontribs','SaveResults']]
		fw1 = Firework([Relax()]
						,name = names[0]
						,spec = {"jobID": 			self.jobid
								,'_fworker': 		clusters[0].fworker
								,"_pass_job_info": 	True
								,"_files_out": 		{"fw1":"inp.gpw"}
								,"_queueadapter": 	clusters[0].qfunc(self.guessTime('Relax'))
								,"_launch_dir": 		clusters[0].launchdir+names[0]+timestamp})
		fw2 = Firework(	[GetXCcontribs()]
						,name 		= names[1]
						,parents 	= [fw1]
						,spec 		= 	{"jobID": 			self.jobid
										,'_fworker': 		clusters[1].fworker
										,"_pass_job_info": 	True
										,"_files_in": 		{"fw1":"inp.gpw"}
										,"_queueadapter": 	clusters[1].qfunc(self.guessTime('GetXCcontribs'))
										,"_launch_dir": 		clusters[1].launchdir+names[1]+timestamp})

		fw3 = Firework(	[SaveResults()]
						,name 		= names[2]
						,parents 	= [fw1,fw2]
						,spec 		= 	{"jobID": 			self.jobid
										,'_fworker': 		clusters[2].fworker #MUST be sherlock
										,"_queueadapter": 	clusters[2].qfunc(self.guessTime('SaveResults'))
										,"_launch_dir": 	clusters[2].launchdir+names[2]+timestamp})

		return Workflow([fw1,fw2,fw3],name='BulkRelaxation_%d'%self.jobid)

	def submitNEB(self): raise NotImplementedError
	def submitDOS(self): raise NotImplementedError

	def submitVib(self): 
		from fireworks   		import Firework,Workflow
		from standardScripts 	import Vibrations,SaveResults
		clusters = self.allocate(2) 
		timestamp = '_'+datetime.now().strftime('%Y_%m_%d_%H_%M')
		names = [x+'_%d'%self.jobid for x in ['Vibrations','SaveResults']]
		fw1 = Firework(	[OptimizeLattice()]
						,name = names[0]
						,spec = {"jobID": 			self.jobid
								,'_fworker': 		clusters[0].fworker
								,"_pass_job_info": 	True
								,"_files_out": 		{"fw1":"vibrations.txt"}
								,"_queueadapter": 	clusters[0].qfunc(self.guessTime('Vibrations'))
								,"_launch_dir": 		clusters[0].launchdir+names[0]+timestamp})


		fw2 = Firework(	[SaveResults()]
						,name 		= 	names[1]
						,parents 	= 	[fw1,fw2,fw3]
						,spec 		= 	{"jobID": 			self.jobid
										,'_fworker': 		clusters[1].fworker #MUST be sherlock
										,"_files_in": 		{"fw1":"vibrations.txt"}
										,"_queueadapter": 	clusters[1].qfunc(self.guessTime('SaveResults'))
										,"_launch_dir": 	clusters[1].launchdir+names[1]+timestamp})

		return Workflow([fw1,fw2],name='Vibration_%d'%self.jobid)


	def submitBulkRelax(self): 
		from fireworks   		import Firework,Workflow
		from standardScripts 	import OptimizeLattice,GetBulkModulus,GetXCcontribs,SaveResults
		clusters = self.allocate(4) 
		timestamp = '_'+datetime.now().strftime('%Y_%m_%d_%H_%M')
		
		names = [x+'_%d'%self.jobid for x in ['OptimizeLattice','GetBulkModulus','GetXCcontribs','SaveResults']]

		fw1 = Firework(	[OptimizeLattice()]
						,name = names[0]
						,spec = {"jobID": 			self.jobid
								,'_fworker': 		clusters[0].fworker
								,"_pass_job_info": 	True
								,"_files_out": 		{"fw1":"inp.gpw"}
								,"_queueadapter": 	clusters[0].qfunc(self.guessTime('OptimizeLattice'))
								,"_launch_dir": 		clusters[0].launchdir+names[0]+timestamp})

		fw2 = Firework(	[GetBulkModulus()]
						,name 		= names[1]
						,parents	= [fw1]
						,spec		= 	{"jobID": 			self.jobid
										,'_fworker': 		clusters[1].fworker
										,"_pass_job_info": 	True
										,"_queueadapter": 	clusters[1].qfunc(self.guessTime('GetBulkModulus'))
										,"_launch_dir": 		clusters[1].launchdir+names[1]+timestamp})

		fw3 = Firework(	[GetXCcontribs()]
						,name 		= names[2]
						,parents 	= [fw1]
						,spec 		= 	{"jobID": 			self.jobid
										,'_fworker': 		clusters[2].fworker
										,"_pass_job_info": 	True
										,"_files_in": 		{"fw1":"inp.gpw"}
										,"_queueadapter": 	clusters[2].qfunc(self.guessTime('GetXCcontribs'))
										,"_launch_dir": 		clusters[2].launchdir+names[2]+timestamp})

		fw4 = Firework(	[SaveResults()]
						,name 		= names[3]
						,parents 	= [fw1,fw2,fw3]
						,spec 		= 	{"jobID":self.jobid
										,'_fworker':clusters[3].fworker #MUST be sherlock
										,"_queueadapter": 	clusters[3].qfunc(self.guessTime('SaveResults'))
										,"_launch_dir": 		clusters[3].launchdir+names[3]+timestamp})


		return Workflow([fw1,fw2,fw3,fw4],name='BulkRelaxation_%d'%self.jobid)



	##############
	# BULK RELATED
	##############

	def fromParams(self,params): 
		""" 
		Params is a list of between 1 and 6 numbers (a,b,c,alpha,beta,gamma). 
		Depending on structure, we can construct the cell from a subset of these parameters
		"""
		atoms = self.atoms().copy()
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
		elif self.structure() =='triclinic': cell = params
		else: raise NotImplementedError, 'job.fromParams() cannot handle unknown structure = '+self.structure()

		atoms.set_cell(cell,scale_atoms=True)
		return atoms

	def iGuess(self):
		a,b,c,alpha,beta,gamma = self.paramsinit()
		if self.structure() in ['fcc','bcc','rocksalt','diamond','cscl','zincblende']: return [a]
		elif self.structure() in ['hexagonal']: return [a,c]
		elif self.structure() in ['triclinic']: return self.paramsinit()
		else: raise ValueError, 'Bad entry in "structure" field for Atoms object info dictionary: '+self.structure()

	####################
	# CALCULATOR RELATED
	####################
	def gpawRestart(self):
		from gpaw import restart,PW,Davidson,Mixer,MixerSum,FermiDirac

		spinpol = any([x>0 for x in self.magmomsinit()])

		return restart('preCalc_inp.gpw',mode    = PW(self.pw)
										,xc      = self.xc
										,kpts    = self.kpt()
										,spinpol = spinpol
										,convergence = {'energy':self.econv()} #eV/electron
										,mixer =  ((MixerSum(beta=self.mixing(),nmaxold=self.nmix(),weight=100)) 
												if spinpol else (Mixer(beta=self.mixing(),nmaxold=self.nmix(),weight=100)))
										,maxiter       = self.maxstep()
										,nbands        = self.nbands()
										,occupations   = FermiDirac(self.sigma())
										,setups        = self.psp #(pspDict[calc.psp]).pthFunc[cluster]
										,eigensolver   = Davidson(5)
										,poissonsolver = None # {'dipolelayer': 'xy'} if isinstance(self,SurfaceJob) else None
										,txt='%d_%s.txt'%(self.pw,self.xc)
										,symmetry={'do_not_symmetrize_the_density': True}) #ERROR IN LI bcc 111 surface
										

	def gpawCalc(self,xc,spinpol):
		from gpaw import GPAW,PW,Davidson,Mixer,MixerSum,FermiDirac
		return GPAW(mode         = PW(self.pw)                         #eV
					,xc          = xc
					,kpts        = self.kpt()
					,spinpol     = spinpol
					,convergence = {'energy':self.econv()} #eV/electron
					,mixer       = ((MixerSum(beta=self.mixing(),nmaxold=self.nmix(),weight=100)) 
									if spinpol else (Mixer(beta=self.mixing(),nmaxold=self.nmix(),weight=100)))
					,maxiter       = self.maxstep()
					,nbands        = self.nbands()
					,occupations   = FermiDirac(self.sigma())
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
						,dw      = self.pw*self.dwratio()
						,xc      = xc
						,kpts    = self.kpt()
						,spinpol = spinpol
						,convergence = {'energy':self.econv()
										,'mixing':self.mixing()
										,'nmix':self.nmix()
										,'maxsteps':self.maxstep()
										,'diag':'david'}
						,nbands = self.nbands()
						,sigma  = self.sigma()
						,dipole = {'status': True if self.kind()=='surface' else False}
						,outdir = 'calcdir'	 # output directory
						,startingwfc='file' if restart else 'atomic+random'
						,psppath= psppath
						,output = {'removesave':True})

	def calc(self,precalc=False,restart=False):
		""" Creates an ASE calculator object """

		xc = self.precalc if precalc else self.xc
		spinpol = any([x>0 for x in self.magmomsinit()])

		if   self.dftcode == 'gpaw':            return self.gpawCalc(xc,spinpol)
		elif self.dftcode == 'quantumespresso': return self.qeCalc(xc,spinpol,restart)
		else: raise ValueError, 'Unknown dftcode '+self.dftcode

	def vibCalc(self):
		from espresso.vibespresso import vibespresso

		spinpol = any([x>0 for x in self.magmomsinit()])
		dipole 	= True if self.kind()=='surface' else False

		return  vibespresso(pw=self.pw        				 # planewave cutoff
				            ,dw=self.pw*self.dwratio()     # density cutoff
							,xc      = self.xc
							,kpts    = self.kpt()
							,spinpol = spinpol
							,convergence = {'energy':self.econv()
											,'mixing':self.mixing()
											,'nmix':self.nmix()
											,'maxsteps':self.maxstep()
											,'diag':'david'}
							,nbands = self.nbands()
							,sigma  = self.sigma()
							,dipole = {'status': dipole}
					        ,outdirprefix = 'vibdir'
							,output = {'removesave':True}
							,psppath= psppath)

	def optimizePos(self,atoms,calc,saveWF=False):
		atoms.set_calculator(calc)

		maxForce = np.amax(abs(atoms.get_forces()))
		if maxForce > self.fmax:
			parprint("max force = %f, optimizing positions"%(maxForce))
			dyn = optimize.BFGS(atoms=atoms, logfile='qn.log', trajectory='qn.traj',restart='qn.pckl')
			dyn.run(fmax=self.fmax())

		if saveWF and self.dftcode=='gpaw': 
			atoms.calc.write('preCalc_inp.gpw', mode='all') #for use in getXCContribs
			atoms,calc = self.gpawRestart()

	########################################################################
	# SQL METHODS
	########################################################################
	def sqlTable(self,tentative=False):  return 'tentativejob' if tentative else 'job'

	def sqlCols(self):   return ['jobkind'
								,'aseidinitial'
								,'vibids'
								,'nebids'
								,'xc'
								,'pw'
								,'kptden'
								,'psp'
								,'xtol'
								,'strain'
								,'convid'
								,'precalc'
								,'dftcode'
								,'comments'
								,'error'
								,'status']

	def sqlInsert(self): return [self.jobkind
								,self.aseidinitial
								,self.vibids
								,self.nebids
								,self.xc
								,self.pw
								,self.kptden
								,self.psp
								,self.xtol
								,self.strain
								,self.convid
								,self.precalc
								,self.dftcode
								,self.comments
								,self.error
								,self.status]

	def sqlEq(self): return ['jobkind'
							,'aseidinitial'
							,'vibids'
							,'nebids'
							,'xc'
							,'pw'
							,'kptden'
							,'psp'
							,'xtol'
							,'strain'
							,'convid'
							,'precalc'
							,'dftcode']



##############################################################################
# OTHER FUNCTIONS
##############################################################################
def db2object(i,table='job',obj='job'): 
	row   = query1all('_rowid_',i,table)
	if 		obj == 'job':	 return Job(*row)
	elif 	obj == 'result': return Result(*row)
	elif 	obj == 'convergence': return Convergence(*row)

def cell2param(cell):
	def angle(v1,v2): 
		return degrees(np.arccos(np.dot(v1,np.transpose(v2))/(np.linalg.norm(v1)*np.linalg.norm(v2))))
	a = np.linalg.norm(cell[0])
	b = np.linalg.norm(cell[1])
	c = np.linalg.norm(cell[2])
	alpha = angle(cell[1],cell[2])
	beta  = angle(cell[0],cell[2])
	gamma = angle(cell[0],cell[1])
	return (a,b,c,alpha,beta,gamma)

def kptdensity2monkhorstpack(atoms, kptdensity, even=True):
    """Convert k-point density to Monkhorst-Pack grid size.
    atoms: Atoms object  ---  Contains unit cell and information about boundary conditions.
    kptdensity: float    ---  Required k-point density.  Default value is 3.5 point per Ang^-1.
    even: bool           ---  Round up to even numbers.
    """
    recipcell = atoms.get_reciprocal_cell()
    kpts = []
    for i in range(3):
        if atoms.pbc[i]:
            k = 2 * 3.14159 * sqrt((recipcell[i]**2).sum()) * kptdensity
            if even:
                kpts.append(2 * int(np.ceil(k / 2)))
            else:
                kpts.append(int(np.ceil(k)))
        else:
            kpts.append(1)
    return kpts

