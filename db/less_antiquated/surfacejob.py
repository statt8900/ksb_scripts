import os,time,ast,random,base64
import numpy as np
from datetime 	import datetime
from copy 		import deepcopy
from math 		import degrees,sqrt

from db  		import getASEdb,query1all,getCluster,updateStatus

from ase  			import optimize,Atoms
from ase.io    		import read
from ase.db     	import connect
from ase.parallel 	import parprint
from cluster 		import sherlock,suncat,suncat2,suncattest
from job 			import Job,kptdensity2monkhorstpack


"""
Defines the object of a surface job and some associated FUNCTIONS
"""

class SurfaceJob(Job):

	def __init__(self
				,surfid
				,xc,psp,pw,dwratio,kptden,kpt,econv,mixing,nmix,maxstep,nbands,sigma
				,fmax,precalcxc,dftcode
				,ID=None,created=time.time(),createdby='Kris',lastmodified=time.time(),kind='surfrelax',comments=None,status='initialized'):
		
		"""Initialize a surface optimization job."""

		db                	= connect(getASEdb())
		row               	= db.get(surfid)

		self.surfid      	= surfid    # Bulk	
		
		atoms         	  	= db.get_atoms(id=surfid)
		self.atoms        	= atoms
		
		self.symbols      	= atoms.get_chemical_symbols()
		self.cell           = atoms.get_cell()
		self.positions    	= atoms.get_positions()
		self.magmoms      	= atoms.get_initial_magnetic_moments()
		self.tags         	= atoms.get_tags()
		self.constraints  	= atoms.constraints
		self.surfname     	= row.get('name')
		self.structure     	= row.get('structure')
		self.sites     		= row.get('sites') 			# String version of picture
		self.facet 		  	= row.get('facet')
		self.xy 			= row.get('xy')     		#(Int,Int)
		self.layers 		= row.get('layers')     	# Int
		self.constrained	= row.get('constrained')    # Int
		self.symmetric 		= row.get('symmetric') 		# Bool
		self.vacuum 		= row.get('vacuum') 		# Int
		self.adsorbates 	= row.get('adsorbates') 	# {String : String} (adsorbate/interstitial : location)
		self.vacancies 		= row.get('vacancies') 		# [Int]
		self.emt          	= row.get('emt') 			# Float
		self.surfcomments 	= row.get('comments') 		# String
		self.relaxed 		= row.get('relaxed')
		self.kind 			= row.get('kind')

		self.xc       	= xc            # String
		self.psp      	= psp           # String
		self.pw       	= pw            # Int
		self.dwratio	= dwratio       # Int
		self.kptden   	= kptden        # Float

		if kpt is None: 	self.kpt = kptdensity2monkhorstpack(atoms,kptden)[0:2]+[1] #Force z kpt to be 1
		else:           	self.kpt = ast.literal_eval(kpt)         

		self.econv        = econv         # Float
		self.mixing       = mixing        # Float
		self.nmix         = nmix          # Int
		self.maxstep      = maxstep       # Float
		self.nbands       = nbands        # Int
		self.sigma        = sigma         # Float
		self.fmax         = fmax          # Float
		
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
			
	def showSites(self):
		with open('/scratch/users/ksb/img/sites.png','wb') as f:
			f.write(self.sites.decode('base64'))
		os.system('display /scratch/users/ksb/img/sites.png &')

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
		
		rand = [subAllocate(rx) for rx in [random.random(),random.random()]]+[sherlock]
		
		#return  [sherlock]*4 if self.dftcode=='gpaw' else rand 

		# Load the deck:
		return [suncat2,suncat,sherlock]

	def guessTime(self,jobkind): 
		mbeef = self.xc == 'mBEEF'
		if jobkind == 'relaxsurface':
			if mbeef: 	return 0.5
			else: 		return 0.5	
		elif jobkind == 'getsurfacexccontribs':
			if mbeef: 	return 0.1
			else: 		return 0.1
		elif jobkind == 'savesurfaceresults': 		return 0.1
		else: raise ValueError, 'unknown Job Kind '+jobkind
	
	def submit(self): 
		from fireworks   		import LaunchPad,Firework,Workflow
		from standardScripts 	import RelaxSurface,GetSurfaceXCcontribs,SaveSurfResults
		
		launchpad=LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')
		
		names = [x+'_%d'%self.ID for x in ['RelaxSurface','GetSurfaceXCcontribs','SaveSurfaceResults']]
		
		clusters = self.allocate() 
		timestamp = '_'+datetime.now().strftime('%Y_%m_%d_%H_%M')
		
		fw1 = Firework(	[RelaxSurface()]
						,name = names[0]
						,spec = {"jobID": 			self.ID
								,'_fworker': 		clusters[0].fworker
								,"_pass_job_info": 	True
								,"_files_out": 		{"fw1":"inp.gpw"}
								,"_queueadapter": 	clusters[0].qfunc(self.guessTime('relaxsurface'))
								,"_launch_dir": 	clusters[0].launchdir+names[0]+timestamp})
		
		fw2 = Firework(	[GetSurfaceXCcontribs()]
						,name 		= 	names[1]
						,parents 	= 	[fw1]
						,spec 		= 	{"jobID": 			self.ID
										,'_fworker': 		clusters[1].fworker
										,"_pass_job_info": 	True
										,"_files_in": 		{"fw1":"inp.gpw"}
										,"_queueadapter": 	clusters[1].qfunc(self.guessTime('getsurfacexccontribs'))
										,"_launch_dir": 	clusters[1].launchdir+names[1]+timestamp})
		
		fw3 = Firework(	[SaveSurfResults()]
						,name 		= 	names[2]
						,parents 	= 	[fw1,fw2]
						,spec 		= 	{"jobID":self.ID
										,'_fworker': 		clusters[2].fworker #MUST be sherlock
										,"_queueadapter": 	clusters[2].qfunc(self.guessTime('savesurfaceresults'))
										,"_launch_dir": 	clusters[2].launchdir+names[2]+timestamp})
		
		wflow=Workflow([fw1,fw2,fw3],name='SurfaceRelaxation_%d'%self.ID)
		
		print "Submitting job with ID = %d"%self.ID
		
		updateStatus('surfacejob',self.ID,'initialized','queued')
		
		launchpad.add_wf(wflow)




	########################################################################
	########################################################################
	# SQL METHODS
	########################################################################
	def sqlTable(self,tentative=False):  return 'tentativesurf' if tentative else 'surfacejob'

	def sqlCols(self):   return ['surfid'
								,'surfname','facet','xy','layers','constrained','symmetric','vacuum','adsorbates','vacancies','emt','surfcomments'
								,'xc','psp','pw','dwratio','kptden','kpt','econv','mixing','nmix','maxstep','nbands','sigma'
								,'fmax'
								,'precalcxc','dftcode'
								,'created','createdby','lastmodified','kind','comments'
								,'status']

	def sqlInsert(self): return [str(self.surfid)
								,self.surfname
								,str(self.facet)
								,str(self.xy)
								,self.layers
								,self.constrained
								,self.symmetric
								,self.vacuum
								,str(self.adsorbates)
								,str(self.vacancies)
								,str(self.emt)
								,self.surfcomments
								,self.xc
								,self.psp
								,self.pw
								,self.dwratio
								,self.kptden
								,str(self.kpt)
								,self.econv
								,self.mixing
								,self.nmix
								,str(self.maxstep)
								,str(self.nbands)
								,str(self.sigma)
								,str(self.fmax)
								,self.precalcxc
								,self.dftcode
								,self.created
								,self.createdby
								,self.lastmodified
								,self.kind
								,self.comments
								,self.status]

	def sqlEq(self): return ['surfid'
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
							,'precalcxc'
							,'dftcode']


#####################################
# OTHER FUNCTIONS
#####################################

def db2SurfaceJob(i,tentative=False): 
	table = 'tentativesurf' if tentative else 'surfacejob'
	row             = query1all(table,'id',i)
	immutableThings = (row[1],)+row[13:28]
	mutableThings   = (row[0],)+row[28:]
	return SurfaceJob(*(immutableThings+mutableThings))

###############################################################################
###############################################################################
###############################################################################

class SurfaceResult(object):
	def __init__(self,jobid,pos,magmom,e0,f0,xcContribs,avgtime,niter): 
		self.jobid 		= jobid
		self.pos 		= pos
		self.magmom 	= magmom
		self.e0 		= e0
		self.f0 		= f0
		self.xcContribs = xcContribs
		self.avgtime 	= avgtime
		self.niter 		= niter

	def __str__(self): return "SURFACERESULT FOR JOB#"+str(self.jobid)
	def addToASEdb(self,job):
		db    = connect('/scratch/users/ksb/db/ase.db')	
		
		info  = {'name':		job.surfname
				,'structure':	job.structure
				,'sites':		job.sites
				,'facet':		job.facet
				,'xy':			job.xy
				,'layers':		job.layers
				,'constrained':	job.constrained
				,'symmetric':	job.symmetric
				,'vacuum':		job.vacuum
				,'adsorbates':	job.adsorbates
				,'vacancies':	job.vacancies
				,'kind':		job.kind
				,'surfcomments':'%s\nInitial structure: %d'%(job.surfcomments,job.surfid)
				,'emt':job.emt
				,'relaxed':True
				,'jobid':job.ID}

		optAtoms = Atoms(symbols=job.symbols,scaled_positions=self.pos,cell=job.cell,magmoms=self.magmom,constraint=job.constraints,tags=job.tags)
		db.write(optAtoms,key_value_pairs=info)

		self.surfname     	= row.get('name')
		self.structure     	= row.get('structure')
		self.sites     		= row.get('sites') 			# String version of picture
		self.facet 		  	= row.get('facet')
		self.xy 			= row.get('xy')     		#(Int,Int)
		self.layers 		= row.get('layers')     	# Int
		self.constrained	= row.get('constrained')    # Int
		self.symmetric 		= row.get('symmetric') 		# Bool
		self.vacuum 		= row.get('vacuum') 		# Int
		self.adsorbates 	= row.get('adsorbates') 	# {String : String} (adsorbate/interstitial : location)
		self.vacancies 		= row.get('vacancies') 		# [Int]
		self.emt          	= row.get('emt') 			# Float
		self.surfcomments 	= row.get('comments') 		# String
		self.relaxed 		= row.get('relaxed')
		self.kind 			= row.get('kind')


	def getASEid(self):
		db    = connect('/scratch/users/ksb/db/ase.db')
		for row in db.select(relaxed=True,kind='surface'):
			if self.jobid == row.jobid: 
				return row.id
		raise ValueError, 'Tried to get ase id before it was even added to database?'

	def sqlTable(self):  return 'surfaceresult'

	def sqlCols(self):   return ['jobid'
								,'aseid'
								,'positions'
								,'magmoms'
								,'energy'
								,'forces'
								,'xccoeffs'
								,'time'
								,'niter']
	
	def sqlInsert(self): 
		return [self.jobid
				,self.getASEid()
				,str(self.pos)
				,str(self.magmom)
				,self.energy
				,str(self.forces)
				,str(self.xcCoeffs)
				,self.avgtime
				,self.niter]

	def sqlEq(self): return ['jobid']


