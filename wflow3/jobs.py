#External Modules
import time,random,pickle
import datetime as d
# Internal Modules
import databaseFuncs,cluster


"""
Description:

Notes:
"Espresso executable doesn't seem to have been started"
	- caused by including /home/vossj/suncat/bin in PATH
"""

################################################
# Job base - must be inherited by all jobs
##############################################

class Job(object):
	"""Methods valid for all jobs"""

	def __init__(self,params):
		self.params  		= params
		self.jobkind 		= self.params['jobkind']		# STRING
		self.inittraj_pckl 	= self.params['inittraj_pckl'] 	# PICKLE
		self.name 			= self.params['name'] 			# STRING
		self.kind 			= self.params['kind'] 			# STRING (molecule/bulk/surface)

		self.pw 	= self.params['pw'] 		# INT
		self.xc 	= self.params['xc'] 		# STRING
		self.kptden = self.params['kptden'] 	# INT
		self.dftcode= self.params['dftcode'] 	# STRING



	def generalCheck(self):
		essential = set(['jobkind','inittraj_pckl','comments','name','kind'
						,'pw','xc','kptden','psp','dwrat','econv','mixing','nmix','maxstep','nbands','sigma','fmax','dftcode']) & set(self.kindEssentials()) & set(self.jobEssentials())
		actual 	= set(self.params.keys())
		assert essential <= actual, 'Missing the following general keys: '+str(essential - actual)

	def kindEssentials(self):
		if   self.kind == 'bulk': 	 return ['structure','bulkvacancy_json','bulkscale_json']
		elif self.kind == 'surface': return ['bulkparent','structure','sites_base64'
											,'facet_json','xy_json','layers','constrained','symmetric','vacuum','adsorbates_json']
		elif self.kind == 'molecule': return []
		else: raise NotImplementedError, 'unknown structure kind ',self.kind

	def jobEssentials(self):
		if   self.jobkind == 'latticeopt': 			return ['xtol']
		elif self.jobkind == 'bulkmod': 			return ['strain']
		elif self.jobkind == 'vib': 				return ['vibids_json','delta']
		elif self.jobkind in ['relax','vcrelax','xc']: return []
		else: raise NotImplementedError, 'unknown job kind '+self.jobkind

	def __str__(self):
		def wrapper(listx): return '_'.join(map(str,listx))
		def strTraj(self): return wrapper([getattr(pickle.loads(self.inittraj_pckl),x)() for x in ['get_atomic_numbers','get_positions','get_cell','get_initial_magnetic_moments']]).replace('\n','')
		def strCalc(self): return wrapper([self.params[x] for x in ['pw','xc','kptden','psp','dwrat','econv','mixing','nmix','maxstep','nbands','sigma','fmax','dftcode']])
		def strKind(self): return wrapper([self.params[x] for x in self.kindEssentials()])
		def strJob(self):  return wrapper([self.params[x] for x in self.jobEssentials()])
		return wrapper([ self.jobkind, strTraj(self), strCalc(self) ,strKind(self) , strJob(self)]).replace(' ','')

	def new(self,listOfIncompleteJobStrs=[]):
		if not databaseFuncs.anyQuery('strjob = \'%s\''%str(self)):
			if str(self) not in listOfIncompleteJobStrs: return True
		return False

	def submit(self,listOfIncompleteJobStrs=[]):
		""" use manageIncompleteJobs.listOfIncompleteJobStrs()"""
		from fireworks  import LaunchPad
		self.generalCheck()
		if self.new(listOfIncompleteJobStrs):
			launchpad 	= LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')
			wflow 		= self.wflow()
			launchpad.add_wf(wflow)
			time.sleep(2) # folder names are unique due to timestamp to the nearest second
			return 1
		else:
			print 'Repeat job!'
			return 0

	def allocate(self,n):
		if self.dftcode == 'gpaw':	return [cluster.sherlockgpaw]*n
		elif self.dftcode == 'quantumespresso':
			rand=[]
			def odds(randomx):
				if randomx < 0.3: 	return cluster.suncat
				elif randomx < 1: 	return cluster.suncat2
			for i in range(n): rand.append(random.random())
			return [odds(x) for x in rand]

	def guessTime(self):
		baseline = 5 # CHANGE THIS
		multFactor = self.pw/400.0 * self.kptden/2*(4 if self.xc=='mBEEF' else 1)
		multFactor *= 3 if 'hcp' in self.name else 1
		return min(40, 1 + baseline * multFactor)

	def standardScript(self):
		from standardScripts 	import OptimizeLattice,GetBulkModulus,VCRelax,GetXCcontribs,Vibrations,Relax,DOS,NEB
		scriptDict = {'latticeopt': OptimizeLattice(),'bulkmod': GetBulkModulus()
					,'vcrelax': VCRelax(),'relax': Relax(),'dos': DOS(),'neb': NEB()
					,'xc': GetXCcontribs(),'vib': Vibrations()}
		return scriptDict[self.jobkind]

	def wflow(self):
		from fireworks   		import Firework,Workflow

		cluster = self.allocate(1)[0]
		timestamp = '_'+d.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')

		fw1 = Firework([self.standardScript()]
						,name = self.jobkind
						,spec = {'params': 			self.params
								,'_fworker': 		cluster.fworker
								,'_queueadapter': 	cluster.qfunc(self.guessTime())
								,'_launch_dir': 	cluster.launchdir+self.jobkind+timestamp})
		return Workflow([fw1],name=self.jobkind)
