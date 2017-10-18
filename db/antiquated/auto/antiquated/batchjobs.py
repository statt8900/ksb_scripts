import os,json,pickle,shutil

from itertools    import product
from cluster      import getCluster,clusterDict
from job          import jobKind,parseBulkResult,parseSurfResult,parseBulkJob,parseSurfaceJob
from miniDatabase import clusterDomain,stoichStructDomain	
from printParse   import printList,printStoich,printAds

from plotting     import 	(xyPlot,xyResult,barResults,plotEq
							,calcLatParamA,exptLatParamA,errorLatParamA # Job   -> Float  ... xy function
							,calcBulkMod,exptBulkMod,errorBulkMod       # Job   -> Float  ... xy function
							,jobPW,jobKPT,jobTime						# Job   -> Float  ... xy function
							,getStoich,getPW,getKPT						# Job   -> String ... label function
							,maeLatParamA,maeBulkMod)                   # [Job] -> Float  ... bar function

from filters import * #really we just want all the filters
from prettytable import PrettyTable
from copy import deepcopy

"""
This is the main script that a user interacts with to manage their workflow
Filters from the filter module are used to identify a subset of jobs to initialize or plot
The scrape function updates a pickle file that places each job (folder in root directory)
into a status bin. 
Other functions do stuff to jobs depending on what status they have.
Plotting relies on the PlotFunc objects, defined in plotting module
"""

######################################
# THINGS TO DO WITH LISTS OF JOBS
######################################
def initializeJobs(jobs):	
	writecluster = getCluster(2).name
	print "Initializing %d jobs: "%(len(jobs))
	for j in jobs:
		print 'Initializing '+str(j)
		for c in clusterDomain: j.initialize(writecluster,c)
	
########################################################
########################################################
def getStatus(index):
	resultRoot = getCluster(2).result # folder to store result pickle file
	with open(resultRoot+'/status.pckl','rb') as f: statusDict = pickle.load(f)
	return statusDict[index]

def clearWF():
	cluster = getCluster(2)
	root    = cluster.root   # folder with jobs as subfolders
	jobs    = os.listdir(root)
	inProgress = map(str,getStatus('inProgress'))
	for j in jobs:
		if j not in inProgress: 
			files = os.listdir(root+'/'+j)
			for f in files:
				if 'gpw' in f: os.remove(f)

def scrape():
	"""
	Scrapes directory, returning dictionary { keywords : [Job] }
	Alternatively, takes a keyword ('complete','timeout') and returns [Job] directly
	"""
	cluster = getCluster(2)
	root       = cluster.root   # folder with jobs as subfolders
	resultRoot = cluster.result # folder to store result pickle file

	statuses = ['complete'
				,'incomplete'
					,'initialized'	
					,'cancelled'
					,'inProgress'
					,'scfFailed'
					,'unconverged'
					,'unknownFunctional'
					,'diagonalize'
					,'pspStates'
					,'typeError','valueError' #... other generic errors as catch all?
					,'timeout']

	jobs   = os.listdir(root)

	status = {}
	for s in statuses: status[s] = [] #initialize

	def parse(x): return parseSurfaceJob(x) if ('Surf' in x) else parseBulkJob(x)

	def updateStatus(statusDict,job,clusterName,err):
		if   cluster.timeout(err):       statusDict['timeout'].append(job)
		elif cluster.unconverged(err):   statusDict['unconverged'].append(job)
		elif cluster.fail(err):          statusDict['scfFailed'].append(job)
		elif 'CANCELLED' in err:         statusDict['cancelled'].append(job)
		elif 'diagonalize error' in err: statusDict['diagonalize'].append(job)
		elif 'Cannot figure out what states should exist on this pseudopotential' in err: statusDict['pspStates'].append(job)
		elif 'Unknown functional' in err: statusDict['unknownFunctional'].append(job)
		elif 'TypeError' in err: statusDict['typeError'].append(job)
		elif 'ValueError' in err: statusDict['valueError'].append(job)
		else: statusDict['inProgress'].append(job)
		return statusDict

	for j in jobs:
		pth = root+'/'+j
		if os.path.isfile(pth+'/result.json'): 	status['complete'].append(parse(j))
		else:
			status['incomplete'].append(parse(j))
			if os.path.isfile(pth+'/myjob.err'): 
				with open(pth+'/myjob.err','r') as f: status = updateStatus(status,parse(j),'sherlock',f.read())
			else: status['initialized'].append(parse(j))
	
	with open(resultRoot+'/status.pckl','wb') as handle: pickle.dump(status, handle,protocol=pickle.HIGHEST_PROTOCOL)
	print "Done scraping "+root
	x = PrettyTable([s for s in statuses if len(status[s]) > 0])
	x.add_row([len(status[z]) for z in statuses if len(status[z]) > 0])
	print (x)

##########################################################################
##########################################################################
def modifyStatus(status,function):
	js = getStatus(status)
	print "Doing something to %d jobs with status %s"%(len(js),status)
	for j in js: function(j)

def printStatus(status): 	
	js = getStatus(status)
	print "%d jobs with status %s"%(len(js),status)
	for j in js: print j

# Functions to modify with
def deleteFolder(job):	shutil.rmtree(job.jobPth(getCluster(2).name))

def initJob(job): 
	print 'Initializing '+str(job)
	for c in clusterDomain: job.initialize(getCluster(2).name,c)

def submitJob(job):
	cluster = getCluster(2)
	os.chdir(cluster.root+'/'+str(job))
	print "In directory "+os.getcwd()
	if cluster.name=='sherlock': 
		print "Submitting "+str(job)+" from sherlock"
		os.system('./sherlock_metarun.sh')
	else: 
		print "Submitting "+str(job)+" from suncat"
		os.system('./suncat2_metarun.sh')

def testRedundancy(job): 
	raise NotImplementedError
	# Search through all jobs, filtering by stoichiometry (necessary condition)
	# compare EMT calculation of initial structures: if same to 0.0000001 or something, call redundant
	# Worry: energetically equivalent but structurally-unique? 

#####################
# CONVERGENCE TACTICS
#####################
def decreaseMixing(n):
	def decreaseFunc(job):
		newJob = deepcopy(job)
		newJob.calc.mixing-=n
		print "Initializing new job with lower mixing: ",newJob
		initJob(newJob)
		print "Deleting old job: ",job
		deleteFolder(job)
	return decreaseFunc

def increaseSigma(n):
	def increaseFunc(job):
		newJob = deepcopy(job)
		newJob.calc.sigma+=n
		print "Initializing new job with higher sigma: ",newJob
		initJob(newJob)
		print "Deleting old job: ",job
		deleteFolder(job)
	return increaseFunc

def increaseEconv(n):
	def increaseFunc(job):
		newJob = deepcopy(job)
		newJob.calc.eConv+=n
		print "Initializing new job with higher eConv: ",newJob
		initJob(newJob)
		print "Deleting old job: ",job
		deleteFolder(job)
	return increaseFun

#################################
# On the fly filters

#[li500MBEEF,li700MBEEF,li900MBEEF,li1100MBEEF,li1300MBEEF,li1500MBEEF,li1700MBEEF]
#[li555MBEEF,li111111MBEEF,li151515MBEEF]
###############################################
# COMMANDS
###############################################
if __name__=='__main__':
	scrape()
	#printStatus('unconverged')
	#jobs = filterBulk([johannesBulkFilter,nonmag,mbeef]) # + filterSurf([testSuiteBulk],[defaultSurf]) #
	#initializeJobs(jobs)
	#scrape()

	#modifyStatus('initialized',decreaseMixing(5))
	#modifyStatus('initialized',submitJob)
	if False:
		import matplotlib.pyplot as plt
		
		#fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2)
		fig, (ax1) = plt.subplots(nrows=1,ncols=1)
		#fig.subplots_adjust(hspace=0.6,wspace=0.4)

		xyPlot(jobPW
			,errorBulkMod
			,[li555MBEEF,li111111MBEEF,li151515MBEEF]
			,ax1
			,['r--','g--','b--','m--','y--','k--','c--'])


		#x,y = xyResults(jobPW,calcLatParamA,LiBccKPT151515,ax1)
		#x,y = xyResults(exptLatParamA,calcLatParamA,pureHighAcc,ax2,getStoich)
		#plotEq(2,6,ax2)
		#barResults(maeLatParamA,[pureHCP,pureFCC,pureBCC],ax3)
		#barResults(maeBulkMod,[pureHCP,pureFCC,pureBCC],ax4)
		plt.show()

	#printStatus('diagonalize')

