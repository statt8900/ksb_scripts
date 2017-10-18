import job
import os




xBulk1 = job.BulkStruct({'Li':1},{},'bcc', (1,1,1),{'Li':[0]},{})
xCalc1 = job.Calc('pbe',500,(5,5,5),0,5)
xConv1 = job.Convergence('sherlock',precalc,'00:20')

xBulk2 = job.BulkStruct({'Au':1,'Pd':1},{},'fcc', (1,1,2),{'Au':[0],'Pd':[1]},{})
xCalc2 = job.Calc('beef',400,(2,2,2),0,10)
xConv2 = job.Convergence('suncat2',precalc,'00:20')


xJob=job.Job('quantumEspresso',xBulk2,xCalc2,xConv2)


def scrape():
	""" [Job] """ 
	path   = os.getcwd() + '/'
	dirs   = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))] 
	for j in dirs:
		try:  
			print j
			print job.parseJob(j)
		except: pass


def update(jobList):
	for j 

def filterJob(filterDict,j):
	""" Map (Job -> a) a -> Job -> Bool """
	return all([f(j)==ans for f,ans in filterDict.items()])




testFilterDict  = {job.get_xc : job.pbe}

print filterJob(testFilterDict,xJob)

#scrape()