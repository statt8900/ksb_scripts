import os,time,subprocess,shutil

from db 			import getJobIDs,allJobIDs
from db 			import countDB,updateDB,updateStatus,deleteRow,insertObject
from db 			import query1,plotQuery
from constraint     import *
from fireworks   	import LaunchPad
from job 			import db2object
"""
Contains useful functions for interacting with ASE, data.db (Sqlite), and FireWorks (MongoDB) databases

Desired Functionality:
	"DO THIS all jobs that meet THESE CONDITIONS"
	- where THIS could be:
		- submit to clusters (Fireworks)
		- create+launch new jobs with new convergence parameters
		- open up ASE GUI to view traj file
		- ? ? ? 
	- and where THESE CONDITIONS could be
		- ... is a Lithium surface job that has not been run yet
		- ... fizzled with an error KohnShamConvergenceError
		- ... completed (but are flagged for adsorbate changing site from initial site)
		- ? ? ? 

"""

archiveErrors = ['KohnShamConvergenceError' #errors that do not require individual attention
				,'TIMEOUT']                 #wflows can be archived, other methods will scan for certain types of errors and resubmit

############################
# Things to do with job sets
############################
lpad = LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')

class lpadIDjob(object):
	"""Kind of superfluous, but a way to store all info about a job in one place """
	def __init__(self,lpadid):
		fw 			= lpad.get_fw_by_id(lpadid)
		jobid 		= fw.spec['jobID']
		self.id 	= lpadid
		self.fw 	= fw
		self.jobid 	= jobid
		self.job 	= db2object(jobid)

	def err(self): 
		action =  self.fw.launches[0].action
		return None if action is None else action.stored_data['_exception']['_stacktrace'] 

# Get all the fws associated with a job
def fromJobID(jid): return [lpadIDjob(l) for l in lpad.get_fw_ids({'spec.jobID':jid})]
################################
#Things to do with lists of fws
################################

def archiveFizzled(): 
	updateFizzled()
	fizzled = [lpadIDjob(x) for x in getStatus("FIZZLED")]
	for f in fizzled: 
		err = f.err() 
		if err is None: 
			updateDB('error','jobid',f.jobid,'TIMEOUT',None,'job')
			updateStatus(f.jobid,'fizzled','archived')
			lpad.archive_wf(f.id)

		elif any([x in err for x in archiveErrors]):
			updateDB('error','jobid',f.jobid,err,None,'job')
			updateStatus(f.jobid,'fizzled','archived')
			lpad.archive_wf(f.id)
		
def increaseConvergence():
	kohnsham = [x[0] for x in plotQuery(['jobid'],[KOHNSHAM])]
	new  = 0
	for k in kohnsham:
		jobj = db2object(k)
		jobj.convid += 1
		jobj.status='initialized'
		jobj.error,jobj.comments=None,None
		id,status=insertObject(jobj)
		if status == 'inserted': new+=1
	print 'added %d new jobs'%new

########
# ERRORS
########
def printErr(lastLines=True):
	for i in getStatus("FIZZLED"):
		try:
			l=lpadIDjob(i)

			try: err=l.err()
			except AttributeError: err='\nNone\n'
			if lastLines: err = err.split('\n')[-2]
			print l.id,err
		except AssertionError: pass
"""
def kohnSham(): 	return [l for l in [lpadIDjob(x) for x in getStatus("FIZZLED")] if "KohnShamConvergenceError" in l.err()]
def sysExit(): 		return [l for l in [lpadIDjob(x) for x in getStatus("FIZZLED")] if "SystemExit: 1" in l.err()]
def pspStates(): 	return [l for l in [lpadIDjob(x) for x in getStatus("FIZZLED")] if "ValueError: Cannot figure out what states should exist on this pseudopotential" in l.err()]
"""

##############################
#Querying to get lists of jobs
##############################

def getStatus(status): return lpad.get_fw_ids({'state':status})
	
def updateFizzled():
	for fw in [lpadIDjob(x) for x in getStatus('FIZZLED')]:
		if query1('status','jobid',fw.jobid,'job')!='fizzled':
			updateStatus(fw.jobid,'queued','fizzled')

def nResults(): print countDB('bulkresult')

def jobIDs(table,x):
	if isinstance(x,int): return [x]    				# input is jobId
	elif isinstance(x,list): 
		if len(x)==0 or isinstance(x[0],int): return x 	# input is [jobID]
		else: return getJobIDs(table,x) 						# input is [Constraint]

def daemon():
	updateFizzled()
	archiveFizzled()
	increaseConvergence()
#############################
def submitAll(): 
	initializedJobs = [x[0] for x in plotQuery(['jobid'],[INITIALIZED])]
	if raw_input('Submit %d jobs?\n(y/n)--> '%len(initializedJobs)).lower() in ['y','yes']:
		for jid in initializedJobs:
			db2object(jid).submit()


#############################
#############################
def main():
	pass # not really necessary to have a main? 
	
if __name__ == '__main__':
	main()
