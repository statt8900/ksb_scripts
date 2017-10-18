import os,time,subprocess

from db 			import getJobIDs,countDB,updateStatus,updateDB,query1,plotQuery,deleteRow
from bulkjob 		import db2BulkJob
from constraint     import *
from fireworks   	import LaunchPad
from fireworks.queue.queue_adapter 	import QueueAdapterBase
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter
from fireworks.core.fworker import FWorker

"""
Contains useful functions for interacting with ASE, data.db (Sqlite), and FW (MongoDB) databases

"""

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

	def err(self): return self.fw.launches[0].action.stored_data['_exception']['_stacktrace']


def fromJobID(jid): return [lpadIDjob(l) for l in lpad.get_fw_ids({'spec.jobID':jid})]
################################
#Things to do with lists of fws
################################
def removeall(idlist):
	for i in idlist: remove(lpadidjob)
def remove(lpadidjob): #arg either int or lpadIDjob
	try:
		lpad.delete_wf(lpadidjob.id)
	except TypeError:
		print "could not remove wf %d from job %d due to typeerror"%(lpadidjob.id,lpadidjob.jobid)

def removeFizzled():
	updateFizzled()
	fizzled = [lpadIDjob(x) for x in getStatus("FIZZLED")]
	for f in fizzled: 
		updateStatus('bulkjob',f.jobid,'fizzled','initialized')
		remove(f.id)

def printAll(jids):
	for j in jids:
		print j.id
		print j.fw

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
def kohnSham(): 	return [l for l in [lpadIDjob(x) for x in getStatus("FIZZLED")] if "KohnShamConvergenceError" in l.err()]
def sysExit(): 		return [l for l in [lpadIDjob(x) for x in getStatus("FIZZLED")] if "SystemExit: 1" in l.err()]
def pspStates(): 	return [l for l in [lpadIDjob(x) for x in getStatus("FIZZLED")] if "ValueError: Cannot figure out what states should exist on this pseudopotential" in l.err()]

def updateFizzled():

	for l in getStatus("FIZZLED"):
		try:
			f=lpadIDjob(l)
			try: err=f.err()
			except AttributeError: err='None'

			if query1('bulkjob','status','id',f.jobid) != 'fizzled':
				updateStatus('bulkjob',f.jobid,'queued','fizzled') #unsafe: updateDB('bulkjob','status','id',f.jobid,'fizzled')
		except AssertionError: pass
		updateDB('bulkjob','comments','id',f.jobid,err)

##############################
#Querying to get lists of jobs
##############################
"""
This is all bad code. Write different functions rather than overloading
"""
def allFWs(): return [lpadIDjob(x) for x in lpad.get_fw_ids()]
def getStatus(status): return lpad.get_fw_ids({'state':status})
def incomplete(): return [lpadIDjob(x) for x in lpad.get_fw_ids() if x not in getStatus('COMPLETED')]
def nResults(): print countDB('bulkresult')

def jobIDs(table,x):
	if isinstance(x,int): return [x]    				# input is jobId
	elif isinstance(x,list): 
		if len(x)==0 or isinstance(x[0],int): return x 	# input is [jobID]
		else: return getJobIDs(table,x) 						# input is [Constraint]


def incompleteDB(table): return [x[0] for x in plotQuery(table,['bulkjob.id'],[incompleted])]

def removeIncomplete(table):
	for i  in incompleteDB(table): 
		lpadidjobs = fromJobID(i)
		for l in lpadidjobs: 
			remove(l)
		deleteRow('bulkjob','id',i)

#############################
def addToLpad(table,x): 
	for jID in jobIDs(table,x): 
		db2BulkJob(jID).submit()

def addSurfToLpad(table,x): 
	for jID in jobIDs(table,x): 
		db2SurfaceJob(jID).submit()


def addInitialized(table):
	addToLpad(table,getJobIDs(table,[initialized]))

def addFunc(table): 
	if 'surf' in table: 
		addToLpad(table,getJobIDs(table,testFilters[table]))
	elif 'bulk' in table:
		addToSurfLpad(table,getJobIDs(table,testFilters[table]))

#############################
#############################
def main():
	pass # not really necessary to have a main? 
	
if __name__ == '__main__':
	main()
