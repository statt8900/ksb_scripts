# External modules
import os,sys,shutil,json,pickle,fireworks,time,ase,itertools,ast,subprocess
import ase.visualize as viz
# Internal Modules
import details,misc,jobs,readJobs
import dbase as db
from printParse import *

############################
# MANAGE DATABASE
############################
lpad = fireworks.LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')
deletedpath = '/scratch/users/ksb/img/deletedfwids.txt'
def resetDB():
	if ask('Do you want to reset all tables?'):
		print 'wiping db'
		db.wipeDB();
		print 'loading details'
		details.load('status err_A err_BM')

def clearAllJobs():
	"""Deletes folders in $ALL_FWS, which should be repopulated w/in 2 minutes. Do this after having deleted something in $FW_PATH (on suncat/sherlock)"""
	pth = os.environ['ALL_FWS']
	shutil.rmtree(pth)
	os.mkdir(pth)
	for i in range(60): 
		print 'Waiting for alljobs to be populated'
		time.sleep(1)
	resetDB()

def deleteArchived():  delete(FWARCHIVED)

def delete(constraint,check=True,comment=None):
	output = db.query(['fwid','launchdir'],constraint)
	question = 'Are you sure you want to delete %d jobs? (selected by constraint: %s)'%(len(output),constraint)
	if (not check) or ask(question):
		for fwid,launchdir in output:
			db.updateDB('deleted','fwid',fwid,1,'job')
			db.updateDB('status','fwid',fwid,'deleted','job')
			if comment is not None: db.updateDB('comments','fwid',fwid,comment,'job')
			if 'scratch' in launchdir: 
				with open(launchdir+'/deleted','w') as f: f.write('')
			else:
				subprocess.Popen(['ssh','ksb@suncatls1.slac.stanford.edu', 'touch %s/deleted'%launchdir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			lpad.archive_wf(fwid)

def rerun(cnstrt='1'):
	fwids = db.queryCol('fwid',AND([FIZZLED,cnstrt]))
	question='Do you want to rerun %d  fizzled jobs selected by %s?'%(len(fwids),cnstrt)
	if ask(question):
		for f in fwids: lpad.rerun_fw(f)
		misc.launch()

def fix():
	output = db.query(['fwid','fwpckl'],AND([LATTICEOPT,GPAW]))
	for fw,fwpckl in output:
		q=pickle.loads(fwpckl)['spec']['_queueadapter']
		q['pre_rocket'] = q['pre_rocket'].replace(';;',';')
		lpad.update_spec([fw],{'_queueadapter':q})


def setWalltimeUnder40():
	""" If too much time was requested accidentally (such that job is stuck in READY), use this to modify the job"""
	details.load(['walltime'])
	output = db.query(['fwid','queueadapter'],QOVER40)
	print "Resetting walltime for %d jobs"%len(output)
	for fw,fwpckl in output:
		q=pickle.loads(fwpckl)['spec']['_queueadapter']
		q['walltime'] = '39:59' + ':00' if q['walltime'].count(':')>1 else ''
		lpad.update_spec([fw],{'_queueadapter':q})
	misc.launch()

def relaunch(cnst = '1',expiration_secs=80000,skip_load=False):
	"""Relaunch timed out jobs and unconverged GPAW jobs"""
	lpad.detect_lostruns(expiration_secs=expiration_secs, fizzle=True)
	lpad.detect_unreserved(expiration_secs=3600*24*7, rerun=True)

	if not skip_load: details.load('status fwpckl status trace',NOTCOMPLETED) 

	unconverged = AND([NOTCOMPLETED,KOHNSHAM,cnst])
	timedout 	= AND([TIMEOUT,NOTKOHN,cnst])
	launchflag	= False
	tOutput 	= db.query(['fwid','fwpckl'],timedout)
	uOutput 	= db.query(['fwid','params_json'],unconverged)
	tQuestion 	= "Do you want to relaunch %d timed out runs?"%len(tOutput)
	uQuestion 	= "Do you want to relaunch %d unconverged jobs?"%len(uOutput)

	if ask(tQuestion):
		launchflag = True
		for fid,fwpckl in tOutput:
			if lpad.get_fw_dict_by_id(fid)['state']=='FIZZLED':
				q 				= pickle.loads(fwpckl)['spec']['_queueadapter']
				wallT  			= q['walltime']
				q['walltime'] 	=  doubleTime(wallT)
				lpad.update_spec([fid],{'_queueadapter':q})
				lpad.rerun_fw(fid)
			else: print ("Wait up to 24 hours to relaunch fwid %d, or change default expiration time for detect_lostruns"%(fid))

	if ask(uQuestion):
		launchflag = True
		for fw,paramstr in uOutput:
			p = json.loads(paramstr)
			p['sigma'] 	+= 0.1 
			p['mixing'] = p['mixing']*0.5
			job= jobs.assignJob(p)
			if job.new():
				job.check();job.submit()
				delete(FWID(fw),check=False)
		readJobs.readJobs()

	if launchflag: 
		misc.launch()
		duplicates()


def duplicates(cnst='1',deleteFlag=False):
	output = db.query(['fwid','strjob'],cnst)
	rptDict={}
	for fwid,strjob in output:
		for f,s in output:
			if strjob == s and f!=fwid: 
				if fwid not in list(itertools.chain.from_iterable(rptDict.values())):
					if fwid in rptDict.keys(): rptDict[fwid].append(f)
					else: rptDict[fwid] = [f]
	print 'FWIDs with equal strjob entries: \n',abbreviateDict(rptDict)
	if deleteFlag:
		delfws = list(itertools.chain.from_iterable(rptDict.values()))
		if ask('Are you sure you want to delete %d duplicates?'%len(delfws)):
			for f in delfws: delete(FWID(f),False)


############
# QUERY DATA
############
def unconverged(var='raw_energy',xc='PBE',kptden=2,dftcode='gpaw'): raise NotImplementedError

def vizTraj(fwid): viz.view(pickle.loads(db.query1('inittraj_pckl','fwid',fwid)))

def viewImg(col='sites',constraint=SURFACE,fast=False):
	# view 'sites' of a SURFACE or 'bulkmodimg' of a BULKMOD
	imgs = db.queryCol(col,constraint)
	for i,img in enumerate(imgs):
		if fast or i==0 or raw_input('Continue?') or True:
			imgroot = os.environ['IMG_PATH']
			with open(imgroot+'/img.png','wb') as f:
				f.write(img.decode('base64'))
			os.system('display %s/img.png &'%imgroot)


##

def toDo():
	params = db.queryCol('params_json',VIB,deleted=True)
	for p in params:
		j = jobs.assignJob(json.loads(p))
		if j.new(): j.check();j.submit()
	misc.launch()
