# External modules
import fireworks,subprocess,os,sys
#internal modules
import jobs,dbase,misc

lpad = fireworks.LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')

def readJobs():
	"""Looks for jobs, adds new ones to job table"""
	fwpathsher,fwpathsunc = '/scratch/users/ksb/fireworks/jobs/','/nfs/slac/g/suncatfs/ksb/fireworks/jobs/'
	existingJobs = [str(x[0]) for x in dbase.sqlexecute('SELECT launchdir from job')]
	ls = subprocess.Popen(['ssh','ksb@suncatls1.slac.stanford.edu', 'cd %s;ls'%fwpathsunc], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	suncout, err 	=  ls.communicate()

	suncatJobs 		= [fwpathsunc + d for d in suncout.split('\n') 		if fwpathsunc+d not in existingJobs and len(d)>1]
	sherlockJobs 	= [fwpathsher + x for x in os.listdir(fwpathsher) 	if fwpathsher+x not in existingJobs]
	tot = len(suncatJobs + sherlockJobs)
	for i,d in enumerate(suncatJobs + sherlockJobs):
		print d
		print '%d/%d'%(i+1,tot) ; sys.stdout.write("\033[F") # Cursor up one line
		fwid = getFWID(d)
		deleted = int(os.path.exists(d+'/deleted'))
		inputDict = misc.mergeDicts([{'fwid':fwid,'launchdir':d,'deleted':deleted},getInitData(fwid)])

		command = "INSERT into job ({0}) values ({1})".format(	','.join(inputDict.keys())
															,','.join(['?']*len(inputDict)))
		try: dbase.sqlexecute(command,inputDict.values())
		except: #remove 'bad keys'
			for k in ['relax','vacancies_json']:
				try: del inputDict[k]
				except KeyError: pass
			command = "INSERT into job ({0}) values ({1})".format(	','.join(inputDict.keys()),','.join(['?']*len(inputDict)))
			dbase.sqlexecute(command,inputDict.values())

def getFWID(jobpth): 	
	"""Only called when adding a new job as a row to the DB"""
	suffix = '/FW_submit.script'
	if '/nfs/' in jobpth:
		cat 	 = subprocess.Popen(['ssh','ksb@suncatls1.slac.stanford.edu', 'cat %s'%(jobpth)+suffix], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err =  cat.communicate()
		out 	 = out.split('\n')
	else:
		with open(jobpth+suffix,'r') as f: out = f.readlines()

	for l in out:
		if '--fw_id' in l: return int(l.split()[-1])
	raise ValueError, 'No fw_id found in FW_submit.script: \n\n%s'%out

def getInitData(fid):
	params 	= lpad.get_fw_dict_by_id(fid)['spec']['params']
	return jobs.assignJob(params).generalSummary()
