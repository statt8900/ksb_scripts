import fireworks,copy
import dbase, jobs,manageIncompleteJobs,manageSharedDatabase

print 'Getting old params'
lpad 		 = fireworks.LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')
allCompleted = lpad.get_fw_ids({'state':'COMPLETED'})
allParams    = [(fwid,lpad.get_fw_dict_by_id(fwid)['spec']['params']) for fwid in allCompleted]
newParams 	 = []
print 'Replacing inittraj_pckl with finaltraj_pckl'
for fwid,p in allParams:
	try: 
		ft = lpad.get_fw_dict_by_id(fwid)['launches'][-1]['action']['stored_data']['finaltraj_pckl'] 
		newp = copy.deepcopy(p)
		newp['inittraj_pckl']=ft
		if newp['jobkind']=='latticeopt': newParams.append(newp)
	except: pass

def sq(sqlcommand, binds = []):
	import sqlite3
	connection = sqlite3.connect('/scratch/users/ksb/db/latticeopt.db')
	cursor = connection.cursor()
	cursor.execute(sqlcommand, binds)
	if sqlcommand.lower()[:6] == 'select' or sqlcommand.lower()[:4] == 'with' or sqlcommand.lower()[:5]=='pragma':
		output_array = cursor.fetchall()
		connection.close()
		return output_array
	else:
		connection.commit()
		connection.close()
		return

if True: sq("CREATE TABLE jobs   (id             integer primary key"
									+',inittraj_pckl varchar'
									+',name varchar'
									+',kind varchar'
									+',comments varchar'
									+',structure varchar'
									+',dftcode varchar'
									+',bulkvacancy_json varchar'
									+',bulkscale_json varchar'
									+',pw numeric'
									+',xc varchar'
									+',kptden numeric'
									+',psp varchar'
									+',dwrat numeric'
									+',econv numeric'
									+',mixing numeric'
									+',nmix integer'
									+',maxstep integer'
									+',nbands integer'
									+',sigma numeric'
									+',fmax numeric'
									+',xtol numeric)')



for np in newParams:
	sq(('insert into jobs' 
	 	+' (inittraj_pckl,name,kind,comments,structure,dftcode,bulkvacancy_json,bulkscale_json,pw,xc,kptden,psp,dwrat,econv,mixing,nmix,maxstep,nbands,sigma,fmax,xtol)'
		+' values (?'+',?'*20+')'),[np['inittraj_pckl'],np['name'],np['kind'],np['comments'],np['structure'],np['dftcode'],np['bulkvacancy_json'],np['bulkscale_json'],np['pw'],np['xc'],np['kptden'],np['psp'],np['dwrat'],np['econv'],np['mixing'],np['nmix'],np['maxstep'],np['nbands'],np['sigma'],np['fmax'],np['xtol']])


"""

#latticeopts =  dbase.query(['fwid','jobkind','name','kind','comments','structure','dftcode','bulkvacancy_json','bulkscale_json','pw','xc','kptden','psp','dwrat','econv','mixing','nmix','maxstep','nbands','sigma','fmax','xtol'],"jobkind='latticeopt' and dftcode='gpaw'")
#params = ['inittraj_pckl','jobkind','name','kind','comments','structure','dftcode','bulkvacancy_json','bulkscale_json','pw','xc','kptden','psp','dwrat','econv','mixing','nmix','maxstep','nbands','sigma','fmax','xtol']

print "updating strjob"
manageSharedDatabase.updateDB('strjob')
print "loading list of jobstrs for incomplete jobs"
strs = manageIncompleteJobs.listOfIncompleteJobStrs()

def sub(p): jobs.Job(p).submit(strs)
count = 0
	if np['jobkind']=='latticeopt': count+=1 #sub(np)
	
	
print count"""