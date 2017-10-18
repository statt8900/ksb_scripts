import os
from db import updateStatus,query1,updateDB

root = '/scratch/users/ksb/db/jobs'

for d in os.listdir(root):
	try:
		print 'trying ',d
		with open('/scratch/users/ksb/db/jobs/'+d+'/myjob.err','r') as f:
			err = f.read()
			if 'TIME LIMIT' in err: 
				if query1('bulkjob','status','id',d) is not 'failed':
					updateStatus(int(d),'running','timeout')
			elif 'CANCELLED AT' in err:
				updateDB('bulkjob','status','id',int(d),'cancelled')

	except IOError: pass
			