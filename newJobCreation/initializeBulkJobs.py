import manageIncompleteJobs as m
import jobs

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

#####
output = sq("select * from jobs where (xc = 'BEEF' or xc = 'mBEEF')") # CHANGE THIS CONSTRAINT TO GET PBE JOBS IN THE FUTURE?
n = len(output)
for i,(ID,pckl,name,kind,comments,structure,dftcode,bulkvac,bulkscale,pw,xc,kptden,psp,dwrat,econv,mixing,nmix,maxstep,nbands,sigma,fmax,xtol) in enumerate(output):
	print '%d/%d'%(i,n)
	p = {'inittraj_pckl':pckl
	,'name':name
	,'kind':kind
	,'jobkind':'latticeopt'
	,'comments':comments
	,'structure':structure
	,'dftcode':dftcode
	,'bulkvacancy_json':bulkvac
	,'bulkscale_json':bulkscale
	,'pw':pw
	,'xc':xc
	,'kptden':kptden
	,'psp':psp
	,'dwrat':dwrat
	,'econv':econv
	,'mixing':mixing
	,'nmix':nmix
	,'maxstep':maxstep
	,'nbands':nbands
	,'sigma':sigma
	,'fmax':fmax
	,'xtol':xtol}
	jobs.Job(p).generalCheck()
	existing = m.listOfIncompleteJobStrs()

	jobs.Job(p).submit(existing)



