import os,sys

from db 					import insertObject,deleteTable,plotQuery,queryManyFromOne
from ase.db 				import connect
from surfacejob 			import SurfaceJob,db2SurfaceJob
from constraint 			import *
from filters 				import *
from initializeCommonJobs 	import * #Domains
"""
Generate a list of jobs to run at some point in the future.

Inputs:
	- DOMAINS for all relevant DFT inputs (xcDomain,pspDomain,pwDomain, etc.)
	- FILTERS to narrow the combinatorial search space generated by the above Domains
	- CONSTRAINTS to more precisely conditions for adding a job 
		(can have complex dependencies between different parameters)

Output: Modifies bulkjob table of data.db
"""

asedb = connect('/scratch/users/ksb/db/ase.db')

"""
####################
# FILTERS TO APPLY #                        
####################                        
"""

filtersDuJour 	= [] 
filters 		= andFilters(essentialSurfJobFilters+filtersDuJour)

"""
########################
# CONSTRAINTS TO APPLY #
########################
"""
constraintsDuJour = [rpbe,qe] 
constraints = essentialSurfJobConstraints + constraintsDuJour

#########################################################
#########################################################

surfs      = [x for x in trajDomain if filters['aseid'](x)]

def main():

	deleteTable('tentativesurf') 

	for i,sID in enumerate(surfs):

		p 		= asedb.get(sID).get('parent')

		table 	= ' bulkjob INNER JOIN bulkresult ON bulkresult.jobid=bulkjob.id '
		cols 	= ['xc','psp','pw','kptden','kpt','econv','mixing','nmix','maxstep','nbands','sigma','fmax','precalcxc','dftcode']
		allInfo = queryManyFromOne(table,cols,'aseid',p)
		xc,psp,pw,kptden,kpt,econv,mixing,nmix,maxstep,nbands,sigma,fmax,precalcxc,dftcode = allInfo

		tentativeInput	 = [sID,xc,psp,pw,10,kptden,None,econv,mixing,nmix,maxstep,nbands,sigma,fmax,precalcxc,dftcode]
		tentativeSurfjob = SurfaceJob(*tentativeInput)
		
		ID,status 		 = insertObject(tentativeSurfjob,tentative=True)
		
	valid = [db2SurfaceJob(x[0],tentative=True) for x in plotQuery('tentativesurf',['tentativesurf.id'],constraints)] # SurfaceJob objects that satisfy constraints
	tot,new = len(valid),0
	question = '%d job(s) meet these criteria. Do you want to add them to surfacejobs in data.db?\n(y/n)--->'%tot
	if raw_input(question).lower() in ['y','yes']: 
		for z in range(len(valid)): 
			ID,status = insertObject(valid[z])
			if status == 'inserted': new +=1

	deleteTable('tentativesurf')
	print '%d jobs considered, %d are new'%(tot,new)
if __name__=='__main__':
	main()

