from db import sqlexecute

"""
Creates databases for bulk jobs and results
"""


bulkjob = ("""CREATE TABLE bulkjob       
					(id 			integer primary key
					,bulkid			integer
					,symbols		varchar
					,cell			varchar
					,a 				numeric
					,b 				numeric
					,c				numeric
					,alpha			numeric
					,beta			numeric
					,gamma			numeric
					,positions		varchar
					,magmoms		varchar
					,tags 			varchar
					,constraints	varchar
					,bulkname 		varchar
					,structure 		varchar
					,emt			numeric
					,bulkcomments	varchar
					,xc 			varchar
					,psp 			varchar
					,pw				integer
					,kptden 		integer
					,kpt 			varchar 
					,econv			numeric
					,mixing			numeric
					,nmix			integer
					,maxstep		integer
					,nbands 		integer
					,sigma 			numeric
					,fmax 			numeric
					,xtol 			numeric
					,strain 		numeric
					,precalcxc	 	varchar
					,dftcode 		varchar
					,created 		numeric
					,createdby	 	varchar
					,lastmodified 	numeric
					,kind 			varchar
					,comments 		varchar
					,status 		varchar)""") 


surfacejob = ("""CREATE TABLE surfacejob       
					(id 			integer primary key
					,surfid			integer
					,surfname 		varchar
					,facet	 		varchar 
					,xy 		 	varchar 
					,layers	 		integer 
					,constrained	integer
					,symmetric	 	varchar
					,vacuum	 		integer
					,adsorbates	 	varchar 
					,vacancies	 	varchar 
					,emt			numeric
					,surfcomments	varchar
					,xc 			varchar
					,psp 			varchar
					,pw				integer
					,dwratio		integer
					,kptden 		integer
					,kpt 			varchar 
					,econv			numeric
					,mixing			numeric
					,nmix			integer
					,maxstep		integer
					,nbands 		integer
					,sigma 			numeric
					,fmax 			numeric
					,precalcxc	 	varchar
					,dftcode 		varchar
					,created 		numeric
					,createdby	 	varchar
					,lastmodified 	numeric
					,kind 			varchar
					,comments 		varchar
					,status 		varchar)""") 


bulkresult = ("""CREATE TABLE bulkresult 
					(id      	    integer primary key
					,jobid	  		integer
					,aseid 			integer
					,cell			varchar
					,a 				numeric
					,b 				numeric
					,c				numeric
					,alpha			numeric
					,beta			numeric
					,gamma			numeric
					,positions		varchar
					,magmoms		varchar
					,energy      	numeric
					,forces      	varchar
					,bulkmodulus 	numeric
					,bfit        	numeric
					,xccoeffs    	varchar
					,time        	numeric
					,niter       	integer
					,foreign key(jobid) references bulkjob(id))""")

surfaceresult = ("""CREATE TABLE surfaceresult 
					(id      	    integer primary key
					,surfjobid	  	integer
					,aseid 			integer
					,positions		varchar
					,magmoms		varchar
					,energy      	numeric
					,forces      	varchar
					,xccoeffs    	varchar
					,time        	numeric
					,niter       	integer
					,foreign key(surfjobid) references surfjob(id))""")

commands = [bulkjob,surfacejob,bulkresult,surfaceresult]

#########################################################################
#########################################################################
#########################################################################
def main():
	for i,command in enumerate(commands): 
		sqlexecute(command)

if __name__ == '__main__':
	main()
