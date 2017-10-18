from db import sqlexecute,dropTable

"""
Creates databases for bulk jobs and results
"""

def jobTable(tableName): return ("""CREATE TABLE %s       
					(jobid 			integer primary key
					,jobkind 		varchar
					,aseidinitial	integer
					,vibids			varchar
					,nebids			varchar
					,xc				varchar
					,pw				integer
					,kptden 		numeric
					,psp	 		varchar
					,xtol 			numeric
					,strain 		numeric
					,convid 		integer
					,precalc	 	varchar
					,dftcode 		varchar
					,comments 		varchar
					,error 			varchar
					,status 		varchar)"""%tableName) 

jobtable 			= jobTable('job')

convergenceTable = ("""CREATE TABLE convergence
					(convergenceid 	integer primary key
					,dwrat			numeric
					,econv			numeric
					,mixing			numeric
					,nmix			integer
					,maxstep		integer
					,nbands 		integer
					,sigma 			numeric
					,fmax 			numeric
					,delta 			numeric
					,climb 			varchar)""")

# include essential result information in table, just in case directory gets lost?
resultTable = ("""CREATE TABLE result 
					(resultid      	integer primary key
					,jobtableid	  	integer
					,launchdir 		varchar 
					,aseid 			integer
					,energy      	numeric
					,forces      	varchar
					,bulkmodulus 	numeric 
					,bfit        	numeric 
					,xccoeffs    	varchar
					,viblist    	varchar 
					,dos 			varchar
					,barrier 		varchar
					,time        	numeric
					,niter       	integer
					,foreign key(jobtableid) references job(jobid))""")

detailstable =  """CREATE TABLE details
								(detailid 		integer primary key
								,jobtableid 	integer
								,foreign key (jobtableid) references job(jobid))"""

def resetDetails(): 
	dropTable('details')
	sqlexecute(detailstable)

tentativedetailstable =  """CREATE TABLE tentativedetails
								(detailid 		integer primary key
								,jobtableid 	integer
								,foreign key (jobtableid) references tentativejob(jobid))"""

commands = [jobtable,convergenceTable,resultTable,detailstable]

#########################################################################
#########################################################################
#########################################################################
def main():
	for command in commands: 
		sqlexecute(command)

if __name__ == '__main__':
	main()
