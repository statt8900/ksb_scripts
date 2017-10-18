import sqlite3,misc,json,printParse,sys,os


def callTest():
	os.system('python /scratch/users/ksb/test/test.py')

def randomPersonsOptPy():
	import datalog 	
	atoms = 'blah'
	print atoms
	datalog.log(atoms)


if __name__ == '__main__': 
	print 'hey'
	randomPersonsOptPy()
"""
Primary functions for interfacing with mai#n database, $SCRATCH/db/data.db 
"""

sherlockDB = '/scratch/users/ksb/db/data.db'

##########################################################
# QUERY FUNCTIONS
##########################################################

def query1(columnOut,testColumn,value,tableName='job'):
	command = "SELECT {1} FROM {0} WHERE {2} = ? ".format(tableName,columnOut,testColumn)
	output= sqlexecute(command,[value])
	if len(output) != 1: raise ValueError, 'query1 did not return one result '+str(output)
	return output[0][0]

def query(coloutputs,constraint='1',qtable='job',order='fwid',deleted=False):
	""" [STRING] -> CONSTRAINT -> [[sqloutput]] """

	select 	= 'select ' 	+ ','.join(coloutputs)
	table 	= ' from ' 		+ qtable
	delete 	= ((('' if deleted else ' not ')+' deleted and ') if qtable=='job' else '')
	where 	= ' where ' 	+ delete + constraint
	orderby = ' order by ' 	+ order
	command = select + table + where + orderby
	return sqlexecute(command)

def queryCol(col,constraint='1',qtable='job',order='fwid',deleted=False): return [x[0] for x in query([col],constraint,qtable,order,deleted)]

def anyQuery(constraint,qtable='job',deleted=False): return  len(queryCol('*',constraint,qtable, deleted=deleted)) > 0

def queryDistinct(coloutputs,groupbylist,constraint='1',order='fwid',qtable='job'):
	def ss(x): return "'"+x.replace("'","''")+"'" if isinstance(x,unicode) or isinstance(x,str) else str(x)

	distinct 	= sqlexecute('select DISTINCT {0} from {1}  where {2}'.format(','.join(groupbylist),qtable,constraint))
	colvals 	= [zip(groupbylist,x) for x in distinct]
	cnstrts 	= [' AND '.join([x+' = '+ss(y) for x,y in cv])  for cv in colvals]
	output 		= []
	for c in cnstrts: 
		output.append(query(coloutputs,constraint+' AND not deleted AND '+c,qtable,order))   
	return output
		

def queryTuple(n,constraint):
	"""
	Goal: generate a tuple of job ids so that in the future, not every combinatorial possibility has to be checked
	A custom WHERE statement filters all possibilities
	Properties of the first tuple element are accessed with j0.PROPERTY, etc.
	E.g. [all triples of jobs that have pw1<pw2<pw3] == queryTuple(3,['j0.pw < j1.pw','j1.pw < j2.pw'])

	These queries must only access unary Details, however any complex relationship can be PROJECTED as a unary detail
	E.g. "the 2nd element must be a molecule that has an adsorption energy on Cu with such-and-such calculator"
		can be projected as a unary property of the job that relaxed the molecule which was used in the adsorption job
	"""
	select 		= 'select ' +  ','.join(['j%d.fwid'%i for i in range(n)])
	table 		= ' from' + 'inner join'.join([' job as j%d '%i for i in range(n)])
	deleted 	= ' and '.join(['not j%d.deleted'%i for i in range(n)])
	where 		= ' where %s and '%deleted +constraint

	command = select + table + where
	#print command
	return sqlexecute(command)


############################
## MODIFY Functions
############################
def sqlexecute(sqlcommand, binds = []):

	connection = sqlite3.connect(sherlockDB,timeout=30)
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
def sqlexecutemany(sqlcommand, binds = []):
	#print sqlcommand
	connection 	= sqlite3.connect(sherlockDB,timeout=60)
	cursor 		= connection.cursor()
	cursor.executemany(sqlcommand, binds)
	if 'select' not in sqlcommand.lower(): connection.commit()
	cursor.close()
	connection.close()
	return

def updateDB(setColumn,idColumn,ID,newValue,tableName='job'):	
	sqlexecute("UPDATE {0} SET {1}= ? WHERE {2} = ?".format(tableName,setColumn,idColumn),[newValue,ID])
def dropTable(tableName): sqlexecute("DROP TABLE %s"%tableName)
def deleteRow(tableName,idColumn,ID): sqlexecute("DELETE FROM {0} WHERE {1} = {2}".format(tableName,idColumn,ID))
def addCol(colname,coltype,tableName): 	sqlexecute("ALTER TABLE {0} ADD COLUMN {1} {2}".format(tableName,colname,coltype))



	
################################################################################
def diffJob(fwid1,fwid2): 
	d1 = json.loads(query1('params_json','fwid',fwid1,'job'))
	d2 = json.loads(query1('params_json','fwid',fwid2,'job'))
	print printParse.abbreviateDict(misc.dict_diff(d1,d2))
	
################################################################################


def wipeDB():
	try: sqlexecute("DROP TABLE job")
	except: pass
	
	sqlexecute("""CREATE TABLE job 	(fwid 			integer not null
									,launchdir 		varchar not null
									,deleted 		integer not null
									,jobkind		varchar not null
									,inittraj_pckl	BLOB 	not null
									,comments		varchar 
									,trajcomments	varchar 
									,name 			varchar not null
									,relaxed		integer not null				
									,kind 			varchar not null
									,pw 			integer not null
									,xc 			varchar not null
									,kptden 		numeric not null
									,psp 			varchar not null
									,dwrat 			numeric not null
									,econv 			numeric not null
									,mixing 		numeric not null
									,nmix			integer not null		
									,maxstep 		integer not null
									,nbands 		integer not null
									,sigma 			numeric not null
									,fmax 			numeric not null
									,dftcode 		varchar not null
									
									,structure 			varchar
									,bulkvacancy_json 	varchar
									,bulkscale_json 	varchar
									,bulkparent 		integer
									,sites_base64 		BLOB						
									,facet_json 		varchar
									,xy_json 			varchar
									,layers 			integer
									,constrained 		integer
									,symmetric 			integer
									,vacuum 			integer

									,adsorbates_json 	varchar
									,vibids_json 		varchar
									,xtol 				numeric
									,strain 			numeric
									,delta 				numeric
									
									,params_json 		varchar not null
									,natoms				integer not null
									,symbols_pckl 		varchar not null
									,symbols_str 		varchar not null
									,species_pckl 		varchar not null
									,species_str 		varchar not null
									,comp_pckl 			varchar not null
									,comp_str 			varchar not null
									,numbers_pckl 		varchar not null
									,numbers_str 		varchar not null
									,metalstoich_pckl 	varchar not null
									,metalstoich_str 	varchar not null
									,metalspecies_pckl 	varchar not null
									,metalspecies_str 	varchar not null
									,metalcomp_pckl 	varchar not null
									,metalcomp_str	 	varchar not null
									,paramsinit_pckl 	varchar not null
									,kpt_pckl			varchar not null
									,kpt_str			varchar not null
									,emt				numeric not null
									,emtsym				integer not null
									,spinpol 			integer not null
									,dipole 			integer not null
									,calclabel			varchar not null
									,strjob				varchar not null

									,bulkvacancy_str	varchar
									,bulkscale_str 		varchar
									,facet_str 			varchar
									,xy_str 			varchar
									,surfacearea 		numeric

									,adsorbates_str 	varchar
									,dof 				integer
									,vibids_str 		varchar
									,nebfinal_pckl 		BLOB 

									,PRIMARY KEY (fwid,launchdir))""")





