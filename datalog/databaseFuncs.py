import sqlite3,json,sys,os
import misc
"""
Primary functions for interfacing with mai#n database, $SCRATCH/db/data.db 
"""

sherlockDB = '/scratch/users/ksb/share/suncatdata.db'
jointTable = 'completed left join derived on jobid = derived_job left join keyword on jobid = keyword_job left join calc on jobcalc = calc_id'

##########################################################
# QUERY FUNCTIONS
##########################################################

def query1(columnOut,testColumn,value,tableName=jointTable):
	command = "SELECT {1} FROM {0} WHERE {2} = ? ".format(tableName,columnOut,testColumn)
	output= sqlexecute(command,[value])
	if len(output) != 1: raise ValueError, 'query1 did not return one result '+str(output)
	return output[0][0]

def query(coloutputs,constraint='1',qtable=jointTable,order='jobid',limit=None,deleted=False):
	""" [STRING] -> CONSTRAINT -> [[sqloutput]] """

	select 	= 'select ' 	+ ','.join(coloutputs)
	table 	= ' from ' 		+ qtable
	delete 	= (('' if deleted else ' not ')+' deleted and ') if 'completed'  in qtable else ''
	where 	= ' where ' 	+ delete + constraint
	orderby = ' order by ' 	+ order
	limit   = '' if limit is None else ' limit %s '%str(limit)
	command = select + table + where + orderby + limit
	#print 'command = ',command
	return sqlexecute(command)

def queryCol(col,constraint='1',qtable=jointTable,order='jobid',deleted=False): return [x[0] for x in query([col],constraint=constraint,qtable=qtable,order=order,deleted=deleted)]

def anyQuery(constraint,qtable=jointTable,deleted=False): return  len(queryCol('*',constraint,qtable, deleted=deleted)) > 0

def queryDistinct(coloutputs,groupbylist,constraint='1',order='jobid',qtable=jointTable):
	def ss(x): return "'"+x.replace("'","''")+"'" if isinstance(x,unicode) or isinstance(x,str) else str(x)
	print 'entering queryDistinct with coloutputs',coloutputs
	print 'groupbylist ',groupbylist
	print 'constraint ',constraint
	distinct 	= sqlexecute('select DISTINCT {0} from {1}  where {2}'.format(','.join(groupbylist),qtable,constraint))

	print 'distinct',distinct
	colvals 	= [zip(groupbylist,x) for x in distinct]
	cnstrts 	= [' AND '.join([x+' = '+ss(y) for x,y in cv])  for cv in colvals]
	output 		= []

	print 'constraints ',cnstrts

	for c in cnstrts: 
		
		output.append(query(coloutputs,constraint+' AND not deleted AND '+c,qtable,order))   

	print 'output to queryDistinct ',output
	return output
		
"""
def queryTuple(n,constraint):
	Goal: generate a tuple of job ids so that in the future, not every combinatorial possibility has to be checked
	A custom WHERE statement filters all possibilities
	Properties of the first tuple element are accessed with j0.PROPERTY, etc.
	E.g. [all triples of jobs that have pw1<pw2<pw3] == queryTuple(3,['j0.pw < j1.pw','j1.pw < j2.pw'])

	These queries must only access unary Details, however any complex relationship can be PROJECTED as a unary detail
	E.g. "the 2nd element must be a molecule that has an adsorption energy on Cu with such-and-such calculator"
		can be projected as a unary property of the job that relaxed the molecule which was used in the adsorption job
	select 		= 'select ' +  ','.join(['j%d.fwid'%i for i in range(n)])
	table 		= ' from' + 'inner join'.join([' job as j%d '%i for i in range(n)])
	deleted 	= ' and '.join(['not j%d.deleted'%i for i in range(n)])
	where 		= ' where %s and '%deleted +constraint

	command = select + table + where
	#print command
	return sqlexecute(command)
	"""


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

def updateDB(setColumn,idColumn,ID,newValue,tableName='derived'):	
	sqlexecute("UPDATE {0} SET {1}= ? WHERE {2} = ?".format(tableName,setColumn,idColumn),[newValue,ID])
def dropTable(tableName): 
		try: sqlexecute("DROP TABLE %s"%tableName)
		except Exception as e: print 'dropTable exception: ',e
def deleteRow(tableName,idColumn,ID): sqlexecute("DELETE FROM {0} WHERE {1} = {2}".format(tableName,idColumn,ID))
def addCol(colname,coltype,tableName): 	sqlexecute("ALTER TABLE {0} ADD COLUMN {1} {2}".format(tableName,colname,coltype))
def diffFWID(f1,f2,inputOnly = True):
	select = 'fwid,job_type_ksb,job_name,dftcode,xc,planewave_cutoff,kptden_ksb,structure_ksb,bulkvacancy_ksb,bulkscale_ksb,psp_ksb,dwrat_ksb,econv_ksb,mixing_ksb,maxstep_ksb,nmix_ksb,nbands_ksb,sigma_ksb,fmax_ksb,xtol_ksb,strain_ksb,delta_ksb,system_type_ksb,strjob ' if inputOnly else '*'
	rows = []
	def dict_from_row(row): return dict(zip(row.keys(), row))       
	for fw in [f1,f2]:
		connection 	= sqlite3.connect(sherlockDB)
		connection.row_factory = sqlite3.Row
		
		cursor = connection.execute('select %s from %s where fwid=%d'%(select,jointTable,fw))
		rows.append(dict_from_row(cursor.fetchone()))
		connection.close()


	for item in misc.dict_diff(rows[0],rows[1]).items(): print item

def createSuncatDataDB():
    #Initializes sqlite database with two tables. Run only once.
    #    'completed' - static information logged during runtime upon job completion
    #    'derived'   - dynamic information derivable from static information, meant to be expanded/modified

    sqlexecute("CREATE TABLE completed    (jobid             integer primary key"
                                        +',deleted           integer not null'
                                        +',user              varchar not null'
                                        +',timestamp         numeric not null'
                                        +',working_directory varchar not null unique'
                                        +',storage_directory varchar not null unique'
                                        +',jobcalc           integer not null '
										+',FOREIGN KEY(jobcalc) REFERENCES calc (calc_id))')

    sqlexecute("CREATE TABLE derived      (derived_id integer primary key"
                                        +',derived_job integer not null'
                                        +',FOREIGN KEY(derived_job) REFERENCES completed (jobid))')

    sqlexecute("CREATE TABLE refeng      (refeng_id integer primary key"
                                        +',element     integer not null '
                                        +',refeng_job  integer not null unique'
                                        +',refeng_calc integer not null'
                                        +',refeng      numeric not null'
                                        +',FOREIGN KEY(refeng_job) REFERENCES completed (jobid)'
                                        +',FOREIGN KEY(refeng_calc) REFERENCES calc (calc))')

    sqlexecute("CREATE TABLE calc         (calc_id integer primary key"
										+',dftcode 	varchar not null'
                                     	+',xc 		 varchar'
                                     	+',pw 		 numeric'
										+',kpts_json varchar'
										+',fmax      numeric '
										+',psp      varchar'
										+',econv    numeric '
										+',xtol     numeric '
										+',strain   numeric '
										+',dw 		numeric'
										+',sigma 	numeric'
										+',nbands 	integer'
										+',mixing 	numeric'
										+',diag 	varchar'
										+',spinpol 	integer'
										+',dipole 	integer'
										+',maxstep 	integer)')

    sqlexecute("CREATE TABLE keyword      (keyword_id integer primary key"
                                      	+',keyword_job integer not null'
                                     	+',FOREIGN KEY(keyword_job) REFERENCES completed (jobid))')
    return 1

def table(colstr='fwid',cnstr='1',order='jobid'):
    import sqlite3,prettytable
    colstr      = colstr.replace(' ',',')
    connection  = sqlite3.connect(sherlockDB)
    cursor      = connection.cursor()
    cursor.execute('select %s from %s where not deleted and %s order by %s'%(colstr,jointTable,cnstr,order))
    mytable     = prettytable.from_db_cursor(cursor)
    cursor.close();connection.close()

    print mytable
    print 'Results: %d'%len(mytable._rows)

