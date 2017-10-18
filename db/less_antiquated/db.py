import sqlite3,time,os

"""
Primary functions for interfacing with main database, $SCRATCH/db/data.db 

To do: automate updating of lastModified field?
e.g. def lastModified(jobID): updateDB('bulkjob','lastmodified','id',jobID,time.time())
"""

sherlockDB = '/scratch/users/ksb/db/data.db'
suncatDB   = '/nfs/slac/g/suncatfs/ksb/db/data.db'

def getDBpath():
	c = getCluster()
	dbLocation = {'sherlock':sherlockDB,'suncat':suncatDB}
	return dbLocation[c]

############################
## MODIFY Functions
############################
def sqlexecute(sqlcommand, binds = []):
	dbpath = getDBpath()
	connection = sqlite3.connect(dbpath,timeout=600)
	cursor = connection.cursor()
	cursor.execute(sqlcommand, binds)
	if sqlcommand.lower()[:6] == 'select' or sqlcommand.lower()[:4] == 'with':
		output_array = cursor.fetchall()
		connection.close()
		return output_array
	else:
		connection.commit()
		connection.close()
		return

def sqlEq(x): 
	"""
	Custom SELECT command for each object in database. 
	Returns boolean if an object is found with a specific set of attributes that constitute equality
	"""
	col_vals = zip(x.sqlCols(),x.sqlInsert())
	def isStr(y): 		return isinstance(y,str) or isinstance(y,unicode)
	def formatStr(z): 	return '"'+z+'"' if isStr(z) else z
	return "SELECT EXISTS(SELECT 1 FROM {0} WHERE {1} LIMIT 1)".format(x.sqlTable(),
			' AND '.join(['{0}={1}'.format(k,formatStr(v)) for k,v in col_vals if k in x.sqlEq()]))

def insertObject(obj,dbpath = sherlockDB,tentative=False):
	"""
	Inserts an object into its proper table after checking that it doesn't already exist 
	For bulkjobs, we want the ability to tentatively add 
	"""
	query_out = [[0]] if tentative else sqlexecute(sqlEq(obj)) 

	if query_out[0][0]==0: #it's new
		connection = sqlite3.connect(dbpath)
		cursor = connection.cursor()
		
		try:				table = obj.sqlTable(tentative)
		except TypeError: 	table = obj.sqlTable()

		cols  = ','.join(obj.sqlCols())
		questionMarks = ','.join(['?']*len(obj.sqlCols()))
		vals = obj.sqlInsert()
		
		cursor.execute('insert into {0}({1}) values({2})'.format(table,cols,questionMarks), vals)
		ID = cursor.lastrowid
		connection.commit()
		connection.close()
		status = 'inserted'
	else:
		assert query_out[0][0]==1, "Weird result for query_out in insertObject: "+str(query_out)
		ID,status = 1,'already in db'; print status

	return ID,status

def updateDB(tableName,setColumn,idColumn,ID,newValue,prevValue=None):
	"""
	Change a particular cell of the database. Optionally, require that previous value be confirmed for what you think it is.
	"""
	realVal = query1(tableName,setColumn,idColumn,ID)
	if prevValue is not None: 
		assert realVal==prevValue, 'updateDB presumed item %s in %s had %s = %s, but instead it is %s'%(ID,tableName,setColumn,prevValue,realVal)
	sqlexecute("UPDATE {0} SET {1}= ? WHERE {2} = ?".format(tableName,setColumn,idColumn),[newValue,ID])

def deleteTable(tableName): sqlexecute("DELETE FROM %s"%tableName) #clear table
def deleteRow(tableName,idColumn,ID): sqlexecute("DELETE FROM {0} WHERE {1} = {2}".format(tableName,idColumn,ID))

##################
# QUERY FUNCTIONS
##################

def query1(tableName,columnOut,testColumn,value):
	command = "SELECT {1} FROM {0} WHERE {2} = ? ".format(tableName,columnOut,testColumn)
	output= sqlexecute(command,[value])
	assert len(output) == 1, 'query1 returned more than one result '+str(output)
	return output[0][0]

def query1all(tableName,testColumn,value):
	command = "SELECT * FROM {0} WHERE {1} = ? ".format(tableName,testColumn)
	output= sqlexecute(command,[value])
	assert len(output) == 1, 'query1all returned more than one result '+str(output)
	return output[0]

def queryOneFromMany(tableName,columnOut,testColumn,value):
	command = "SELECT {1} FROM {0} WHERE {2} = ? ".format(tableName,columnOut,testColumn)
	output= sqlexecute(command,[value])
	assert len(output[0])==1, 'queryMany doesn''t output list of singleton tuples?'
	return [o[0] for o in output]

def queryManyFromOne(tableName,columnsOut,testColumn,value):
	command =  "SELECT {1} from {0} WHERE {2} = {3}".format(tableName,','.join(columnsOut),testColumn,value)
	print command
	output = sqlexecute(command)
	return output[0]

def countDB(tableName): 
	"""
	returns [1,2,3, ... #number of rows] for a given table
	"""
	connection = sqlite3.connect(sherlockDB)
	cursor = connection.cursor()
	rowsQuery = "SELECT Count() FROM %s" % tableName
	cursor.execute(rowsQuery)
	numberOfRows = cursor.fetchone()[0]
	connection.commit() # necessary for select statement?
	connection.close()
	return numberOfRows


##############################
### BULKJOB SPECIFIC FUNCTIONS
##############################

def updateStatus(table,ID,status,newStatus): # = updateDB('bulkjob','status','id',ID,newstatus,status)
	realstatus = query1(table,'status','id',ID)
	assert realstatus==status, 'updateStatus presumed job#%d had status %s, but instead it has status %s'%(ID,status,realstatus)
	sqlexecute("UPDATE bulkjob SET status= ? WHERE id = ?",[newStatus,ID])


##################################
# Constraint using functions
##################################
def plotQuery(table,columnsOut,constraints):
	tableDict = {'bulkresult':	 	" FROM bulkjob LEFT JOIN bulkresult on bulkresult.jobid=bulkjob.id"
				,'tentative': 		" FROM tentative LEFT JOIN bulkresult on bulkresult.jobid=tentative.id" 
				,'tentativesurf': 	" FROM tentativesurf LEFT JOIN surfaceresult on surfaceresult.surfjobid=tentativesurf.id" 
				,'surfaceresult': 	" FROM surfacejob LEFT JOIN surfaceresult on surfaceresult.surfjobid=surfacejob.id" }
	tableName = tableDict[table]
	command =  "SELECT "+','.join(columnsOut)+tableName+' WHERE 1 '+ ''.join([' AND '+str(x) for x in constraints]) # 1==True in case there are no constraints
	print command
	return sqlexecute(command)

def getJobIDs(table,constraints): 
		if 'bulk' in table: 	tableInput,col 	= 'bulkresult',		['bulkjob.id']
		elif 'surf' in table: 	tableInput,col 	= 'surfaceresult',	['surfacejob.id']
		else: raise ValueError ,'Bad table argument for getJobIDs: ',table
		return [x[0] for x in plotQuery(tableInput,col,constraints)]

######################
# SUPPORTING FUNCTIONS
######################

def getCluster():
	hostname = os.environ['HOSTNAME'].lower()
	if 'sh' in hostname: return 'sherlock'
	elif 'su' in  hostname: return 'suncat' #important to distinguish suncat2 and 3?
	elif 'kris' in hostname: return 'kris'
	else: raise ValueError, "getCluster did not detect SH or SU in %s"%hostname

def getASEdb():
	aseLocation = {'sherlock':'/scratch/users/ksb/db/ase.db','suncat':'/nfs/slac/g/suncatfs/ksb/db/ase.db'}
	return aseLocation[getCluster()]

####################
# Function graveyard
####################

