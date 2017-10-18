import sqlite3,time,os

"""
Primary functions for interfacing with main database, $SCRATCH/db/test.db 
"""

sherlockDB = '/scratch/users/ksb/db/test.db'





def sqlexecutemany(sqlcommand, binds = []):
	print sqlcommand
	connection 	= sqlite3.connect(sherlockDB,timeout=60)
	cursor 		= connection.cursor()
	cursor.executemany(sqlcommand, binds)
	connection.commit()
	connection.close()
	return

def sqlexecute(sqlcommand, binds = []):
	connection = sqlite3.connect(sherlockDB,timeout=60)
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

def updateDB(setColumn,idColumn,ID,newValue,tableName):	
	sqlexecute("UPDATE {0} SET {1}= ? WHERE {2} = ?".format(tableName,setColumn,idColumn),[newValue,ID])

def dropTable(tableName): sqlexecute("DROP TABLE %s"%tableName)
def deleteRow(tableName,idColumn,ID): sqlexecute("DELETE FROM {0} WHERE {1} = {2}".format(tableName,idColumn,ID))
def addCol(colname,coltype,tableName): sqlexecute("ALTER TABLE {0} ADD COLUMN {1} {2}".format(tableName,colname,coltype))

def main():

	try: dropTable('detail0')
	except: pass
	sqlexecute("CREATE TABLE detail0 (id0 integer primary key)")
	addCol('col1','numeric','detail0')
	addCol('col2','numeric','detail0')
	addCol('col3','numeric','detail0')
	addCol('col4','numeric','detail0')
	for i in range(5): sqlexecute('insert into detail0 (col1) values (0,initialized)')

	rows = [2,3]
	cols = ['col3','col4']
	outDict = {'col3':[2.5,3.5],'col4':[1.0,-3]}
	print zip(*([outDict[x] for x in cols]+[rows]))
	sqlexecutemany('update detail0 SET %s where id0 = ?'%(', '.join([x+' = ?' for x in cols])),zip(*([outDict[x] for x in cols]+[rows])))

if __name__ == '__main__':
	main()


def addKnowledge(colname,knowns):
	""" adds col to dictionary, assumes all dependencies already in dict """
	d,n = dDict[colname],len(knowns.items()[0])
	knowns[colname]=[d.func(*[knowns[arg][i] for arg in d.inputcols] for i in range(n))]
	return knowns

