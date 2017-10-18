"""
Perform a hard reset on everything
"""

import os

if raw_input('are you sure?') in ['y','yes']:

	print "Deleting data.db"
	try: os.remove('/scratch/users/ksb/db/data.db')
	except OSError: pass

	print "Deleting ase.db"
	try: os.remove('/scratch/users/ksb/db/ase.db')
	except OSError: pass

	print "Initializing new data.db from createDB.py"
	import createDB
	createDB.main()

	print "Populating convergence table of data.db"
	import convergence
	convergence.main()

	print "initializing ase.db from createTrajs.py, using data from Keld's data_solids_wPBE"
	import createBulks
	createBulks.main()
	
	#print "adding new initialized jobs to db from initializeJobs.py"
	#import initializeBulkJobs
	#initializeBulkJobs.main()
