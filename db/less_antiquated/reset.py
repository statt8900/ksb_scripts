"""
Perform a hard reset on everything
"""

import os

check = raw_input('are you sure?')

if check in ['y','yes']:
	print "Deleting data.db"
	try: os.remove('/scratch/users/ksb/db/data.db')
	except OSError: pass

	print "Deleting ase.db"
	try: os.remove('/scratch/users/ksb/db/ase.db')
	except OSError: pass

	print "Initializing new data.db from createDB.py"
	import createDB
	createDB.main()

	print "creating ase.db from createTrajs.py"
	import createTrajs
	createTrajs.main()

	print "adding new initialized jobs to db from initializeJobs.py"
	import initializeJobs
	initializeJobs.main()
