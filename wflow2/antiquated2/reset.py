"""
Perform a hard reset on everything
"""
#External Modules
import os,shutil,datetime,printParse,fireworks
#Internal Modules
import dbase,createBulks,details

if True:#printParse.ask('are you sure?'):

	print "Deleting data.db"
	try: os.remove(os.environ['DATA_PATH'])
	except OSError: pass

	print "Deleting ase.db"
	try: os.remove(os.environ['ASE_PATH'])
	except OSError: pass


	if printParse.ask("remove FW_PATH as well?"):
		try:
			print "removing FW_PATH folder"
			assert os.environ['FW_PATH'] == '/scratch/users/ksb/fireworks/jobs/'
			shutil.rmtree(os.environ['FW_PATH'])
			os.mkdir(os.environ['FW_PATH'])
			print "removing ALL_FWS folder"
			assert os.environ['ALL_FWS'] == '/scratch/users/ksb/fireworks/alljobs/'
			shutil.rmtree(os.environ['ALL_FWS'])
			os.mkdir(os.environ['ALL_FWS'])

			print "removing contents of /nfs/slac/g/suncatfs/ksb/fireworks/jobs/*"
			os.system('ssh ksb@suncatls1.slac.stanford.edu rm -r /nfs/slac/g/suncatfs/ksb/fireworks/jobs/*')

		except OSError: pass
		

	if printParse.ask('remove FW database as well?'):

		lpad = fireworks.LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')
		if printParse.ask('really sure?'):
			lpad.reset(datetime.datetime.now().strftime('%Y_%m_%d').replace('_','-'))

	dbase.wipeDB()


	print "initializing ase.db from createTrajs.py,using data from Keld's data_solids_wPBE"
	createBulks.main()

	print "adding column headers"
	details.addCols()