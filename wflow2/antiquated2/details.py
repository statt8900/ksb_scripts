#External Modules
import subprocess,re,sys,os,json,pickle,sys,datetime,ase,ase.db,glob,sqlite3,ast,fireworks,prettytable,pymatgen,copy
import numpy 						as np
import ase.thermochemistry 			as asethermo
import pymatgen.symmetry.analyzer  	as psa
import pymatgen.io.ase 				as pmgase
import ase.visualize 				as viz
# Internal Modules
import data_solids_wPBE,constraint,jobs,misc,readJobs
import dbase 		as db		
from printParse import *

# Constant Terms
nELEM		= len(ase.data.chemical_symbols)
toAtomicNum = {symb:num for symb,num in zip(ase.data.chemical_symbols,range(nELEM))}
lpad 		= fireworks.LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')

#
class Detail(object):
	def __init__(self,colname,datatype,rank,none,inputcols,func,axislabel=None,derivConv=None):
		if axislabel is None: axislabel = colname # Default axis label
		self.colname 		= colname 			# column name in database
		self.datatype 		= datatype 			# SQLITE datatype
		self.rank 			= rank 				# order in which columns should be calculated
		self.none 			= none 				# Boolean: skip calculation if any arguments are none
		self.inputcols 		= inputcols 	  	# Dependencies on other columns
		self.func 			= func 				# Function applied to input columns to get result
		self.axislabel 		= axislabel 		# For plotting
		self.derivConv 		= derivConv 		# Dictionary with convergence criteria (only for plotting)
	def __str__(self): return self.colname

	def dependencies(self,existing=[]):
		global dDict
		newdeps = [d for d in self.inputcols if d not in existing+misc.initcols]
		for d in newdeps:
			newdeps+= dDict[d].dependencies(existing+newdeps)
		return newdeps

#####################
# Auxillary Functions
#####################
def identity(x): 			return x
def div(x,y): 				return x / y
def diff(f1=identity,f2=identity): 	return lambda x,y:  None if f1(x) is None or  f2(y) is None else f1(x) - f2(y)
def intersect(x,y): 		return list(set(x) & set(y))
def getFWpckl(fwid): 		return pickle.dumps(lpad.get_fw_dict_by_id(fwid))


def keyList(xs,postFunc=identity,preFunc=identity): 
	def subFunc(dic):
		if dic is None: return None
		if not isinstance(dic,dict): print 'dic ',dic
		out = preFunc(dic)
		for x in xs: 
			try:
				out=out[x]
			except KeyError: return None
		return postFunc(out)
	return subFunc

def key(x,postFunc=identity,preFunc=identity): return keyList([x],postFunc,preFunc)
			
def fileExists(filename): 	return lambda pth: os.path.exists(os.environ['ALL_FWS']+pth+'/'+filename)

def getFromPckl(f=identity): return lambda pckl: f(pickle.loads(pckl))	



def getFromJob(f=identity):  
	def fromJob(params):
		try: 
			return f((jobs.assignJob(json.loads(params))))
		except AttributeError: return None
	return fromJob
def handleError(jsn,status): return json.loads(jsn)['error'] if status == 'failed' else None
def handleTrace(jsn,status): return json.loads(jsn)['trace'] if status == 'failed' else None

def lastLaunch(fwpckl):
	launches = pickle.loads(fwpckl)['launches']
	if len(launches) > 0:
		return launches[-1]
	else: return None

def getResult(launchdir): 
	try:
		if 'scratch' in launchdir:
			with open(launchdir+'/result.json','r') as f: return f.read()
		else: 
			cat 	 = subprocess.Popen(['ssh','ksb@suncatls1.slac.stanford.edu', 'cat %s'%(launchdir)+'/result.json'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err =  cat.communicate()
			
			return out if len(out)>0 else None
	except IOError: return None
	""" if fwpckl is input:
	ll = lastLaunch(fwpckl)
	if ll is not None:
		return json.dumps(ll['action']['stored_data']) if ll['action'] is not None else None
	else: return None
	"""

def getStoichDirectly(launchdir,deleteFlag):
	fwid 	= getFWID(launchdir,deleteFlag)
	fwpckl 	= getFWpckl(fwid)
	params 	= json.loads(pickle.loads(fwpckl)['spec']['params_json'])
	job 	= jobs.assignJob(params)
	return job.numbers()

def getEnergyDirectly(launchdir,deleteFlag):
	fwid 	= getFWID(launchdir,deleteFlag)
	fwpckl 	= getFWpckl(fwid)
	return pickle.loads(fwpckl)['launches'][0]['action']['stored_data']['e0']


def molarfrac(symb):
	def subFunc(strlist): 
		symbs = ast.literal_eval(strlist)
		return symbs.count('H')/float(len(symbs))
	return subFunc

def getG(eform,kind,vibeng,inittraj,pointgroup,spin):
	if kind == 'molecule':
		if pointgroup in ['D*h','C*v']: geometry = 'linear'
		elif pointgroup == 'K*h': geometry = 'monatomic'
		else: geometry = 'nonlinear'
		symnum = misc.pointGroupToSymnum(pointgroup)

		thermo = asethermo.IdealGasThermo(vib_energies=vibeng,atoms=pickle.loads(inittraj,potentialenergy=eform)
										,geometry=geometry,symmetrynumber=symnum,spin=spin)
	else:
		thermo = asethermo.HarmonicThermo(vib_energies=vibeng,potentialenergy=eform)

	G = thermo.get_free_energy(temperature=298.15, pressure=101325.)
	return G


################
# Data interface
################

def keldData(k): 	
	#Presumes only col requested is 'name'. Accesses key of data entry in Keld's dataset
	def getKeldData(name):
		try: 				return data_solids_wPBE.data[data_solids_wPBE.getKey[name] ][k]
		except KeyError: 	
			try: return data_solids_wPBE.data[data_solids_wPBE.getKey[name+' kittel'] ][k] #for 'bulkmodulus' / 'bulkmodulus kittel'
			except KeyError: return None
	return getKeldData

def errAFunc(a,aexpt,name):
	"""Calculations done on primitive cells - need to convert to unit cell"""
	if   'hcp' in name: multFactor = 1
	elif 'bcc' in name: multFactor = 3**(0.5)/2.
	else:				multFactor = 2**(-0.5) #trigonal-shaped primitive cells
	return a - multFactor*aexpt

########
# Status
########

def getStatus(fwstatus,fwpckl):
	if fwstatus.lower() not in ['running','fizzled']: return fwstatus.lower()
	
	else:
		ll = lastLaunch(fwpckl)

		starttime   = ll['time_start']
		format 		= '%Y-%m-%dT%H:%M:%S.%f'
		start 		= datetime.datetime.strptime(starttime,format)
		now 	    = datetime.datetime.utcnow()
		delta 	    = (now - start).total_seconds()
		walltime 	= pickle.loads(fwpckl)['spec']['_queueadapter']['walltime']
		walltimelist= walltime.split(':')
		wallseconds = sum([int(x)*y for x,y in zip(walltimelist,[3600,60,1])])  # HH:MM
		if delta > wallseconds: return 'timeout'
		else: 					return fwstatus.lower()

def getTrace(fw_pckl,status):
	ll = lastLaunch(fw_pckl)
	if status == 'COMPLETED' or ll is None: return None
	else:
		try: return ll['action']['stored_data']['_exception']['_stacktrace']
		except (TypeError, KeyError) as e: return None

def getHFrac(kind,symbols_pckl,metalstoich_pckl):
	if kind=='bulk':
		symbols,metalstoich=pickle.loads(symbols_pckl),pickle.loads(metalstoich_pckl)
		h = symbols.count('H')
		if metalstoich==[]: return None
		return h/float(len(metalstoich))
	else: return None



def getRefEnergy(fwid,jobkind,pw,xc,kptden,psp,dwrat,econv,fmax,dftcode,xtol,speciespckl,symbolspckl):#,deleteFlag):
	""" FILEPATH  -> FLOAT """ #other inputs just to make sure they're updated
	""" relies on all references having 'fwid'"""
	if jobkind not in ['latticeopt','relax','vcrelax']: return None
	def lookupRef(symb): 
		global refEngDict		

		refName 		= misc.refNames[symb]

		calcConstraint 	= AND([RELAXORLAT,PW(pw),XC(xc),PSP(psp),DWRAT(dwrat),ECONV(econv),FMAX(fmax),DFTCODE(dftcode)]
							+([] if xtol is None else [XTOL(xtol)])
							+([] if refName in misc.gasRefs.values() else [KPTDEN(kptden)]))

		try:
			normalizedEngs = refEngDict[symb][calcConstraint] #memoized
			return normalizedEngs
		except KeyError: 

			output = db.query(['fwid','raw_energy','symbols_pckl','kind','numbers_pckl','vibengs_pckl','inittraj_pckl','pointgroup','spin'],AND([COMPLETED,NAME(refName),calcConstraint]))
			if len(output)==0 or any([x[1] is None for x in output]):
				print 'Reference for %s not yet calculated with calc %s'%(refName,calcConstraint)
				return None
			elif len(output)>1: 
				print 'Bizarre: multiple references for %s with calc %s'%(refName,calcConstraint)
				fw1,fw2 = output[0][0],output[1][0]
				print 'First two fwids are %d and %d'%(fw1,fw2)
				print 'Using first reference found. Difference betweeen jobs is: '
				db.diffJob(fw1,fw2)


			eng,symbs 		= output[0][1],pickle.loads(output[0][2])
			normalize 		= symbs.count(symb)	
			kind 			= output[0][3]
			nums 			= pickle.loads(output[0][4])
			vibengs_pckl 	= output[0][5] if any([x in misc.nonmetals for x in nums]) else []
			inittraj 		= output[0][6]
			pointgroup 		= output[0][7]
			spin 			= output[0][8]
			g 				= getG(eng,kind,vibengs_pckl,inittraj,pointgroup,spin)

			normEG 	= eng/normalize,g/normalize

			refEngDict[symb][calcConstraint] = normEG
			return normEG

	refEngs	= {e:lookupRef(e) for e in pickle.loads(speciespckl)} # {'H':(-32.1,','Pt-fcc','Pt-fcc','Pt-fcc']
	if refEngs.values().count(None)>0: 	return None
	else: 
		e,g=0,0
		for sym in pickle.loads(symbolspckl):
			e+=refEngs[sym][0]
			g+=refEngs[sym][1]
		return json.dumps(e,g)

def getPointgroup(kind,inittraj): 
	if kind == 'molecule':
		atoms 	= pickle.loads(inittraj)
		cm 		= atoms.get_center_of_mass()
		m 		= pymatgen.Molecule(atoms.get_chemical_symbols(),atoms.get_positions()-cm)
		try:
			return psa.PointGroupAnalyzer(m).sch_symbol
		except ValueError:
			viz.view(atoms)
			sys.exit()
	else: return None

def getSpacegroup(kind,inittraj): 
	if kind == 'bulk':
		pmg = pmgase.AseAtomsAdaptor().get_structure(pickle.loads(inittraj))
		return psa.SpacegroupAnalyzer(pmg).get_space_group_number()
	else: return None

def getSpin(name): return 1 if name=='O2' else 0
def getEform(e,egref):
	return e - json.loads(egref)[0]
def getGform(e,egref):
	return e - json.loads(egref)[1]


#######################################################################################################################################################################
#######################################################################################################################################################################
def getTPS(fwpckl,resultjson):
	try: 				n = json.loads(resultjson)['niter'] 
	except KeyError: 	return None
	ll = lastLaunch(fwpckl)
	if ll is not  None: return ll['runtime_secs']/float(n)
	else: 				return None

allDetails = [
	Detail('blank',				'varchar',	0,False,['fwid'],				lambda x: '')

	,Detail('pointgroup',		'numeric',	0,False,['kind','inittraj_pckl'], 	getPointgroup)
	,Detail('spacegroup',		'numeric',	0,False,['kind','inittraj_pckl'], 	getSpacegroup)
	,Detail('spin',				'numeric',	0,False,['name'], 		getSpin)
	,Detail('fwpckl',			'varchar',	1,False,['fwid'],				getFWpckl)

	,Detail('hfrac',			'numeric',	1,False,['kind'
											,'symbols_pckl'
											,'metalstoich_pckl'],			getHFrac)

	,Detail('dof',				'integer',	1,False,['params_json'],		getFromJob(lambda x: x.dof()))

	,Detail('a_expt',			'numeric',	1,False,['name'],				keldData('lattice parameter'),	r'Expt Lattice A, $\AA$')

	,Detail('bulkmod_expt',		'numeric',	1,False,['name'],				keldData('bulk modulus'),		'Expt Bulk Mod, GPa')

	,Detail('isref',			'integer',	1,False,['name'],				lambda name: name in misc.refNames)

	,Detail('fwstatus',			'varchar',	2,False,['fwpckl'],				key('state',preFunc=pickle.loads))

	,Detail('result',			'varchar',	3,False,['launchdir'],			getResult)

	,Detail('raw_energy',		'numeric',	4,False,['result'],		 	key('raw_energy',preFunc=json.loads),'Raw Electronic Energy, eV',{'pw':(0.001,200)})
	,Detail('niter', 			'integer',	4,False,['result'],			key('niter',preFunc=json.loads),'Number of ionic steps')
	,Detail('finaltraj_pckl',	'varchar',	4,False,['result'],			key('finaltraj_pckl',preFunc=json.loads))
	,Detail('bfit',				'numeric',	4,False,['result'],			key('bfit',preFunc=json.loads))
	,Detail('bulkmodimg',		'varchar',	4,False,['result'],			key('bulkmodimg',preFunc=json.loads))
	
	,Detail('xccontribs',		'varchar',	4,False,['result'],			key('xcContribs',preFunc=json.loads))
	,Detail('uvib',				'numeric',	4,False,['result'],			key('uvib',preFunc=json.loads))
	,Detail('svib',				'numeric',	4,False,['result'],			key('svib',preFunc=json.loads))
	,Detail('vibengs_pckl',		'varchar',	4,False,['result'],			key('vibengs_pckl',preFunc=json.loads))
	,Detail('vibsummary',		'varchar',	4,False,['result'],			key('vibsummary',preFunc=json.loads))
	,Detail('zpe',				'numeric',	4,False,['result'],			key('zpe',preFunc=json.loads))
	,Detail('timeperstep',		'numeric',	4,False,['fwpckl'
													,'result'],			getTPS,'Time per ionic step, min')
	,Detail('uvib_minus_zpe',	'numeric',	5,False,['uvib,zpe'],		diff)
	,Detail('err_BM',			'numeric',	5,False,['result'
													,'bulkmod_expt'],	diff(f1=key('bulkmod',preFunc=json.loads)),'Bulk Modulus Error, GPa',{'pw':(100,200)})

	,Detail('status',			'varchar',	5,False,['fwstatus'
													,'fwpckl'],		getStatus)
	,Detail('calc_A',			'numeric',	5,False,['result'],			key('finalcell_pckl',lambda cell: misc.cell2param(pickle.loads(cell))[0],json.loads),	r'Calculated Lattice A, $\AA$')

	,Detail('trace',			'varchar',	6,False,['fwpckl'
													,'status'],			getTrace)
	,Detail('err_A',		 	'numeric',	6,False,['calc_A'
													,'a_expt'
													,'name'],			errAFunc,	r'Error in Lattice A, $\AA$',{'pw':(0.001,300)})


	##### REFERENCE DETAILS (Requiring database to have entries already for energy, etc.

	,Detail('refeng_json','numeric',100,False,['fwid','jobkind','pw','xc','kptden','psp','dwrat','econv','fmax','dftcode','xtol','species_pckl','symbols_pckl'],getRefEnergy)
	,Detail('eform','numeric',101,False,['raw_energy','refeng_json'],getEform,r'$\Delta E_{form},eV$')
	,Detail('raw_g',	'numeric',102,False,['eform','kind','vibengs_pckl'
											,'inittraj_pckl','pointgroup','spin'], getG)


	,Detail('gform','numeric',103,False,['raw_g','refeng_json'],getGform,r'$\Delta E_{form},eV$',r'$\Delta G_{form},eV$')
	,Detail('eformperatom','numeric',102,False,['eform','natoms'],div,r'$\Delta H_{form},eV/atom$')
	,Detail('surfaceenergy','numeric',103,False,['eform','surfacearea'],div,r'$\Delta H_{form},eV/\AA^2$')

	,Detail('pw',None,None,None,None,None,'PW cutoff')
	,Detail('xc',None,None,None,None,None,'Exchange-Correlation Functional')
	] 


global dDict
dDict = {x.colname:x for x in allDetails}
defaults = [d.colname for d in allDetails if d.rank<100]

def testRank():
	for d in allDetails:
		for dd in d.dependencies(): assert d.rank>dDict[dd].rank, '%s (%d) should have a higher rank than %s (%d)'%(d.colname,d.rank,dd,dDict[dd].rank)

####################################################################################
# Functions that USE details #######################################################
####################################################################################

def addKnowledge(colname,knowns,load_dependencies,colnames):
	""" adds col to output dictionary, assumes all dependencies already in dict """
	d,n = dDict[colname],len(knowns.values()[0])

	print 'adding '+colname

	def args(i): return [knowns[arg][i] for arg in d.inputcols]

	knowns[colname]=[]

	noskip = load_dependencies or colname in colnames

	if not noskip: print 'reading %s from table'%colname

	for i in range(n):
		print '%d/%d'%(i+1,n) ; sys.stdout.write("\033[F") # Cursor up one line
		if noskip:
			knowns[colname].append(d.func(*args(i)) if d.none or args(i).count(None) == 0 else None)	
		else:
			knowns[colname].append(db.query1(colname,'fwid',knowns['fwid'][i]))
	return knowns

def preload(verbose=False):
	global refEngDict
	refEngDict = {k:{} for k in misc.refNames.keys()}

	if verbose: print "adding colnames"

	for d in allDetails:
		try:  db.addCol(d.colname,d.datatype,'job')
		except sqlite3.OperationalError: pass

	if verbose: print "add jobs from directories"

	readJobs.readJobs()



def load(colnames=' '.join(defaults),cnst='1',loadDeps=True):
	"""For a certain query of some columns, required information is agglomerated then added to to db in one connection""" 

	print 'entering preload'
	preload()

	print 'entering load'

	colnames = [c for c in colnames.split() if c not in misc.initcols] #don't need to load launchdir

	outDict,allcols = {c:db.queryCol(c,cnst) for c in misc.initcols},copy.deepcopy(colnames)

	for c in colnames: 
		try: allcols+=dDict[c].dependencies(allcols) #get all dependencies  
		except KeyError as e: print 'keyerror: ',e
	if len(colnames)==0: return 0

	allcols.sort(key=lambda x: dDict[x].rank) 
	print 'allcols ',allcols	

	for c in allcols: outDict=addKnowledge(c,outDict,loadDeps,colnames)

	print "executing sql..."
	sqlinputs = zip(*[outDict[x] for x in colnames+['fwid']])
	for i,s in enumerate(sqlinputs):
		print '%d/%d'%(i+1,len(sqlinputs)) ; sys.stdout.write("\033[F") # Cursor up one line
		command = 'update job SET %s where fwid = ?'%(', '.join([x+' = ?' for x in colnames]))
		db.sqlexecute(command,s)
		#db.sqlexecutemany('update job SET %s where fwid = ?'%(', '.join([x+' = ?' for x in allcols])),zip(*[outDict[x] for x in allcols+['fwid']]))
	print "updated %d columns for %d rows"%(len(allcols),len(outDict.values()[0]))


def table(colstr='fwid',cnstr='1',update=False,tabstr='job'):
	if update: load(colstr,cnstr)
	colstr = colstr.replace(' ',',')
	connection 	= sqlite3.connect(db.sherlockDB)
	cursor 		= connection.cursor()
	cursor.execute('select %s from %s where %s'%(colstr,tabstr,cnstr))
	mytable 	= prettytable.from_db_cursor(cursor)
	cursor.close()
	connection.close()

	print mytable
	print 'Results: %d'%len(mytable._rows)



""" OTHER OPTIONS FOR DETAILS
	#,Detail('node',			'varchar',101,['lastlaunch_pckl'],getFromPckl(key('host')))
	#,Detail('qtime',			'numeric',	3,['lastlaunch_pckl'],	getFromPckl(key('reservedtime_secs')))
	#,Detail('created_on',		'varchar',	3,['fwpckl'],	getFromPckl(key('created_on')))
	#,Detail('updated_on',		'varchar',	3,['fwpckl'],	getFromPckl(key('updated_on')))
	#,Detail('nlaunches', 		'integer',	3,['fwpckl'],	getFromPckl(key('launches',len)))
	#,Detail('nodes',			'integer',5,['queueadapter'],	getFromJson(key('nodes')))
	#,Detail('ntasks',			'integer',5,['queueadapter'],	getFromJson(key('ntasks_per_node')))
	#,Detail('queue',			'varchar',5,['queueadapter'],	getFromJson(key('queue')))
	#,Detail('starttime',		'varchar',101,['lastlaunch_pckl'],getFromPckl(key('time_start')))
	#,Detail('nelec',			'integer',	1,False,['params_json'],	getFromJob(lambda x: misc.symbols2electrons(x.symbols())))
	#,Detail('lastlaunch_pckl',	'BLOB',		2,False,['fwpckl'],			getFromPckl(key('launches',lambda x: pickle.dumps(x[-1]) if x!=[] else None)))
	#,Detail('spec',		 	'varchar', 	2,False,['fwpckl'],			getFromPckl(key('spec',json.dumps)))
	#,Detail('queueadapter',	'varchar',	3,False,['spec'],			getFromJson(key('_queueadapter',json.dumps)))
	#,Detail('fworker',			'varchar',	3,False,['spec'],			getFromJson(key('_fworker')))
	#,Detail('runtime',			'numeric',	3,False,['lastlaunch_pckl'],getFromPckl(key('runtime_secs')))
	#,Detail('walltime',		'varchar',	4,False,['queueadapter'],	getFromJson(key('walltime')))
	#,Detail('forces_pckl',		'varchar',	4,False,['result'],			getFromJson(key('forces_pckl')))
	#,Detail('latticeopt_pckl',	'varchar',	4,False,['result'],			getFromJson(key('latticeopt_pckl')))
	#,Detail('finalcell_pckl',	'varchar',	4,False,['result'],			getFromJson(key('finalcell_pckl')))
	#,Detail('finalpos_pckl',	'varchar',	4,False,['result'],			getFromJson(key('finalpos_pckl')))
	#,Detail('bulkmod',			'numeric',	4,False,['result'],			getFromJson(key('bulkmod')))
"""

