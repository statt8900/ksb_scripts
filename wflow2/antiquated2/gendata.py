#external modules
import json,pickle,copy
import ase.visualize as viz
# internal modules
import jobs,dbase,details,readJobs,surfFuncs,interstitialFuncs,misc,gas
from printParse import *
"""
from jobs 				import jobs.assignJob
from printParse 		import *
from dbase 				import dbase.query
from details	 		import load
from surfFuncs 			import bulk2surf,adsorbedSurface
from interstitialFuncs 	import getInterstitials
from misc 				import nonmetals
from gas 				import aseMolDict
"""
#################
# CREATE NEW JOBS
#################
############################################################################
############################################################################



def varies(constraint,param,rang):
	readJobs.readJobs()
	jbs 		= []
	defaults 	= RELAXORLAT
	cnst 		= AND([defaults,constraint])
	output 		= dbase.queryCol('params_json',cnst)

	if ask('Do you want to vary %s over range %s for %d jobs?'%(param,str(rang),len(output))):
		for params in map(json.loads,output):
			for v in rang:
				ps = copy.deepcopy(params)
				assert not (jobs.assignJob(ps).new()),'varying a job that isn\'t complete or doesn\'t exist'
				ps[param] = v
				job = jobs.assignJob(ps)
				if job.new(): jbs.append(job)
		
		
		if ask('Do you want to launch %d new jobs?'%len(jbs)):
			for j in jbs: j.check();j.submit()
			misc.launch()

def vary(constraint,param,rang,check=True):
	# verify that constraint uniquely picks out a job
	readJobs.readJobs()
	defaults = RELAXORLAT
	output = dbase.queryCol('params_json',AND([defaults,constraint]))
	assert len(output)==1, 'Constraint %s is not specific enough, the following sets of params are returned: %s'%(constraint,'\n\n'.join([abbreviateDict(z) for z in output]))
	params = json.loads(output[0])
	assert param in params.keys(), 'param %s not in keys: %s'%(abbreviateDict(param),params.keys())
	print 'Param dict with %s being varied over the range %s: %s'%(param,str(rang),abbreviateDict(params))

	jbs=[]
	for v in rang:
		ps = copy.deepcopy(params)
		assert not (jobs.assignJob(ps).new()),'varying a job that isn\'t complete or doesn\'t exist'
		ps[param] = v
		job = jobs.assignJob(ps)
		if job.new(): jbs.append(job)

	if not check or ask('Do you want to launch %d new jobs?'%len(jbs)):
		for j in jbs: j.check();j.submit();
		if check: misc.launch()
############################################################################
############################################################################
def molecule(molname=None):
	readJobs.readJobs()
	for name,m in gas.aseMolDict.items():
		if molname is None or name == molname:
			mol = pickle.dumps(m)
			variableParams = {'pw': 		500
							,'xc': 			'RPBE'
							,'psp':			'gbrv15pbe'
							,'fmax': 		0.05
							,'dftcode':		'quantumespresso'}

			defaultParams = {'jobkind': 'relax'
							,'inittraj_pckl': 	mol
							,'comments': 		'created by gendata.molecule()'
							,'trajcomments':	''
							,'name': 			name
							,'relaxed': 		0
							,'kind': 			'molecule'
							,'kptden': 		1 #doesn't matter, will be ignored. Use 1 to be safe.
							,'dwrat': 10,'econv':5e-4,'mixing':0.1,'nmix':10,'maxstep':500,'nbands':-12,'sigma':0.1}
			p = misc.mergeDicts([variableParams,defaultParams])
			if molname is None or ask('Do you want to run a gas phase calculation with these params?\n%s'%(abbreviateDict(p))):
				j = jobs.assignJob(p)
				if j.new(): j.check();j.submit();
	misc.launch()
############################################################################
############################################################################

def calculateBulkModulus(constraint="1"):

	details.load(['finaltraj_pckl'],incomplete=False)

	cons 	= AND([COMPLETED,LATTICEOPT,constraint]) # necessary for bulkmodulus calc to be valid


	output 	= dbase.query(['fwid','params_json','finaltraj_pckl'],cons)
	question = 'Are you sure you want to calculate bulk modulus for %d structures?'%len(output)
	if ask(question):
		newjbs=[]
		for fwid,paramStr,ftraj in output:
			param 	  = json.loads(paramStr)
			del param['xtol']
			param['jobkind']		= 'bulkmod'
			param['strain']			= 0.03
			param['relax'] 			= fwid
			param['inittraj_pckl']	= misc.restoreMagmom(ftraj)
			job = jobs.assignJob(param)
			if job.new(): newjbs.append(job)
		if ask("launch %d new jobs?"%len(newjbs)):
			for j in newjbs: j.check();j.submit()
			misc.launch()

############################################################################
############################################################################

def getXCcontribs(constraint="1"):
	details.load(['finaltraj_pckl'],incomplete=False)
	cons 	= AND([COMPLETED,GPAW,RELAXORLAT,constraint]) 
	output 	= dbase.query(['fwid','params_json','finaltraj_pckl'],cons)
	question = 'Are you sure you want to calculate XC contributions for %d structures?'%len(output)
	if ask(question):
		newjbs=[]
		for fwid,paramStr,ftraj in output:
			param 	  = json.loads(paramStr)
			param['jobkind'] 		= 'xc'
			param['inittraj_pckl']	= misc.restoreMagmom(ftraj)
			param['relax'] 			= fwid

			job = jobs.assignJob(param)
			if job.new(): newjbs.append(job)
		if ask("launch %d new jobs?"%len(newjbs)):
			for j in newjbs: j.check();j.submit()
			misc.launch()


############################################################################
############################################################################

def getVibs(constraint="1"):
	cons 	= AND([QE,COMPLETED,RELAXORLAT,NONMETAL,constraint])
	details.load('finaltraj_pckl',cons)
	output 	= dbase.query(['fwid','params_json','finaltraj_pckl','natoms','numbers_pckl'],cons)
	
	question 		= 'Are you sure you want to calculate vibrations for %d structures?'%len(output)
	if ask(question):
		newjbs=[]
		for fwid,paramStr,ftraj,n,ns in output:
			param 	= json.loads(paramStr)
			nums 	= zip(range(n),pickle.loads(ns))
			nonmetal= [x for x,y in  nums if y in misc.nonmetals] #default do vibrations only on nonmetal atoms
			
			if len(nonmetal) > 0:
				param['jobkind'] 		= 'vib'
				param['relax'] 			= fwid
				param['inittraj_pckl']	= misc.restoreMagmom(ftraj)
				param['vibids_json'] 	= json.dumps(nonmetal)
				param['delta']			= 0.04
				job = jobs.assignJob(param)
			
				print abbreviateDict(param)

				if job.new(): newjbs.append(job)

		if ask("launch %d new jobs?"%len(newjbs)):
			for j in newjbs: j.check();j.submit()
			misc.launch()

############################################################################
############################################################################

def getBareSlab(constraint="1",facet=[1,1,1],xy=[1,1],layers=4,constrained=2,symmetric=0,vacuum=10,vacancies=[]):

	cons 	= AND([COMPLETED,LATTICEOPT,constraint])
	details.load('finaltraj_pckl',cons)
	output 	= dbase.query(['fwid','params_json','finaltraj_pckl'],cons)
	
	question = 'Are you sure you want to create bare slabs for %d structures?'%len(output)
	if ask(question):
		newjbs=[]
		for fwid,paramStr,ftraj in output:
			param 	 = json.loads(paramStr)
			
			surf,img = surfFuncs.bulk2surf(ftraj,facet,xy,layers,constrained,symmetric,vacuum,vacancies)
			
			param['jobkind'] 		= 'relax'
			param['kind'] 			= 'surface'
			param['name'] 			+= '_'+','.join(map(str,facet))+'_'+'x'.join(map(str,(xy+[layers])))
			param['relax'] 			= 0
			param['bulkparent'] 	= fwid
			param['inittraj_pckl']	= pickle.dumps(surf)
			param['sites_base64'] 	= img
			param['facet_json'] 	= json.dumps(facet)
			param['xy_json'] 		= json.dumps(xy)
			param['layers'] 		= layers
			param['constrained'] 	= constrained
			param['symmetric'] 		= symmetric
			param['vacuum'] 		= vacuum
			param['vacancies_json'] = json.dumps(vacancies)
			param['adsorbates_json']= json.dumps({})
			
			job = jobs.assignJob(param)
			if job.new():
				viz.view(surf)
				question = 'Does this structure look right?\n'+abbreviateDict(param)
				if ask(question):
					job.check();job.submit()
		misc.launch()


############################################################################
############################################################################

def addAdsorbate(constraint="1"):
	constrnts 		= AND([COMPLETED,RELAX,SURFACE,SYMMETRIC(False),constraint])

	ads 			= {'H':['O1']}
	details.load('finaltraj_pckl',constrnts)
	output 			= dbase.query(['fwid','params_json','finaltraj_pckl'],constrnts)
	question 		= 'Are you sure you want to add adsorbates to %d slabs?'%len(output)

	if ask(question):
		for fw,paramStr,ftraj in output:
			params = json.loads(paramStr)
			
			newsurf = surfFuncs.adsorbedSurface(ftraj,json.loads(params['facet_json']),ads)
			if jobs.assignJob(params).spinpol(): newsurf.set_initial_magnetic_moments([3 if e in misc.magElems else 0 for e in newsurf.get_chemical_symbols()])
			ase.visualize.view(newsurf)
			
			params['name']				+= '_'+printAds(ads)
			params['surfparent'] 		= fw
			params['inittraj_pckl'] 	= pickle.dumps(newsurf)
			params['adsorbates_json'] 	= json.dumps(ads)
			job 						= jobs.assignJob(params)
			if job.new():
				viz.view(newsurf)
				question = 'Does this structure look right?\n'+abbreviateDict(params)
				if ask(question):
					job.check();job.submit()
		misc.launch()


############################################################################
############################################################################
def addInterstitial(constraint="1",load=True,emttol=0.2):

	inter,num = 'H',2

	cons 	= AND([COMPLETED,LATTICEOPT,QE,constraint])
	
	if load: details.load(['finaltraj_pckl'],cons)
	else: readJobs.readJobs()

	output 	= dbase.query(['fwid','params_json','finaltraj_pckl'],cons)
	
	question 		= 'Are you sure you want to add interstitials to %d structures?'%len(output)
	if ask(question):
		jbs = []
		for fwid,paramStr,ftraj in output:
			param 	= json.loads(paramStr)
			spnpl	= jobs.assignJob(param).spinpol()

			param['structure'] 		= 'triclinic'
			param['jobkind'] 		= 'vcrelax'
			param['relaxed'] 		= fwid
			if 'xtol' in param.keys(): del param['xtol']

			trajs = [[(pickle.loads(ftraj),'')]]

			for i in range(num):


				lastround = trajs[-1]
				trajs.append([])
				for inputtraj,strname in lastround:

					for newtraj,newname in interstitialFuncs.getInterstitials(inputtraj,inter,spnpl):

						trajs[-1].append((newtraj,strname+newname))


			def modParams(par,trj,nam):
				p = copy.deepcopy(par)
				p['name'] 			+= nam
				p['inittraj_pckl']	= pickle.dumps(trj)
				return p
			
			onelevel = [item for sublist in trajs[1:] for item in sublist]

			tentativejoblist = [jobs.assignJob(modParams(param,t,n)) for t,n in onelevel]
			addJobs,efilt = [],0
			for j in tentativejoblist:
				if any([ abs(j.emt() - x.emt()) < emttol for x in addJobs]): efilt+=1
				else: addJobs.append(j)
			jbs.extend(addJobs)
		
		newjbs = [j for j in jbs if j.new()]
		for j in newjbs:
			print j.symbols(),j.emt()
		check = ask('Do you want to check %d/%d new jobs? (%d filtered by emt)'%(len(newjbs),len(jbs)+efilt,efilt))
		if not check and ask('Do you want to exit?'): return 0
		for jb in newjbs:

			question = 'Does this structure look right?\n'+abbreviateDict(j.params)
			if check: viz.view(newtraj)
			if not check or ask(question): jb.check();jb.submit()
		misc.launch()
		