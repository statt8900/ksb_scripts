#external modules
import warnings
warnings.filterwarnings("ignore", message="Moved to ase.neighborlist")
warnings.filterwarnings("ignore", message="matplotlib.pyplot")

import json,pickle,copy
import ase.visualize as viz
# internal modules
import jobs,surfFuncs,interstitialFuncs,misc,gas,emt
import databaseFuncs as db
import manageIncompleteJobs as manage
from sqlShortcuts import *

#################
# CREATE NEW JOBS
#################
def ask(x): return raw_input(x+'\n(y/n)---> ').lower() in ['y','yes']

############################################################################
############################################################################

def varies(constraint,param,rang):
	"""
	varies("xc='BEEF' and name like '%Pd%'",'pw',range(500,1000,100)) 
		--> create jobs varying planewave_cutoff (500,600,...,1000) for all BEEF jobs with Pd
	"""
	jbs 		= []
	defaults 	= RELAXORLAT
	cnst 		= AND(defaults,constraint)
	output 		= db.queryCol('storage_directory',cnst)
	existing 	= manage.listOfIncompleteJobStrs()

	if ask('Do you want to vary %s over range %s for %d jobs?'%(param,str(rang),len(output))):
		for stordir in output:
			params = json.loads(misc.readOnSherOrSlac(stordir+'params.json'))
			for v in rang:
				ps = copy.deepcopy(params)
				assert not (jobs.Job(ps).new()),'varying a job that isn\'t complete or doesn\'t exist'
				ps[param] = v
				job = jobs.Job(ps)
				if job.new(existing): jbs.append(job)
		
		if ask('Do you want to launch %d new jobs?'%len(jbs)):
			for j in jbs: j.submit()
			misc.launch()

############################################################################
############################################################################

def molecule(molname=None):
	"""
	Give molecule name to choose initial structure, modify dictionary parameters for relaxation job
	"""
 
	existing = manage.listOfIncompleteJobStrs()

	for name,m in gas.aseMolDict.items():
		if molname is None or name == molname:
			mol = pickle.dumps(m)
			params = 		{'pw': 				500
							,'xc': 				'BEEF'
							,'psp':				'gbrv15pbe'
							,'fmax': 			0.05
							,'dftcode':			'quantumespresso'
							,'jobkind': 		'relax'
							,'inittraj_pckl': 	mol
							,'name': 			name
							,'relaxed': 		0
							,'kind': 			'molecule'
							,'kptden': 			1 #doesn't matter, will be ignored. Use 1 to be safe.
							,'dwrat': 10,'econv':5e-4,'mixing':0.1,'nmix':10,'maxstep':500,'nbands':-12,'sigma':0.1}

			if molname is None or ask('Do you want to run a gas phase calculation with these params?\n%s'%(misc.abbreviateDict(params))):
				jobs.Job(params).submit(existing)
			
	misc.launch()
############################################################################
############################################################################

def calculateBulkModulus(constraint="1"):

	cons 	= AND(LATTICEOPT,constraint) # necessary for bulkmodulus calc to be valid
	output 	= db.query(['fwid','storage_directory','job_name','structure_ksb','bulkvacancy_ksb','bulkscale_ksb','system_type_ksb','planewave_cutoff','xc','kptden_ksb','psp_ksb','dwrat_ksb','econv_ksb','mixing_ksb','nmix_ksb','maxstep_ksb','nbands_ksb','sigma_ksb','fmax_ksb','dftcode'],cons)
	existing = manage.listOfIncompleteJobStrs()

	question = 'Are you sure you want to calculate bulk modulus for %d structures?'%len(output)
	if ask(question):
		newjbs=[]
		for fwid,stor_dir,name,structure,bvj,bsj,kind,pw,xc,kptden,psp,dwrat,econv,mixing,nmix,maxstep,nbands,sigma,fmax,dftcode in output:

			params = {'jobkind':'bulkmod','strain':0.03
					,'inittraj_pckl': misc.restoreMagmom(misc.storageDir2pckl(stor_dir)) ,'name': name+'_bulkmod' 
					, 'structure':structure ,'bulkvacancy_json':bvj,'bulkscale_json':bsj,'kind': kind ,'dftcode': dftcode
					,'pw': pw ,'xc': xc ,'kptden': kptden ,'psp': psp
					,'dwrat': dwrat ,'econv': econv ,'mixing': mixing ,'nmix': nmix ,'maxstep': maxstep ,'nbands': nbands ,'sigma': sigma ,'fmax': fmax
					,'parent': fwid}
			job = jobs.Job(params)
	
			if job.new(existing): newjbs.append(job)
		if ask("launch %d new jobs?"%len(newjbs)):
			for j in newjbs: j.submit()
			misc.launch()

############################################################################
############################################################################

def getXCcontribs(constraint="1"):
	"""NOT YET DEBUGGED"""
	cons 	= AND(GPAW,RELAXORLAT,constraint) 
	output 	= db.query(['fwid','storage_directory','job_name','system_type_ksb','planewave_cutoff','xc','kptden_ksb','psp_ksb','dwrat_ksb','econv_ksb','mixing_ksb','nmix_ksb','maxstep_ksb','nbands_ksb','sigma_ksb','fmax_ksb','dftcode'],cons)
	question = 'Are you sure you want to calculate XC contributions for %d structures?'%len(output)
	if ask(question):
		newjbs=[]
		for fwid,stor_dir,name,kind,pw,xc,kptden,psp,dwrat,econv,mixing,nmix,maxstep,nbands,sigma,fmax,dftcode in output:
			params = {'jobkind':'xc' ,'inittraj_pckl': storageDir2pckl(stor_dir) ,'name': name+'_xc' ,'kind': kind ,'dftcode': dftcode
					,'pw': pw ,'xc': xc ,'kptden': kptden ,'psp': psp
					,'dwrat': dwrat ,'econv': econv ,'mixing': mixing ,'nmix': nmix ,'maxstep': maxstep ,'nbands': nbands ,'sigma': sigma ,'fmax': fmax
					,'parent': fwid}
			job = jobs.Job(params)
			if job.new() and ask('does this look right %s'%str(params)): newjbs.append(job)
		if ask("launch %d new jobs?"%len(newjbs)):
			for j in newjbs: j.submit()
			misc.launch()


############################################################################
############################################################################

def getVibs(constraint="1"):
	"""
	Launch vibration calculation for all latticeopt or (vc-)relax jobs that meet some (SQL) constraint in the shared database
	Only (and all) nonmetal atoms will be vibrated (defined in misc.py)
	"""

	existing = manage.listOfIncompleteJobStrs()
	
	cons 	= AND(QE,RELAXORLAT,constraint) #could we add a column to quickly filter things with nonmetals?
	
	output 	= db.query( ['fwid','storage_directory','job_name','system_type_ksb' 																									# generic parameters
						,'planewave_cutoff','xc','kptden_ksb','psp_ksb','dwrat_ksb','econv_ksb','mixing_ksb','nmix_ksb','maxstep_ksb','nbands_ksb','sigma_ksb','fmax_ksb','dftcode' # calc parameters
						,'structure_ksb','bulkvacancy_ksb','bulkscale_ksb' 																											# possible bulk parameters
						,'facet_ksb','xy_ksb','layers_ksb','constrained_ksb','symmetric_ksb','vacuum_ksb','vacancies_ksb','adsorbates_ksb','sites_ksb'],cons) 						# possible surf parameters
	
	question = 'Are you sure you want to calculate vibrations for %d structures?'%len(output)
	
	if ask(question):
		newjbs=[]
		for fwid,stor_dir,name,kind,pw,xc,kptden,psp,dwrat,econv,mixing,nmix,maxstep,nbands,sigma,fmax,dftcode,structure,bv,bs,facet,xy,lay,const,sym,vac,vacan,ads,sites in output:
			inittraj_pckl = storageDir2pckl(stor_dir)
			atoms 		  = pickle.loads(inittraj_pckl)
			
			nonmetal_inds = [i for i,x in enumerate(atoms.get_chemical_symbols()) if x in misc.nonmetalSymbs]
			
			if len(nonmetal_inds) > 0: 					
				params = {'jobkind':'vib' ,'inittraj_pckl': inittraj_pckl ,'name': name+'_vib' ,'kind': kind ,'dftcode': dftcode,'structure':structure
						,'pw': pw ,'xc': xc ,'kptden': kptden ,'psp': psp
						,'dwrat': dwrat ,'econv': econv ,'mixing': mixing ,'nmix': nmix ,'maxstep': maxstep ,'nbands': nbands ,'sigma': sigma ,'fmax': fmax
						,'structure':structure,'bulkvacancy_json':bv,'bulkscale_json':bs
						,'facet_json':facet,'xy_json':xy,'layers':lay,'constrained':const,'symmetric':sym,'vacuum':vac,'vacancies_json':vacan,'adsorbates_json':ads,'sites_base64':sites
						,'parent': fwid ,'vibids_json':json.dumps(nonmetal_inds) ,'delta': 0.04}
				
				job = jobs.Job(params)
				
				if job.new() and ask('does this look right %s'%str(params)): newjbs.append(job)
		
		if ask("launch %d new jobs?"%len(newjbs)):
			for j in newjbs: j.submit(existing)
			misc.launch()

############################################################################
############################################################################

def getBareSlab(constraint="1",check=True):
	"""Create a bare slab for some set of relaxed bulk structures. Alter parameters below as necessary"""
	cons 	= AND(RELAXORLAT,BULK,MBEEF,PAWPSP,BCC,PW(1500),KPTDEN(2),constraint)

	facet 		= [1,0,0]
	xy 			= [1,1]
	layers 		= 6
	constrained = 2
	symmetric 	= 1
	vacuum 		= 10
	vacancies 	= []

	output 	= db.query(['fwid','storage_directory','job_name','structure_ksb'
						,'planewave_cutoff','xc','kptden_ksb','psp_ksb','dwrat_ksb','econv_ksb','mixing_ksb','nmix_ksb','maxstep_ksb','nbands_ksb','sigma_ksb','fmax_ksb','dftcode'],cons)
	question = 'Are you sure you want to create bare slabs for %d structures?'%len(output)

	if ask(question):
		existing = manage.listOfIncompleteJobStrs()

		newjbs=[]
		for fwid,stor_dir,name,structure,pw,xc,kptden,psp,dwrat,econv,mixing,nmix,maxstep,nbands,sigma,fmax,dftcode in output:
			ftraj 	 = misc.readOnSherOrSlac(stor_dir+'raw.pckl')
			
			surf,img = surfFuncs.bulk2surf(ftraj,facet,xy,layers,constrained,symmetric,vacuum,vacancies)
			name+='_'+','.join(map(str,facet))+'_'+'x'.join(map(str,(xy+[layers])))
			params = {'jobkind':'relax' ,'inittraj_pckl': pickle.dumps(surf) ,'name': name ,'kind': 'surface' ,'dftcode': dftcode,'structure':structure # don't need bulk vacancies/scale?
					,'pw': pw ,'xc': xc ,'kptden': kptden ,'psp': psp,'bulkparent': fwid 
					,'dwrat': dwrat ,'econv': econv ,'mixing': mixing ,'nmix': nmix ,'maxstep': maxstep ,'nbands': nbands ,'sigma': sigma ,'fmax': fmax
					,'sites_base64': img,'facet_json': json.dumps(facet),'xy_json': json.dumps(xy),'layers': layers,'constrained':constrained,'symmetric': symmetric
					,'vacuum':vacuum,'vacancies_json': json.dumps(vacancies),'adsorbates_json': json.dumps({})}
			
			job = jobs.Job(params)
			if job.new(existing):
				if check: viz.view(surf)
				question = 'Does this structure look right?\n'+misc.abbreviateDict(params)
				if not check or ask(question): job.submit()
		misc.launch()

############################################################################
############################################################################

def addAdsorbate(constraint="1"):
	constrnts 		= AND(COMPLETED,RELAX,SURFACE,SYMMETRIC(False),constraint)

	ads 			= {'H':['O1']}

	output 			= db.query(['fwid','params_json','finaltraj_pckl'],constrnts)
	question 		= 'Are you sure you want to add adsorbates to %d slabs?'%len(output)

	if ask(question):
		for fw,paramStr,ftraj in output:
			params = json.loads(paramStr)
			
			newsurf = surfFuncs.adsorbedSurface(ftraj,json.loads(params['facet_json']),ads)
			if jobs.Job(params).spinpol(): newsurf.set_initial_magnetic_moments([3 if e in misc.magElems else 0 for e in newsurf.get_chemical_symbols()])
			ase.visualize.view(newsurf)
			
			params['name']				+= '_'+printAds(ads)
			params['surfparent'] 		= fw
			params['inittraj_pckl'] 	= pickle.dumps(newsurf)
			params['adsorbates_json'] 	= json.dumps(ads)
			job 						= jobs.Job(params)
			if job.new():
				viz.view(newsurf)
				question = 'Does this structure look right?\n'+abbreviateDict(params)
				if ask(question):
					job.check();job.submit()
		misc.launch()


############################################################################
############################################################################
def addInterstitial(constraint="1",emttol=0.2,vischeck=True,limit=1):
	"""
	Specify an interstitial and the number (up to which) you want to add them to the cell
	Filter duplicates via EMT-calculated energy (stipulate a tolerance)
	Assume final structure is triclinic
	"""

	inter,num = 'H',2

	cons 	= AND(LATTICEOPT,QE,PW(500),KPTDEN(2),constraint)

	output 	= db.query( ['fwid','storage_directory','job_name','system_type_ksb' 																									# generic parameters
						,'planewave_cutoff','xc','kptden_ksb','psp_ksb','dwrat_ksb','econv_ksb','mixing_ksb','nmix_ksb','maxstep_ksb','nbands_ksb','sigma_ksb','fmax_ksb','dftcode' # calc parameters
						,'bulkvacancy_ksb','bulkscale_ksb'],cons,limit=limit) 	# bulk parameters
	
	existing = manage.listOfIncompleteJobStrs() 	# check currently running/fizzled jobs to avoid duplicates
	
	if ask('Are you sure you want to create interstitials for %d structures?'%len(output)):
		newjbs,totJbs,totStruct=[],0,0 #initialize variables
		
		for fwid,stor_dir,name,kind,pw,xc,kptden,psp,dwrat,econv,mixing,nmix,maxstep,nbands,sigma,fmax,dftcode,bv,bs in output:
			
			params = {'jobkind':'vcrelax','kind': 'bulk' ,'dftcode': dftcode,'name':name
					,'pw': pw ,'xc': xc ,'kptden': kptden ,'psp': psp,'bulkparent': fwid 
					,'dwrat': dwrat ,'econv': econv ,'mixing': mixing ,'nmix': nmix ,'maxstep': maxstep ,'nbands': nbands ,'sigma': sigma ,'fmax': fmax
					,'structure':'triclinic','bulkvacancy_json':bv,'bulkscale_json':bs} #bulk
			
			ftraj  	= misc.storageDir2pckl(stor_dir)
			spnpl	= any([x>0 for x in pickle.loads(ftraj).get_initial_magnetic_moments()])
			
			trajs = [[(pickle.loads(ftraj),'')]] # initial state
			
			for i in range(num):
				
				lastround = trajs[-1] 				# for all structures with n - 1 interstitials (note: initially there is only one, with 0 interstitials)
				trajs.append([]) 					# initialize container for all new structures with n interstitials
				for inputtraj,strname in lastround:
					for newtraj,newname in interstitialFuncs.getInterstitials(inputtraj,inter,spnpl):
						trajs[-1].append((newtraj,strname+newname)) # add all new trajs found for all previous structures (append suffix to name)
			
			def modParams(par,trj,nam): 					# For a given new traj and name, create input parameters
				p = copy.deepcopy(par) 						# All parameters common to all of these jobs
				p['name'] 			+= nam 					# Unique job name
				p['inittraj_pckl']	= pickle.dumps(trj) 	# Initial structure
				return p
			
			def getEMT(x): 
				a = pickle.loads(x['inittraj_pckl'])
				a.set_calculator(emt.EMT()); 
				return a.get_potential_energy()

			onelevel = [modParams(params,item[0],item[1]) for sublist in trajs[1:] for item in sublist] #collapse list of lists to a single list of input parameters

			emtPairs = sorted(zip(map(getEMT,onelevel),onelevel)) # pairs of (Energy,Params) ordered by energy
			
			emtCounter,jbs = 0,[]
			for e,p in emtPairs:
				if e - emtCounter > emttol: 
					jbs.append(jobs.Job(p))
					emtCounter=e
			
			newjbs.extend([x for x in jbs if x.new(existing)])
			totJbs+=len(jbs) ; totStruct+=len(onelevel)
		
		check = ask('Do you want to check %d/%d new jobs? (%d filtered by emt)'%(len(newjbs),totStruct,totStruct-totJbs))
		if not check and ask('Do you want to exit?'): return 0
		for jb in newjbs:
			question = 'Does this structure look right?\n'+misc.abbreviateDict(jb.params)
			if vischeck: viz.view(pickle.loads(jb.params['inittraj_pckl']))
			if not vischeck or ask(question): jb.submit()
		misc.launch()
		