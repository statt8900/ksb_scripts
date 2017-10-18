from misc import readOnSherOrSlac,parseLine,cell2param
import databaseFuncs as db
##################
# Helper Functions
##################

def getAtoms(storedpath):
	import pickle
	return pickle.loads(readOnSherOrSlac(storedpath+'raw.pckl'))

def result(storagepath,key):
	import json
	return json.loads(readOnSherOrSlac(storagepath+'result.json'))[key]


##################
# From path
#################
def getJobtypeKSB(user,storagepath): 	return result(storagepath,'jobkind')if user=='ksb' else None
def getKindKSB(user,storagepath): 		return result(storagepath,'kind')   if user=='ksb' else None
def getKptdenKSB(user,storagepath): 	return result(storagepath,'kptden')   if user=='ksb' else None

def getStructureKSB(user,systype,storagepath):  return result(storagepath,'structure') 			if user=='ksb' and systype in ['bulk','surface'] else None
def getBscaleKSB(user,systype,storagepath): 	return result(storagepath,'bulkscale_json') 	if user=='ksb' and systype=='bulk' else None
def getBvacKSB(user,systype,storagepath): 		return result(storagepath,'bulkvacancy_json') 	if user=='ksb' and systype=='bulk' else None

def getXyKSB(user,systype,storagepath):  		return result(storagepath,'xy_json') 			if user=='ksb' and systype == 'surface' else None
def getSitesKSB(user,systype,storagepath):  	return result(storagepath,'sites_base64') 		if user=='ksb' and systype == 'surface' else None
def getFacetKSB(user,systype,storagepath):  	return result(storagepath,'facet_json') 		if user=='ksb' and systype == 'surface' else None
def getLayersKSB(user,systype,storagepath):  	return result(storagepath,'layers') 			if user=='ksb' and systype == 'surface' else None
def getVacuumKSB(user,systype,storagepath):  	return result(storagepath,'vacuum') 			if user=='ksb' and systype == 'surface' else None
def getSymmetricKSB(user,systype,storagepath):  return result(storagepath,'symmetric') 			if user=='ksb' and systype == 'surface' else None
def getVacanciesKSB(user,systype,storagepath):  return result(storagepath,'vacancies_json') 	if user=='ksb' and systype == 'surface' else None
def getAdsorbatesKSB(user,systype,storagepath): return result(storagepath,'adsorbates_json') 	if user=='ksb' and systype == 'surface' else None
def getConstrainedKSB(user,systype,storagepath):return result(storagepath,'constrained') 		if user=='ksb' and systype == 'surface' else None

def getQuadfitKSB(user,jobtype,storagepath):  return result(storagepath,'bfit') 	if user=='ksb' and jobtype=='bulkmod' 		else None
def getBulkmodKSB(user,jobtype,storagepath):  return result(storagepath,'bulkmod') 	if user=='ksb' and jobtype=='bulkmod' 		else None

def getBulkparentKSB(user,jobtype ,systype,storagepath): 	return result(storagepath,'bulkparent') if user=='ksb' and (jobtype =='vcrelax' or systype=='surface') else None

def getParentKSB(user,jobtype,storagepath): 
	if user!='ksb' or (jobtype not in ['vib','xc','dos']): return None
	try: 			 return result(storagepath,'parent')
	except KeyError: return None

def getNatoms(storedpath): return len(getAtoms(storedpath))

def getNumbersStr(storedpath): 
	import json
	return json.dumps(list(getAtoms(storedpath).get_atomic_numbers()))


def pawElem(numbrstr):
	import json
	return any([x in json.loads(numbrstr) for x in [3,4,11,12,19,20,30,37,38,55,56]])

def makeStrJob(user,storedpath):
	if user!='ksb': return None
	import jobs,json
	params = json.loads(readOnSherOrSlac(storedpath+'params.json'))
	return str(jobs.Job(params))

############
# FW Related
############

def getFWID(user,storedpath): 	
	if user != 'ksb': return None
	out = readOnSherOrSlac(storedpath+'FW_submit.script')
	return int(parseLine(out,'--fw_id').split()[-1])
	raise ValueError, 'No fw_id found in FW_submit.script: \n\n%s'%out

def fwdict(fwid): 
	import fireworks
	lpad = fireworks.LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')
	return lpad.get_fw_dict_by_id(fwid)

def lastLaunch(fwid): return fwdict(fwid)['launches'][-1]

def getFWstatus(fwid):
	if fwid is None: return None
	return fwdict(fwid)['state']

def getTPS(user,jobkind,storedpath,fwid):
	if user!='ksb' or jobkind not in ['latticeopt','relax']: return None
	n = readOnSherOrSlac(storedpath+'qn.log').count('BFGS')
	return lastLaunch(fwid)['runtime_secs']/float(n)


##########
# From log
##########
def getDFTcode(strpth):
    if   'espresso' in readOnSherOrSlac(strpth+'log'): return 'quantumespresso'
    elif 'gpaw'     in readOnSherOrSlac(strpth+'log'): return 'gpaw'
    else: raise NotImplementedError, 'Vasp?'

def getEng(strpth):    
    def getENG_GPAW(log): return float(parseLine(log,'Extrapolated:',-1).split()[-1])
    def getENG_QE(log):   return float(parseLine(log,'total energy',-2).split()[-2])
    dftcode = getDFTcode(strpth)
    if   dftcode == 'gpaw':            return getENG_GPAW(readOnSherOrSlac(strpth+'log'))
    elif dftcode == 'quantumespresso': return getENG_QE(readOnSherOrSlac(strpth+'log'))
    else: return None

def getPSPpath(storagepath,dftcode):
    def getPSPpath_QE(log):     return parseLine(log,'pseudo_dir').split("='")[1][:-2]
    def getPSPpath_GPAW(log):   return parseLine(log,'setups').split(': ')[-1].strip()
    if   dftcode == 'gpaw':            return getPSPpath_GPAW(readOnSherOrSlac(storagepath+'log'))
    elif dftcode == 'quantumespresso': return getPSPpath_QE(readOnSherOrSlac(storagepath+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'

#################################
# Linking with external databases
#################################

def keldData(k): 	
	import data_solids_wPBE
	
	#Presumes only col requested is 'name'. Accesses key of data entry in Keld's dataset
	def getKeldData(name):
		name = name.split('_')[0] #sometimes suffix ('bulkmod') is appended to name
		try: return data_solids_wPBE.data[data_solids_wPBE.getKey[name] ][k]
		except KeyError: 	
			try: return data_solids_wPBE.data[data_solids_wPBE.getKey[name] ][k+' kittel'] #for 'bulkmodulus' / 'bulkmodulus kittel'
			except KeyError: return None
	return getKeldData

getExptA = keldData('lattice parameter')
getExptBM = keldData('bulk modulus')




##################
### Analyzing cell
##################

def systemType(storedpath):
	from ase.data import covalent_radii
	atoms = getAtoms(storedpath)
	atomicVolume = sum([1.333*3.1415*covalent_radii[x]**2 for x in atoms.get_atomic_numbers()])
	fillFraction = atomicVolume / atoms.get_volume()
	if fillFraction > 0.3:    return 'bulk'
	elif fillFraction > 0.001: return 'surface'
	else: 					   return 'molecule'

def surfaceArea(storedpath,system_type,symmetric):
	if system_type == 'surface' and symmetric is not None:
		import numpy as np
		cell = getAtoms(storedpath).get_cell()
		multFactor = 2 if symmetric else None
		return  np.linalg.norm(np.cross(cell[0],cell[1]))
	else: return None

def getLatticeA(system_type,storedpath):
	if system_type == 'bulk':
		atoms = getAtoms(storedpath)
		return  cell2param(atoms.get_cell())[0]
	else: return None

def errBM(bm,bmexpt): 
	if bm is None or bmexpt is None: return None
	return bm - bmexpt

def errAFunc(a,aexpt,crystal):
	"""Calculations done on primitive cells - need to convert to unit cell"""
	if any([x is None for x in [a,aexpt,crystal]]): return None
	if    crystal in ['hcp','hexagonal']: 	multFactor = 1
	elif  crystal == 'bcc': 				multFactor = 3**(0.5)/2.
	else:									multFactor = 2**(-0.5) #trigonal-shaped primitive cells
	return a - multFactor*aexpt

def kx(stored_path,kpts_json):
	import json
	atoms = getAtoms(storedpath)
	return  json.loads(kpts_json)[0]/cell2param(atoms.get_cell())[0]
def ky(stored_path,kpts_json):
	import json
	atoms = getAtoms(storedpath)
	return  json.loads(kpts_json)[1]/cell2param(atoms.get_cell())[1]
def kz(stored_path,kpts_json):
	import json
	atoms = getAtoms(storedpath)
	return  json.loads(kpts_json)[2]/cell2param(atoms.get_cell())[2]

# Pymatgen
def getCrystalSystem(system_type,storedpath):
	"""triclinic,monoclinic,orthorhombic,tetragonal,trigonal,hexagonal, or cubic"""
	if system_type != 'bulk': return None
	from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
	from pymatgen.io.ase 			import AseAtomsAdaptor

	atoms = AseAtomsAdaptor().get_structure(getAtoms(storedpath))
	return  SpacegroupAnalyzer(atoms,symprec=0.1).get_crystal_system()

def getPointgroup(kind,storagepath): 
	if kind != 'molecule': return None
	import pymatgen
	import pymatgen.io.ase 			  as pmgase
	import pymatgen.symmetry.analyzer as psa
	atoms 	= getAtoms(storagepath)
	cm 		= atoms.get_center_of_mass()
	m 		= pymatgen.Molecule(atoms.get_chemical_symbols(),atoms.get_positions()-cm)
	return psa.PointGroupAnalyzer(m).sch_symbol

def getMolShape(pg): 
	if   pg is None: 		  return None
	elif pg == 'Dh': 		  return 'monatomic'
	elif pg in ['C*v','D*h']: return 'linear'
	else: 					  return 'nonlinear'

def getSpacegroup(kind,storagepath): 
	if kind != 'bulk':return None
	import pymatgen.io.ase 			  as pmgase
	import pymatgen.symmetry.analyzer as psa

	pmg = pmgase.AseAtomsAdaptor().get_structure(getAtoms(storagepath))
	return psa.SpacegroupAnalyzer(pmg).get_space_group_number()

def pointGroupToSymnum(pg):
	if pg is None: return None
	if pg in ['C1','Ci','Cs','C*v']: return 1
	elif pg == 'D*h': return 2
	elif pg in ['T','Td']: return 12
	elif pg == 'Oh': return 24
	elif pg == 'Ih': return 60
	else:
		letter,number = pg[0],pg[1]
		assert letter in ['S,C,D']
		assert number.isdigit()
		multdict = {'S':0.5,'C':1,'D':2}
		return multdict[letter]*int(number)

def surfaceEnergy(user):
	if user !='ksb': return None
	raise NotImplementedError

######################
# Equivalence of Calcs
#---------------------


def getEform(raw_energy,dftcode,xc,pw,kx,ky,kz,fmax,xtol,strain,psp,spinpol,numbstr):
	def getRefEng(symb,dftcode,xc,pw,kx,ky,kz,fmax,xtol,strain,psp,spinpol):
		c = "element = {0} and dftcode='{1}' and xc ='{2}' and pw={3} and fmax={4} and xtol={5} and strain={6} and psp={7} and spinpol={8}".format(dftcode,xc,pw,fmax,xtol,strain,psp,spinpol).replace('=None',' is null')
		if system_type != 'molecule': 
			c+="(system_type = 'molecule' or ((kptden_x BETWEEN {0} and {1}) AND (kptden_y BETWEEN {2} AND {3}) ".format(kx-1,kx+1,ky-1,ky+1)
			if system_type == 'bulk': c+= " AND (kptden_z BETWEEN {0} and {1})))".format(kz-1,kz+1)
			else: c+="))"

		print 'constraint ',c
			
		engs = db.sqlexecute('select refeng from refeng left join calc on refeng_calc = calc_id left join derived on refeng_job = derived_job where '%c)
		print engs
		sys.exit()
	refengDict = {symb:getRefEng(symb,dftcode,xc,pw,kpts_json,fmax,xtol,strain,psp,spinpol) for symb in list(set(json.loads(numbstr)))}
	return raw_energy - sum([refengDict[x] for x in numbstr])


