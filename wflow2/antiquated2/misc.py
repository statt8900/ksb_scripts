#External Modules
import os,itertools,ase,difflib,pprint,pickle
import numpy as np
#############################
"""Miscellaneous functions"""
##############################

def restoreMagmom(trajpckl):
	traj = pickle.loads(trajpckl)
	try:
		mags = traj.get_magmoms()
		if any([x>0 for x in mags]): traj.set_initial_magnetic_moments([3 if e in magElems else 0 for e in traj.get_chemical_symbols()])
	except: pass
	return pickle.dumps(traj)
def launch(): os.system('/home/ksb/scripts/fireworks/launcher.sh')

def dict_diff(first, second):
	""" Return a dict of keys that differ with another config object.  If a value is
		not found in one fo the configs, it will be represented by KEYNOTFOUND.
		@param first:   Fist dictionary to diff.
		@param second:  Second dicationary to diff.
		@return diff:   Dict of Key => (first.val, second.val)
	"""
	KEYNOTFOUND = '<KEYNOTFOUND>'       # KeyNotFound for dictDiff
	diff = {}
	# Check all keys in first dict
	for key in first.keys():
		if (not second.has_key(key)):
			diff[key] = (first[key], KEYNOTFOUND)
		elif (first[key] != second[key]):
			diff[key] = (first[key], second[key])
	# Check all keys in second dict to find missing
	for key in second.keys():
		if (not first.has_key(key)):
			diff[key] = (KEYNOTFOUND, second[key])
	return diff


def mergeDicts(listDicts): return dict(itertools.chain.from_iterable([x.items() for x in listDicts])) #presumes no overlap in keys

def getCluster():
	hostname = os.environ['HOSTNAME'].lower()
	if 		'sh' 	in hostname: return 'sherlock'
	elif   'gpu-15' in hostname: return 'sherlock'
	elif 	'su' 	in hostname: return 'suncat' #important to distinguish suncat2 and 3?
	elif 	'kris' 	in hostname: return 'kris'
	else: raise ValueError, "getCluster did not detect SH or SU in %s"%hostname

def cell2param(cell):
	"""ANGLES ARE IN RADIANS"""
	def angle(v1,v2): 
		return np.arccos(np.dot(v1,np.transpose(v2))/(np.linalg.norm(v1)*np.linalg.norm(v2)))
	a = np.linalg.norm(cell[0])
	b = np.linalg.norm(cell[1])
	c = np.linalg.norm(cell[2])
	alpha = angle(cell[1],cell[2])
	beta  = angle(cell[0],cell[2])
	gamma = angle(cell[0],cell[1])
	return (a,b,c,alpha,beta,gamma)

def partition(domain, relation):
	known = set()
	classes = {}
	for x in domain:
		for representative, values in classes.items():
			if relation(x, representative):
				values.add(x)
				break
		else:
			classes[x] = set([x])
	return classes.values()


def invertDict(d):
	d2 = {}
	print ' d ',d
	for k,v in d.items():
		if v in d2.keys(): d2[v].append(k)
		else: d2[v] = [k]
	print ' \n\n\nd2 ',d2
	return d2

nonmetalSymbs 	= ['H','C','N','P','O','S','Se','F','Cl','Br','I','At','He','Ne','Ar','Kr','Xe','Rn']
nonmetals 		= [ase.data.chemical_symbols.index(x) for x in nonmetalSymbs]
magElems 		= ['Fe','Mn','Cr','Co','Ni']


gasRefs = {'H':'H2','C':'CO','O':'H2O','N':'N2','F':'F2','Cl':'Cl2','Br':'Br2'}
refNames = 	mergeDicts([gasRefs,
			{'Li': u'Li-bcc','Be': u'Be-hcp','C': u'C-diamond'
			,'Na': u'Na-bcc','Mg': u'Mg-hcp','Al':u'Al-fcc','Si': u'Si-diamond'
			,'K': u'K-bcc','Ca': u'Ca-fcc','Sc': u'Sc-hcp','Ti': u'Ti-hcp','Fe': u'Fe-bcc','Co': u'Co-hcp','Ni': u'Ni-fcc','Cu': u'Cu-fcc','Zn': u'Zn-hcp','Ge': u'Ge-diamond'
			,'Rb': u'Rb-bcc','Sr': u'Sr-fcc','Zr': u'Zr-hcp','Nb': u'Nb-bcc','Mo': u'Mo-bcc','Ru': u'Ru-hcp','Rh': u'Rh-fcc','Pd': u'Pd-fcc','Ag': u'Ag-fcc','Cd': u'Cd-hcp','Sn': u'Sn-diamond'
			,'Ba': u'Ba-bcc','Os': u'Os-hcp','Ir': u'Ir-fcc','Pt':u'Pt-fcc','Au': u'Au-fcc','Pb': u'Pb-fcc'}])

def symbols2electrons(symbols,psp='gbrv'):
	symbdict = {'gbrv':{'Ag':19,'Al':3,'As':5,'Au':11,'Ba':10,'Br':7,'B':3,'Be':4,'Ca':10,'Cd':12,'Co':17,'Cr':14,'Cs':9,'C':4,'Cu':19,'Fe':16,'F':7,'Ga':19,'Ge':14,'Hf':12,'Hg':12,'H':1,'In':13,'Ir':15,'I':7,'K':9,'La':11,'Li':3,'Mg':10,'Mn':15,'Mo':14,'Na':9,'Nb':13,'Ni':18,'N':5,'Os':16,'O':6,'Pb':14,'Pd':16,'Pt':16,'P':5,'Rb':9,'Re':15,'Rh':15,'Ru':16,'Sb':15,'Sc':11,'Se':6,'Si':4,'Sn':14,'Sr':10,'S':6,'Ta':13,'Tc':15,'Te':6,'Ti':2,'Tl':13,'V':13,'W':14,'Y':11,'Zn':20,'Zr':12}}
	# added Be, Cd on my own
	return sum([symbdict[psp][x] for x in symbols])

def makeEquivRelation(pairList,domain):
	""" Take any arbitrary relation and make the minimal equivalence relation that contains it"""
	def applyRule1(plist,domain): return set([(x,x) for x in domain]) - set(plist) #new elems
	def applyRule2(plist): 	return set([(y,x) for x,y in plist]) - set(plist) # new elems
	def applyRule3(plist):	return set([(x,z) for x,y in plist for y_,z in plist if y == y_]) - set(plist) 
	pairList += list(applyRule1(pairList,domain))
	pairList += list(applyRule2(pairList))
	n = 1 #set above 0 to enter while loop
	while n > 0: #were any new elements added by rule 3?
		pairList += applyRule3(pairList)
		n = len(applyRule3(pairList))
	return pairList 
# print partition(xrange(-3, 5), lambda x,y: (x-y) % 4 == 0)
def allElems(listoflists):
	return list(set(itertools.chain(*listoflists))) 

def gcd(args):
	"""Greatest common denominator of a list"""
	if len(args) == 1:
		return args[0]
	L = list(args)
	while len(L) > 1:
		a = L[len(L) - 2]
		b = L[len(L) - 1]
		L = L[:len(L) - 2]
		while a:
			a, b = b%a, a
		L.append(b)
	return abs(b)

def normalizeList(l):
	""" [a,a,a,a,b,b,c,c] => [a,a,b,c] """
	if len(l)==0: return l
	d 	= {x:l.count(x) for x in l}
	div = gcd(d.values())
	norm= [[k]*(v/div) for k,v in d.items()]
	return [item for sublist in norm for item in sublist]

"""Update This if generalSummary changes in jobs.py"""
initcols =  (['jobkind','inittraj_pckl','comments','trajcomments'
			,'name','relaxed','kind','pw','xc','kptden','psp','dwrat','econv','mixing'
			,'nmix','maxstep','nbands','sigma','fmax','dftcode']
			+ ['structure','bulkvacancy_json','bulkscale_json','bulkparent','sites_base64'
			, 'facet_json','xy_json','layers','constrained','symmetric','vacuum']
			+ ['adsorbates_json','vibids_json','xtol','strain','delta']
			+ ['params_json','params_json','natoms','symbols_pckl','symbols_str'
			,'comp_pckl','comp_str','species_pckl','species_str','numbers_pckl','numbers_str'
			,'metalstoich_pckl','metalstoich_str','metalspecies_str','metalcomp_str'
			,'metalspecies_pckl','metalcomp_pckl','paramsinit_pckl','kpt_pckl','kpt_str','emt','emtsym'
			,'spinpol','dipole','calclabel','strjob']
			+ ['bulkvacancy_str','bulkscale_str','facet_str','xy_str'
				,'surfaceArea','adsorbates_str']
			+['fwid','launchdir'])

def pointGroupToSymnum(pg):
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


		
###