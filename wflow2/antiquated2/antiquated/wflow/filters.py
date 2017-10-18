from ase.db 	import connect
import os

"""
Filters are cruder versions of Constraints, insofar as they are limited to independently limiting the domain of some SINGLE parameter.

Their advantages are twofold:
	- They require no interface with ase.db or data.db
	- Arbitrary functions can be used (rather than being forced to take the form of a sqlQuery)

They are used (currently) in prefiltering the domains for BulkJob,SurfaceJob, and Surface creation on the fly.

This allows much fewer tentative Job objects to have to be created and tested with Constraints.
"""

asedb = connect('/scratch/users/ksb/db/ase.db')

############################
# Auxillary Filter Functions
############################
def true(y):   return True    #trivially true function
def eq(y):     return lambda x: x == y
def neq(y):    return lambda x: x != y
def LT(y):     return lambda x: x < y
def GT(y):     return lambda x: x > y
def fAND(x,y): return lambda z: x(z) and y(z)
def fOR(x,y):  return lambda z: x(z) or y(z)
def fIN(z):    return lambda x: x in z


def andFilters(fList): 															# Merge two filter dictionaries, combining with logical AND
	def merge(A, B, f):  														# Combine two dictionaries, applying the rule 'f' when the two share a key
		merged = {k: A.get(k, B.get(k)) for k in A.viewkeys() ^ B.viewkeys()}  	# Start with symmetric difference; keys either in A or B
		merged.update({k: f(A[k], B[k]) for k in A.viewkeys() & B.viewkeys()}) 	# Update with `f()` applied to the intersection
		return merged
	outFilter = fList[0]
	for f in fList[1:]:	outFilter = merge(f,outFilter,fAND)
	return outFilter

#################
# Trivial Filters
#################
defaultBulkJobFilter 	= {x:true for x in ['aseid','xc','precalcxc','pw','kpt','dft','psp']}
defaultBareSurfFilter 	= {x:true for x in ['aseid','facet','xy','layers','constrained','symmetric','vacuum','vacancy']}
defaultAdsSurfFilter 	= {x:true for x in ['aseid','facet','xy','layers','constrained','symmetric','vacuum','vacancy','adsorbate']}
defaultSurfJobFilter 	= {x:true for x in ['aseid']}# {x:true for x in ['surfid']} #...

####################
# Calculator Filters
#####################
def xcFilter(x): 	return {'calc': eq(x)} #specify XC functional
def dftFilter(x): 	return {'dft': eq(x)} #specify dft code
def pwFilterLT(x): 	return {'pw':LT(x)}
def pwFilterGT(x): 	return {'pw':GT(x)}
def kptFilter(x): 	return {'kpt':eq(x)}

preCalcRule = {'precalc':fOR(eq('None'),eq('PBE'))} # ESSENTIAL

#################
# BulkJob Filters
#################


def selectElems(elemList):
	"""ONLY traj's with one of these elements are included"""
	return {'aseid':fIN([item for sublist in [[x.id for x in asedb.select(e)] for e in elemList] for item in sublist])}

def singleElemFunc(ID): 
	l = asedb.get(ID).get('numbers')
	return all([x==l[0] for x in l]) #every element in element list is equal

singleElems = {'aseid':singleElemFunc}



relaxedFilter 	= {'aseid':lambda ID: asedb.get(ID).get('relaxed')}
unrelaxedFilter = {'aseid':lambda ID: not asedb.get(ID).get('relaxed')}

def isBulkFunc(ID): 
	try: return asedb.get(ID).get('kind') 	== 'bulk'
	except KeyError: return False #no longer in database
def isSurfFunc(ID): 
	try: return asedb.get(ID).get('kind') 	== 'surface'
	except KeyError: return False #no longer in database

isBulkFilter = {'aseid':isBulkFunc} 			# ESSENTIAL
isSurfFilter = {'aseid':isSurfFunc} 			# ESSENTIAL

essentialBulkFilters = 	[defaultBulkJobFilter
						,isBulkFilter
						,preCalcRule
						,unrelaxedFilter]   # change with caution

######################
# Bare Surface Filters
######################
def facetFilter(x): 		return {'facet':	eq(x)}
def xyFilter(x): 			return {'xy':		eq(x)}
def layerFilter(x): 		return {'layers':	eq(x)}
def constrainedFilter(x): 	return {'constrained':	eq(x)}
def slabSymFilter(x): 		return {'symmetric':eq(x)}
def vacuumFilter(x): 		return {'vacuum':	eq(x)}
def vacancyFilter(x): 		return {'vacancy':	eq(x)}
def adsorbateFilter(x): 	return {'adsorbate':eq(x)}

hasAdsorbateFilter 	= {'adsorbate':neq({})}
hydrogenFilter 		= {'adsorbate':lambda x: len(x)==1 and x.keys()==['H']}

essentialBareSurfFilters 	= 	[defaultBareSurfFilter
								,isBulkFilter
								,relaxedFilter]

essentialSurfJobFilters 	= [defaultSurfJobFilter
								,isSurfFilter]	