# External Modules
import itertools
#Internal Modules
import misc
import dbase  		as db
import printParse 	as pp

########################################################
########################################################
class Equivalence(object):
	def __init__(self,detail2constraint,detail1constraint='1'):
		self.detail2constraint 	= detail2constraint 			# All jobs that are linked by the relation will be grouped together
		self.detail1constraint 	= detail1constraint 	# 
	def output(self,extraConstraint='1'): 
		dom = db.queryCol('fwid',pp.AND([extraConstraint,self.detail1constraint]))
		pairs = db.query(['id0','id1'],self.detail2constraint,'pairs')
		if self.detail1constraint == '1' and extraConstraint == '1': return pairs
		return [(x,y) for (x,y) in pairs if (x in dom) and (y in dom)]
	
	def classes(self,extraConstraint='1'): 
		output 				= self.output(extraConstraint)
		#print 'output = ',output
		domain 				= misc.allElems(output)
		#outputEquivRelation = misc.makeEquivRelation(output,domain) # output # assume it is already an equivalence relation? #
		return misc.partition(domain,lambda x,y: (x,y) in output)
	def getClass(self,launchdir): 
		return [x for x in self.classes() if launchdir in x][0]	


def equalcalcs(cst='1'): 	return	Equivalence('calcs',cst)
def names(cnstrt='1'): 		return	Equivalence('names',cnstrt)
def structs(cnstrt='1'): 	return	Equivalence('structures',cnstrt)
def metalstoichs(cst='1'): 	return  Equivalence('metalstoich_strs',cst)
def metalspecies(cst='1'): 	return  Equivalence('metalspecies_strs',cst)
def kptdens(cnstrt='1'): 	return	Equivalence('kptdens',cnstrt)
def xcs(cnstrt='1'): 		return	Equivalence('xcs',cnstrt)
def latoptdof(cnstrt='1'): 	return	Equivalence('latoptdof',cnstrt)

########################################################
########################################################
class Detail2a(object):
	def __init__(self,colname,rank,equal,querystrs):
		self.colname 		 = colname
		self.rank 			 = rank
		self.equal 			 = equal 		# if the query is inherently an equivalence relation (which is expensive to enforce)
		self.querystrs 		 = querystrs

	def addToPairs(self): 

		try: db.addCol(self.colname,'integer','pairs')
		except: pass
		pairs = (db.queryTuple(2,self.querystrs))

		if not self.equal:
			reflex,sym,trans = set([]),set([]),set([])
			for p0 in pairs: reflex ^= set([(p0[0],p0[0]),(p0[1],p0[1])])
			pairs ^= reflex
			for p0 in pairs: sym^=set([(p0[1],p0[0])])
			pairs ^= sym
			for (p0,p1) in itertools.product(pairs,pairs):
				if p0[1] == p1[0]:	trans^=set([(p0[0],p1[1])])
			pairs ^= trans
			pairs = list(pairs)
		

		db.sqlexecutemany('UPDATE pairs SET {0} = 1 where id0=? and id1= ?'.format(self.colname),pairs)
		db.sqlexecutemany('INSERT OR IGNORE INTO pairs (id0,id1,{0}) VALUES (?,?,1)'.format(self.colname),pairs)


class Detail2b(object):
	def __init__(self,colname,d2as):
		self.colname 	= colname
		self.d2as		= [d2Dict[x] for x in d2as]
	def rank(self): return 1+max([d.rank for d in self.d2as])
	def equal(self): return all([d.equal for d in self.d2as])
	def const(self): return ' AND '.join([d.colname for d in self.d2as])
	def addToPairs(self): 
		db.sqlexecute('UPDATE pairs SET {0} = 1 where {1}'.format(self.colname,self.const()))



def equalD2(col): return Detail2a(col+'s',0,True,pp.EQUAL(col))
detail2as = [equalD2(x) for x in ['name'
								,'structure'
								,'jobkind'
								,'dof'
								,'pw','xc','kptden','psp','dwrat','econv','fmax'
								,'dftcode'
								,'metalspecies_str'
								,'metalstoich_str']]
detail2as =[Detail2a('kptflex',0,True,"(j0.kptden=j1.kptden OR j0.kind='molecule' OR j1.kind='molecule')")] + detail2as
			
d2Dict = {d2.colname:d2 for d2 in detail2as}

detail2bs = [Detail2b('xc_kptden_name',	['xcs','kptdens','names'])
			,Detail2b('latoptdof',		['jobkinds','dofs'])
			,Detail2b('calcs',['pws','xcs','kptflex','dftcodes','psps','dwrats','econvs','fmaxs'])
			]
detail2s = detail2as + detail2bs

def addD2Names():
	for d in detail2s:
		try: 
			print 'adding colname: ',d.colname
			db.addCol(d.colname,'integer','pairs')
		except: pass


def addPairCols():
	for d2 in detail2as: print 'adding ',d2.colname;d2.addToPairs()
	for d2 in detail2bs: print 'adding ',d2.colname;d2.addToPairs()
	print 'added %d new cols to pairs'%len(detail2s)

