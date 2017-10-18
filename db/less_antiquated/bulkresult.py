from bulkjob import cell2param
from ase     import Atoms
from db      import queryManyFromOne,getCluster
from ast     import literal_eval
import os

"""
Defines result object of a bulk calculation
"""

if getCluster()=='sherlock':
	from ase.db  import connect

class BulkResult(object):
	def __init__(self,jobid,cell,pos,magmom,eng,forces,bulkmodulus,bfit,xc,t,niter):
		a,b,c,alpha,beta,gamma = cell2param(cell) #WHY WAS IT BOXED 

		self.jobid         = jobid
		self.cell          = cell
		self.a             = a
		self.b             = b
		self.c             = c
		self.alpha         = alpha
		self.beta          = beta
		self.gamma         = gamma
		self.pos           = pos
		self.magmom        = magmom
		self.energy        = eng 
		self.forces        = forces 
		self.bulkModulus   = bulkmodulus
		self.bfit          = float(bfit) # WHY IS IT BOXED
		self.xcCoeffs      = xc 
		self.avgtime       = t
		self.niter         = niter

	def __str__(self): return "BULKRESULT FOR JOB#"+str(self.jobid)

	def getASEid(self):
		db    = connect('/scratch/users/ksb/db/ase.db')
		for row in db.select(relaxed=True,kind='bulk'):
			if self.jobid == row.jobid: 
				return row.id
		raise ValueError, 'Tried to get ase id before it was even added to database?'

	def addToASEdb(self,job): 
		db    = connect('/scratch/users/ksb/db/ase.db')	
		
		info  = {'name':job.bulkname,'structure':job.structure,'comments':'%s\nInitial structure: %d'%(job.bulkcomments,job.bulkid),'emt':job.emt,'relaxed':True,'jobid':job.ID}
		optAtoms = Atoms(symbols=job.symbols,scaled_positions=self.pos,cell=self.cell,magmoms=self.magmom,constraint=job.constraints)
		db.write(optAtoms,key_value_pairs=info)

	#############################################################################
	def sqlTable(self):  return 'bulkresult'

	def sqlCols(self):   return ['jobid'
								,'aseid'
								,'cell'
								,'a','b','c','alpha','beta','gamma'
								,'positions'
								,'magmoms'
								,'energy'
								,'forces'
								,'bulkmodulus'
								,'bfit'
								,'xccoeffs'
								,'time'
								,'niter']

	def sqlInsert(self): 
		return [self.jobid
				,self.getASEid()
				,str(self.cell)
				,self.a,self.b,self.c,self.alpha,self.beta,self.gamma
				,str(self.pos)
				,str(self.magmom)
				,self.energy
				,str(self.forces)
				,self.bulkModulus
				,self.bfit
				,str(self.xcCoeffs)
				,self.avgtime
				,self.niter]

	def sqlEq(self): return ['jobid']

