import json
from objectClass import Object
from lattice     import structDict
from printParse  import underscoreSep,printStoich,printAds,printNp
from bulk import cell2param
class BulkResult(Object):
	def __init__(self,jobid,cell,pos,magmom,eng,bulkmodulus,bFit,xc,t,nIter):
		a,b,c,alpha,beta,gamma = cell2param(cell)
		
		self.jobid         = jobid
		self.a             = a
		self.b             = b
		self.c             = c
		self.alpha         = alpha
		self.beta          = beta
		self.gamma         = gamma
		self.pos           = pos
		self.magmom        = magmom
		self.energy        = eng 
		self.bulkModulus   = bulkmodulus
		self.bFit          = float(bFit)
		self.xcCoeffs      = xc 
		self.time          = t
		self.nIter         = nIter

	def __str__(self): return "BULKRESULT FOR JOB#"+str(self.jobid)
				
	def makeStructure(self):
		from pymatgen import Structure
		return Structure(self.cell,self.job.bulk.atomList(),self.pos,site_properties={'magmom':self.magmom})

	def makeUnitCell(self):
		from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
		#self.makeStructure().to(filename="bulkUnit"+str(self.job)+".cif")
		out =  SpacegroupAnalyzer(self.makeStructure(),symprec=0.1, angle_tolerance=5).get_conventional_standard_structure()
		#out.to(filename="conventional_"+str(self.job)+".cif")
		return out 
	#######################################
	#######################################
	def sqlTable(self):  return 'bulkresult'

	def sqlCols(self):   return ['jobid'
								,'a','b','c','alpha','beta','gamma'
								,'pos'
								,'magmom'
								,'energy'
								,'bulkmodulus'
								,'bfit'
								,'xccoeffs'
								,'time'
								,'niter']

	def sqlInsert(self): 
		return [self.jobid
				,self.a,self.b,self.c,self.alpha,self.beta,self.gamma
				,str(self.pos)
				,str(self.magmom)
				,self.energy
				,self.bulkModulus
				,self.bFit
				,str(self.xcCoeffs)
				,self.time
				,self.nIter]
	def sqlEq(self): return ['jobid']


def jobKind(strJob): return 'surface' if '___' in strJob else 'bulk' #surface jobs have '___' in their string representation

