from objectClass import Object
from printParse  import underscoreSep,derscoreSep,printList,parseList
from db   import query1all
"""
This module defines the Calc and PSP objects
Additionally, a function to parse str(Calc) into a Calc is here.
Lastly, a dictionary mapping strings to PSP objects is found here.
"""

class Calc(Object):
	def __init__(self,ID,xc,psp,pw,kpt,eConv,fMax,xtol,mixing,nmix,maxsteps,nbands,sigma,magmom):
		"""Calculator""" 
		if isinstance(kpt,unicode) or isinstance(kpt,str):
			kpt = parseList(kpt) if '-' in kpt else int(kpt)
		else: pass

		self.ID       = ID             # Int
		self.xc       = xc             # String (XC name) 
		self.psp      = psp            # String
		self.pw       = pw             # Int
		self.kpt      = kpt            # (Int,Int,Int) OR Int (kpt density)
		self.eConv    = eConv         # Double
		self.fMax     = fMax          # Double
		self.xtol     = xtol           # Double
		self.mixing   = mixing         # Double
		self.nmix     = nmix           # Int
		self.maxsteps = maxsteps       # Int 
		self.nbands   = nbands         # Int, gets multiplied by -1
		self.sigma    = sigma          # Double
		self.magmom   = magmom         # Int, initial guess for magnetic moment 
	def __str__(self): return underscoreSep(1,	[self.xc
												,str(self.pw)
												,printList(self.kpt) if isinstance(self.kpt,tuple) else str(self.kpt)
												,str(self.eConv)
												,str(self.mixing)
												,str(self.nmix)
												,str(self.maxsteps)
												,str(self.nbands)
												,str(self.sigma)
												,str(self.xtol)
												,str(self.magmom)
												,self.psp])


	def sqlTable(self):  return 'calc'
	def sqlCols(self):   return ['xc','psp','pw','kpt','econv','fmax','xtol','mixing','nmix','maxsteps','nbands','sigma','magmom']
	def sqlInsert(self): return [self.xc,self.psp,self.pw,printList(self.kpt) if isinstance(self.kpt,tuple) else str(self.kpt)
								,self.eConv,self.fMax,self.xtol,self.mixing,self.nmix,self.maxsteps
								,self.nbands,self.sigma,self.magmom]
	def sqlEq(self):     return self.sqlCols()



def db2Calc(calcID): return Calc(*query1all('calc','id',calcID))

	