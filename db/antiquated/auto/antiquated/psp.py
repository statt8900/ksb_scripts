from objectClass import Object

class PSP(Object):
	def __init__(self,name,pth,nElec):
		self.name   = name
		self.pth    = pth # Map Cluster FilePath
		self.nElec  = nElec   # Map ElementSymbol Int  ---- THIS IS ACCURATE FOR GBRV15PBE, not sg15....
	
	#################################################################################
	#################################################################################
	#################################################################################

	def sqlTable(self):  return 'psp'
	def sqlCols(self):   return ['name','pth','nelec']
	def sqlInsert(self): return [self.name,self.pth,str(self.nElec)]
	def sqlEq(self):     return ['pth']