from objectClass import Object

"""
The XC object. Currently just a string; could later incorporate internal structure
"""

class XC(Object):
	def __init__(self,name,coeffs=[[]]):
		self.name = name     # String
		self.coeffs = coeffs

	#################################################################################
	#################################################################################
	#################################################################################
	def sqlTable(self):  return 'xc'
	def sqlCols(self):   return ['name','coeffs']
	def sqlInsert(self): return [self.name,str(self.coeffs)]
	def sqlEq(self):     return ['name']

