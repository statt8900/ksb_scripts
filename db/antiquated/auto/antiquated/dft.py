from objectClass import Object
from calc        import pspDict
from cluster      import getCluster

"""
This module defines the class of "DFTcode".
Furthermore, functions to create a GPAW or QE calculator are found in this module
Lastly, a dictionary mapping strings to DFTcode objects is found here
"""


class DFTcode(Object):
	def __init__(self,name,calcFunc,metarunLine): 
		self.name=name                 # String
		self.calcFunc=calcFunc         # Job -> ASE Calculator 

dftGPAW            = DFTcode('gpaw',           gpawCalc)
dftQuantumEspresso = DFTcode('quantumEspresso',qeCalc)


dftDict = {'quantumEspresso':dftQuantumEspresso,'gpaw':dftGPAW}
