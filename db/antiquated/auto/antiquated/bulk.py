import warnings
warnings.filterwarnings("ignore", message="Moved to ase.neighborlist")

import numpy as np
import shutil
from math         import degrees
from objectClass  import Object
from ase          import io
from emt          import EMT
from copy         import deepcopy
from db           import query1all
magElems = ['Fe','Mn','Cr','Co','Ni']

class BulkTraj(Object):
	def __init__(self,id,pth,name=None,stoich=None,bravais=None,
				inita=None,initb=None,initc=None,initalpha=None,
				initbeta=None,initgamma=None,initpos=None,mag=None,
				emt=None,comments=None): 
		if id is None:
			# Creating object for first time
			traj = io.read(pth)
			traj.set_calculator(EMT())
			cell = traj.get_cell()

			atomList = traj.get_chemical_symbols()
			self.atomList = atomList
			self.initPos  = traj.get_scaled_positions()

			self.pth = pth
			self.name = (pth.split('.traj')[0]).split('traj/')[1]

			self.traj = traj

			try:             self.bravais  = traj.info['bravais']
			except KeyError: self.bravais  = 'triclinic' # default is most general (and expensive) case
			try:             self.mag      = traj.info['mag']
			except KeyError: self.mag      = False       # Default nonmag
			try:             self.comments = traj.info['comments']
			except KeyError: self.comments = None

			self.initEMT = traj.get_potential_energy()

			a,b,c,alpha,beta,gamma=cell2param(cell)
			self.inita = a
			self.initb = b
			self.initc = c
			self.initalpha = alpha
			self.initbeta  = beta
			self.initgamma = gamma
		else:
			self.id=id;self.pth=pth;self.name=name;self.stoich=stoich
			self.bravais=bravais;self.inita=inita;self.initb=initb;
			self.initc=initc;self.initalpha=initalpha;self.initbeta=initbeta
			self.initgamma=initgamma;self.initpos=initpos;self.mag=mag;
			self.emt=emt;self.comments=comments
			traj = io.read(pth)
			atomList = traj.get_chemical_symbols()
			self.atomList = atomList

	def __str__(self): return self.name
	
	def fromParams(self,params): 
		atoms = io.read(self.pth)
		a=params[0]
		if   self.bravais =='cubic':    cell = [a,a,a,90,90,90]
		elif self.bravais =='trigonal': cell = [a,a,a,60,60,60]
		elif self.bravais =='hexagonal':cell=[a,a,params[1],90,90,120]
		elif self.bravais =='triclinic':cell=params
		else: raise NotImplementedError
		atoms.set_cell(cell,scale_atoms=True)
		return atoms
	
	def iGuess(self):
		if self.bravais in ['cubic','trigonal']:     return [self.inita]
		elif self.bravais in ['hexagonal']:  return [self.inita,self.initc]
		elif self.bravais in ['triclinic']: return [self.inita,self.initb,self.initc,self.initalpha,self.initbeta,self.initgamma]
		else: raise ValueError, 'Bad entry in "bravais" field for Atoms object info dictionary: '+self.bravais()

	def writeInit(self,path): shutil.copyfile(self.pth,path+'/'+self.name+'.traj')

	def magInit(self,magmom): [(magmom if (ele in magElems) else 0) for ele in self.atomList]

	#################################################################################
	#################################################################################
	#################################################################################
	def sqlTable(self):  return 'bulk'

	def sqlCols(self):   return ['pth'
								,'name'
								,'stoich'
								,'bravais'
								,'inita','initb','initc'
								,'initalpha','initbeta','initgamma'
								,'initpos'
								,'mag'
								,'emt'
								,'comments']

	def sqlInsert(self): return [self.pth
								,self.name
								,str(self.atomList)
								,self.bravais
								,self.inita,self.initb,self.initc
								,self.initalpha,self.initbeta,self.initgamma
								,str(self.initPos)
								,self.mag
								,self.initEMT
								,self.comments]
	def sqlEq(self): return ['pth']
#################################################################################
#################################################################################
#################################################################################
def db2Bulk(pth): return BulkTraj(*query1all('bulk','pth',pth)) # fix this
def cell2param(cell):
	def angle(v1,v2): 
		return degrees(np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))))
	a = np.linalg.norm(cell[0])
	b = np.linalg.norm(cell[1])
	c = np.linalg.norm(cell[2])
	alpha = angle(cell[1],cell[2])
	beta  = angle(cell[0],cell[2])
	gamma = angle(cell[0],cell[1])
	return (a,b,c,alpha,beta,gamma)
