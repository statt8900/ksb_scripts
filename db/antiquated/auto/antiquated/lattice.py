from objectClass import Object
from math import pi,sin,cos,tan,degrees,radians

"""
This module defines the Struct class and its supporting class, Bravais.
Structs are set forms for taking a certain stoichiometry and creating a 3D structure.
The current list is nowhere near exhaustive, but these can be made on the fly as needed.
To make the struct class exhaustive, and to not rely on unique string for every structure, 
	possibly its components could be parameterized by a new class? (GenericStruct?)
"""

class Bravais(Object):
	def __init__(self,name,a,b,c,alpha,beta,gamma):
		self.name = name
		self.a = a
		self.b = b
		self.c = c
		self.alpha = alpha #angle between v1 and v2
		self.beta  = beta  #angle between v2 and v3
		self.gamma = gamma #angle between v3 and v1
		self.params = [a,b,c,alpha,beta,gamma]
		

	def primVec(self): 
		a=self.a;b=self.b;c=self.c;
		alpha=radians(self.alpha);beta=radians(self.beta);gamma=radians(self.gamma)
		v1  = [a,0,0]
		v2  = [b*cos(alpha),b*sin(alpha),0]

		v31 = c*cos(gamma)
		v32 = c*(cos(beta)/sin(alpha) - cos(gamma)/tan(alpha))
		v33 = (c**2 - v31**2 - v32**2)**(0.5)
		v3  = [v31,v32,v33]
		return ([v1,v2,v3]) #round3?

	def __str__(self): return self.name

################################################################
# Bravais Functions
################################################################
def makeCubic(a):                          return Bravais('cubic',a,a,a,90,90,90)
def makeHexagonal(a,c):                    return Bravais('hexagonal',a,a,c,120,90,90)
def makeTetragonal(a,c):                   return Bravais('tetragonal',a,a,c,90,90,90)
def makeTrigonal(a,alpha):                 return Bravais('trigonal',a,a,a,alpha,alpha,alpha)
def makeOrthorhombic(a,b,c):               return Bravais('orthorhombic',a,b,c,90,90,90)
def makeMonoclinic(a,b,c,gamma):           return Bravais('monoclinic',a,b,c,90,90,gamma)
def makeTriclinic(a,b,c,alpha,beta,gamma): return Bravais('triclinic',a,b,c,alpha,beta,gamma)

trajBravaisDict={'cubic':        makeCubic
				,'hexagonal':    makeHexagonal
				,'tetragonal':   makeTetragonal
				,'trigonal':     makeTrigonal
				,'orthorhombic': makeOrthorhombic
				,'monoclinic':   makeMonoclinic
				,'triclinic':    makeTriclinic}

#########################################################################################################

class Struct(object):
	def __init__(self,name,elems,n,mapping,pos,bravais,iGuess):
		self.name    = name      # String
		self.elems   = elems     # Int  # of elements required
		self.n       = n         # Int, # of points in basis
		self.map     = mapping   # [Int], maps elements from a stoich onto the positions in pos
		self.pos     = pos       # [[Int]], basis coordinates
		self.bravais = bravais   # [Double]->Bravais
		self.iGuess  = iGuess    # BulkTraj -> [Double], initial guess for parameters

	def allParams(self,sufficientParams):
		return self.bravais(sufficientParams).params

bcc        = Struct('bcc',       1,2,[0,1],[[0,0,0],[1./2,1./2,1./2]],lambda x: makeCubic(x),         lambda x: [x.guessLat()])
sc         = Struct('sc',        1,1,[0],  [[0,0,0]],                 lambda x: makeCubic(x),         lambda x: [x.guessLat()])
cscl       = Struct('cscl',      2,2,[0,1],[[0,0,0],[1./2,1./2,1./2]],lambda x: makeCubic(x),         lambda x: [x.guessLat()])
fcc        = Struct('fcc',       1,1,[0],  [[0,0,0]],                 lambda x: makeTrigonal(x,60),   lambda x: [x.guessLat()*0.707])
rocksalt   = Struct('rocksalt',  2,2,[0,1],[[0,0,0],[1./2,1./2,1./2]],lambda x: makeTrigonal(x,60),   lambda x: [x.guessLat()*0.707])
zincblende = Struct('zincblende',2,2,[0,1],[[0,0,0],[1./4,1./4,1./4]],lambda x: makeTrigonal(x,60),   lambda x: [x.guessLat()*0.707])
diamond    = Struct('diamond',   1,2,[0,1],[[0,0,0],[1./4,1./4,1./4]],lambda x: makeTrigonal(x,60),   lambda x: [x.guessLat()*0.707])
hcp        = Struct('hcp',       1,2,[0,1],[[0,0,0],[2./3,1./3,1./2]],lambda x,y: makeHexagonal(x,y),   lambda x: [y*x.guessLat() for y in [1,1.8]])

structDict = {'fcc':fcc,'bcc':bcc,'sc':sc,'hcp':hcp,'cscl':cscl,'rocksalt':rocksalt,'zincblende':zincblende,'diamond':diamond}

