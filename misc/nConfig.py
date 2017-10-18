from espresso import espresso
from ase import Atoms
from emt import EMT
import numpy as np
import itertools
"""
positions = np.array([[0,0,0],[0,0,0.5],[0,0.5,0],[0.5,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0],[0.5,0.5,0.5]])
cell      = np.array([[0,2,2],[2,0,2],[2,2,0]])
elements  = ['Pd','Au']

configs = itertools.product(elements,repeat=8)

eDict={}
for i in range(9): eDict[i]=[]

for i,c in enumerate(list(configs)):
	print i
	a   = Atoms(c,positions=positions,cell=cell,pbc=1,calculator=EMT())	
	e   = a.get_potential_energy()
	nPd = c.count('Pd')
	eDict[nPd].append(e)
"""


positions = np.array([[0,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
cell      = np.array([[0,2,2],[2,0,2],[2,2,0]])
elements  = ['Pd','Au']

configs = itertools.product(elements,repeat=4)

eDict={}
for i in range(5): eDict[i]=[]

for i,c in enumerate(list(configs)):
	print i
	a   = Atoms(c,positions=positions,cell=cell,pbc=1,calculator=EMT())	
	e   = a.get_potential_energy()
	nPd = c.count('Pd')
	eDict[nPd].append(e)

print eDict

# plt.hist(eDict[3],nbins=256)