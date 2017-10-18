#!/usr/bin/env python



import sys
from ase.all import Atom, Atoms
import math
import os
import numpy as np
from numpy import sqrt
from numpy import pi
from numpy import arctan as atan
from copy import deepcopy

from ase.lattice.surface import *
from ase.data.molecules import molecule

from ase.all import view
from ase.io import read, write
from glob import glob
##################################################
try: traj = glob(sys.argv[1])[0]
except: traj = 'init.traj'
try: script = glob(sys.argv[2])[0]
except: script = 'qn_job.py'
##################################################

s_CO = 1.43 # A, C-O bond length
s_OH = 0.94 # A, O-H bond length
s_CH = 1.1 # A, C-H bond length
d_CO = 1.217 # A, C=O bond length
##################################################

H = Atoms([Atom('H',[0.0,0.0,0.0],tag=1)])

OH = Atoms([Atom('O',[0.0,0.0,0.0],tag=1),
             Atom('H',[0.0,s_OH**0.5,s_OH**0.5],tag=1)])

CO = Atoms([Atom('C',[0.0,0.0,0.0],tag=1),
             Atom('O',[0.0,0.0,d_CO],tag=1)])

COOH  = Atoms([Atom('C',[0.0,0.0,0.0],tag=1),
             Atom('O',[d_CO/sqrt(2),0.0,d_CO/sqrt(2)],tag=1),
             Atom('O',[-s_CO/sqrt(2),0.0,s_CO/sqrt(2)],tag=1),
             Atom('H',[-s_CO/sqrt(2),0.0,s_CO/sqrt(2)+s_OH],tag=1)])

CHO = Atoms([Atom('C',[0.0,0.0,0.0],tag=1),
             Atom('H',[-s_CH/sqrt(2),0.0,s_CH/sqrt(2)],tag=1),
			Atom('O',[d_CO/sqrt(2),0.0,d_CO/sqrt(2)],tag=1)])

CH2 = Atoms([Atom('C',[0.0,0.0,0.0]),
             Atom('H',[s_CH/sqrt(2),0.0,-s_CH/sqrt(2)]),
	         Atom('H',[-s_CH/sqrt(2),0.0,-s_CH/sqrt(2)])])
CH2.rotate('x', pi)

CH2O = molecule('H2CO')
CH2O.rotate('y',pi/2.);CH2O.translate(-CH2O[1].get_position());CH2O.rotate('x', pi)

CH3 = molecule('CH4')
CH3.pop(4)
CH3.rotate('x', 2*pi/3);CH3.rotate('y', pi);CH3.rotate('z', pi/3)

CH2OH = molecule('CH4')
CH2OH[1].symbol='B'
CH2OH.rotate('x',atan(CH2OH[1].get_position()[1]/CH2OH[1].get_position()[2]))
CH2OH.rotate('y',-atan(CH2OH[1].get_position()[1]/CH2OH[1].get_position()[2]))
CH2OH.rotate('x',pi)
position = CH2OH[3].get_position()
CH2OH.pop(3)
CH2OH.append(Atom('O',position))
CH2OH.append(Atom('H',position + [0.,0.,s_OH]))
CH2OH.pop(1)
CH2OH.rotate('x',pi);CH2OH.rotate('y',-pi/4-pi/16);CH2OH.rotate('y', pi)
CH2OH.translate([0.3,0,0])

CH2CO = molecule('CH4')
CH2CO[1].symbol='B'
CH2CO.rotate('x',atan(CH2CO[1].get_position()[1]/CH2CO[1].get_position()[2]))
CH2CO.rotate('y',-atan(CH2CO[1].get_position()[1]/CH2CO[1].get_position()[2]))
CH2CO.rotate('x',pi)
position = CH2CO[3].get_position()
CH2CO.pop(3)
CH2CO.append(Atom('C',position + [0.2,0.,0.]))
CH2CO.append(Atom('O',position + [0.2,0.,s_CO]))
CH2CO.pop(1)
CH2CO.rotate('x',pi);CH2CO.rotate('y',-pi/4-pi/16);CH2CO.rotate('y', pi)
CH2CO.translate([1.,0,0])


CHOCO =molecule('CH4')
CHOCO[1].symbol='B'
CHOCO.rotate('x',atan(CHOCO[1].get_position()[1]/CHOCO[1].get_position()[2]))
CHOCO.rotate('y',-atan(CHOCO[1].get_position()[1]/CHOCO[1].get_position()[2]))
CHOCO.rotate('x',pi)
position = CHOCO[3].get_position()
CHOCO.pop(3)
CHOCO.append(Atom('C',position + [0.2,0.,0.]))
CHOCO.append(Atom('O',position + [0.2,0.,s_CO]))
CHOCO.pop(1)
CHOCO.rotate('x',pi)
CHOCO.rotate('y',-pi/4-pi/16)
CHOCO.rotate('y', pi)
CHOCO[1].symbol='O'
CH2CO.translate([2.,0,0])
##################################################

# define what you want to make here
adsorbate = H    # [CH2, CH3, CH2OH, CH2CO, CH2O, CH2OH]
adsname  = 'H' # ['CH2', 'CH3', 'CH2OH', 'CH2CO', 'CH2O', 'CH2OH']


sites  = [5,6,11]
n      = len(sites)
ztrans = [5,5,3]
xtrans = [0]*n ; ytrans = [0]*n; 
yrot   = [0]*n

#####################################################

slab = read(traj)

xs,ys,zs=[],[],[]
for i in range(n):
	xs.append(slab.positions[sites[i]][0]+ xtrans[i])
	ys.append(slab.positions[sites[i]][1]+ ytrans[i])
	zs.append(slab.positions[sites[i]][2]+ ztrans[i])
#####################################################


def add_adsorbate(slab,ads,x,y,z):
	"""Adds the adsorbate to the slab."""
	ads.rattle()
	ads.translate([x,y,z])
	for atom in ads: slab.append(atom)
	return slab
#####################################################

for j in range(n):
	path = adsname+"_"+str(j+1)
	if not os.path.exists(path): os.makedirs(path)

	freshslab = read(traj)	

	for k in range(j+1): 
		freshslab = add_adsorbate(freshslab,deepcopy(adsorbate),xs[k],ys[k],zs[k])

	write(path+'/'+path+'.traj', freshslab)
	os.system('cp '+script+' '+ path)



