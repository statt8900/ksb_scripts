#!/usr/bin/env python



import sys
from ase import Atom, Atoms
import math
import os
import numpy as np
from numpy import sqrt
from numpy import pi
from numpy import arctan as atan


from ase.lattice.surface import *
from ase.build import molecule

from ase.constraints import FixAtoms
from ase.io import read, write

try: traj = sys.argv[1]
except: traj = 'init.traj'
try: script = sys.argv[2]
except: script = 'qn_job.py'

trajname = traj.split('/')[-1]

s_CO = 1.43 # A, C-O bond length
s_OH = 0.94 # A, O-H bond length
s_CH = 1.1 # A, C-H bond length
d_CO = 1.217 # A, C=O bond length

H = Atoms([Atom('H',[0.0,0.0,0.0],tag=1)])

OH = Atoms([Atom('O',[0.0,0.0,0.0],tag=1),
             Atom('H',[0.0,s_OH**0.5,s_OH**0.5],tag=1),
            ])


CO = Atoms([Atom('C',[0.0,0.0,0.0],tag=1),
             Atom('O',[0.0,0.0,d_CO],tag=1),
            ])
#CO.rotate('y',pi/2.)

COOH  = Atoms([Atom('C',[0.0,0.0,0.0],tag=1),
             Atom('O',[d_CO/sqrt(2),0.0,d_CO/sqrt(2)],tag=1),
             Atom('O',[-s_CO/sqrt(2),0.0,s_CO/sqrt(2)],tag=1),
             Atom('H',[-s_CO/sqrt(2),0.0,s_CO/sqrt(2)+s_OH],tag=1),
            ])
#COOH.rotate('y',pi/4.)
#COOH.rotate('z',pi/4.)


CHO = Atoms([Atom('C',[0.0,0.0,0.0],tag=1),
             Atom('H',[-s_CH/sqrt(2),0.0,s_CH/sqrt(2)],tag=1),
			Atom('O',[d_CO/sqrt(2),0.0,d_CO/sqrt(2)],tag=1),
            ])
#CHO.rotate('y',pi/4.)
#CHO.rotate('z',pi/4.)

CH2 = Atoms([Atom('C',[0.0,0.0,0.0]),
             Atom('H',[s_CH/sqrt(2),0.0,-s_CH/sqrt(2)]),
	     Atom('H',[-s_CH/sqrt(2),0.0,-s_CH/sqrt(2)])])

#CH2.rotate('z',-pi/4)
CH2.rotate('x', pi)

CH2O = molecule('H2CO')
CH2O.rotate('y',pi/2.)
CH2O.translate(-CH2O[1].position)
#ads.rotate('z',-pi/4)
CH2O.rotate('x', pi)


CH3 = molecule('CH4')
CH3.pop(4)
#ads.rotate('z',-pi/4)
#ads.rotate('y', pi/4)
CH3.rotate('x', 2*pi/3)
CH3.rotate('y', pi)
CH3.rotate('z', pi/3)


CH2OH = molecule('CH4')
CH2OH[1].symbol='B'
CH2OH.rotate('x',atan(CH2OH[1].position[1]/CH2OH[1].position[2]))
CH2OH.rotate('y',-atan(CH2OH[1].position[1]/CH2OH[1].position[2]))
CH2OH.rotate('x',pi)
position = CH2OH[3].position
CH2OH.pop(3)
CH2OH.append(Atom('O',position))
CH2OH.append(Atom('H',position + [0.,0.,s_OH]))
CH2OH.pop(1)
CH2OH.rotate('x',pi)
CH2OH.rotate('y',-pi/4-pi/16)
CH2OH.rotate('y', pi)
CH2OH.translate([0.3,0,0])

CH2CO = molecule('CH4')
CH2CO[1].symbol='B'
CH2CO.rotate('x',atan(CH2CO[1].position[1]/CH2CO[1].position[2]))
CH2CO.rotate('y',-atan(CH2CO[1].position[1]/CH2CO[1].position[2]))
CH2CO.rotate('x',pi)
position = CH2CO[3].position
CH2CO.pop(3)
CH2CO.append(Atom('C',position + [0.2,0.,0.]))
CH2CO.append(Atom('O',position + [0.2,0.,s_CO]))
CH2CO.pop(1)
CH2CO.rotate('x',pi)
CH2CO.rotate('y',-pi/4-pi/16)
CH2CO.rotate('y', pi)
#CH2CO.rotate('z', -pi/4)
CH2CO.translate([1.,0,0])


CHOCO =molecule('CH4')
CHOCO[1].symbol='B'
CHOCO.rotate('x',atan(CHOCO[1].position[1]/CHOCO[1].position[2]))
CHOCO.rotate('y',-atan(CHOCO[1].position[1]/CHOCO[1].position[2]))
CHOCO.rotate('x',pi)
position = CHOCO[3].position
CHOCO.pop(3)
CHOCO.append(Atom('C',position + [0.2,0.,0.]))
CHOCO.append(Atom('O',position + [0.2,0.,s_CO]))
CHOCO.pop(1)
CHOCO.rotate('x',pi)
CHOCO.rotate('y',-pi/4-pi/16)
CHOCO.rotate('y', pi)
CHOCO[1].symbol='O'
CH2CO.translate([2.,0,0])

# define what you want to make here
ads_list = [COOH]    # [CH2, CH3, CH2OH, CH2CO, CH2O, CH2OH]
adsname  = ['COOH'] # ['CH2', 'CH3', 'CH2OH', 'CH2CO', 'CH2O', 'CH2OH']


sites    = [31,21,27,31]
sitename = ['top','hcp','fcc','bridge']

ztrans = [2,3.7,5.7,1.7]  # note that there is no underscore this is the list
ytrans = [0,0,  0,0]
xtrans = [0,0,  0,-1.5]

yrot = [0,0,0,0]
###

def add_adsorbate():
    """Adds the adsorbate to the slab."""
    ads.rattle()
    ads.translate([x_trans,y_trans,z_trans])
    for atom in ads:
        slab.append(atom)



for j in range(len(ads_list)):
	ads = ads_list[j]
	ads_name = adsname[j]  #''.join(ads_list[j].get_chemical_symbols())
	path = ads_name
	if not os.path.exists(path):
			os.makedirs(path)
	for i in range(len(sites)):
		site = sitename[i]
		#site = str(sites[i])
		subpath = ads_name+'/'+site+'/'
		if not os.path.exists(subpath):
			os.makedirs(subpath)
		ads = ads_list[j].copy()  # otherwise it will shift the adsorbate.
		ads.rotate('y',yrot[i])
		slab = read(traj)
		z_trans = slab.positions[sites[i]][2]+ ztrans[i]   # the underscores are the specificcs
		y_trans = slab.positions[sites[i]][1]+ ytrans[i]
		x_trans = slab.positions[sites[i]][0]+ xtrans[i]
# Custom translation or rotation parameters
		add_adsorbate()
		#mask = [a.position[2]<16.802 for a in slab]
		#slab.set_constraint(FixAtoms(mask=mask))
		write(subpath+'init.traj', slab)
		os.system('cp '+script+' '+ subpath)