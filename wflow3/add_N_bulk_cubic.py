#!/usr/bin/env python
import sys
from ase.all import *
from ase.lattice import bulk
from ase.visualize import view
import numpy as np
import os, fnmatch
from Global_Vars import Materials
def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return mathd.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def hcpcoord(x,y,a0):
  return (a0*(x+y*np.cos(np.pi/3)),a0*y*np.cos(np.pi/6))

# sub = '/scratch/users/mstatt/'
path = os.getcwd()
# if len(sys.argv)>1:
# 	symbol = sys.argv[1]
# else:
# 	symbol = path[len(sub):path.find('/',len(sub))]
# if len(sys.argv)>2:
# 	struc = sys.argv[2]
# else:
# 	struc = Materials[symbol]

#atoms = read('/scratch/users/mstatt/'+path[len(sub):path.find('/',len(sub))]+'/lat_opt/out.traj')
atoms = read('init.traj')
symbol = atoms[0].symbol
struc = Materials[symbol]
path = os.getcwd()
struc = 'hcp'
if not os.path.exists(path+'/Empty'):
        os.makedirs(path+'/Empty')
if not os.path.exists(path+'/Osite1'):
	os.makedirs(path+'/Osite1')
if not os.path.exists(path+'/Tsite1'):
	os.makedirs(path+'/Tsite1')
if not os.path.exists(path+'/Osite2'):
        os.makedirs(path+'/Osite2')
if not os.path.exists(path+'/Tsite2') and not struc == 'fcc':
        os.makedirs(path+'/Tsite2')


scale = (3,3,3)
if struc == 'bcc':

	#BCC Lattice
	cell = atoms.get_cell()
	a0 = cell[0,0]

	atoms *= scale
	write(path+'/Empty/init.traj',atoms)
	#Tetrahedral Interstitial Site
	atoms.append(Atom('N',position = (a0/2,a0/4,2 * a0)))
	write(path+'/Tsite1/init.traj',atoms)
	atoms[-1].position[1] = 0
	atoms[-1].position[2] -= a0/4
	write(path+'/Tsite2/init.traj',atoms)
	del atoms[-1]
	#Octahedral Interstitial Site
	atoms.append(Atom('N',position = (a0/2,0,2 * a0)))
	write(path+'/Osite1/init.traj',atoms)
	atoms[-1].position[2] -= a0/2
	write(path+'/Osite2/init.traj',atoms)
elif struc == 'fcc':

	#FCC Lattice
	cell = atoms.get_cell()
	a0 = cell[0,0]
	atoms *= scale
	write(path+'/Empty/init.traj',atoms)
	#Tetrahedral Site
	atoms.append(Atom('N',position = (a0/4,a0/4,a0/4)))
	write(path+'/Tsite1/init.traj',atoms)

	del atoms[-1]
	#Octahedral Site
	atoms.append(Atom('N',position = (a0/2,a0/2,a0/2)))
	write(path+'/Osite1/init.traj',atoms)

	atoms[-1].position[0] += a0/2
	atoms[-1].position[1] -= a0/2
	write(path+'/Osite2/init.traj',atoms)
elif struc == 'hcp':
	if not os.path.exists(path+'/Osite3'):
		os.makedirs(path+'/Osite3')

	#HCP Lattice
	cell = atoms.get_cell()
	a0 = cell[0,0]
	c0 = cell[2,2]
	pos = atoms.get_positions()[1]
	pos[2] *= 0.
	pos[2] += c0
	pos_2 = np.array([a0 * 1/2, a0 * 1/3 * np.cos(np.pi/6), c0/4+c0])
	atoms *= scale

	write(path+'/Empty/init.traj',atoms)
	#Tetrahedral Site
	atoms.append(Atom('N',position = pos))
	write(path+'/Tsite1/init.traj',atoms)
	atoms[-1].position=[0.0,0.0,c0 * 1/2]
	write(path+'/Tsite2/init.traj',atoms)
	del atoms[-1]
	#Octahedral Site
	atoms.append(Atom('N',position = pos_2))
	write(path+'/Osite1/init.traj',atoms)
	atoms[-1].position[0] += a0
	write(path+'/Osite2/init.traj',atoms)
	atoms[-1].position[0] += -a0
	atoms[-1].position[2] += -c0/2
	write(path+'/Osite3/init.traj',atoms)
