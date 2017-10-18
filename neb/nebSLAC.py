#!/usr/bin/env python

#LSF -q suncat2
#LSF -n 60
#LSF -W 50:00
##LSF -q psnehidleq
##LSF -n 1
#LSF -o opt.log
#LSF -e err.log
#LSF -sp 100
#LSF -N

from ase import io
import os
import numpy as np
import glob
from ase.optimize import FIRE
from ase.neb import NEB
from espresso.multiespresso import multiespresso

xc = 'BEEF'

# neb options
climb=False
k=[0.1,0.2,0.2,0.2,0.2,0.1] #remember to tighten spring constants as you go

kpts = (8,8,8)
pw_cutoff=800.
dw_cutoff=8000.
spinpol= False
dipole = {'status':False}
output = {'avoidio':True,'removewf':True,'wf_collect':False}
convergence = {'energy':1e-5,'mixing':0.1,'mixing_mode':'plain','maxsteps':300}
nbands = -10
smearing = 'mv'
#these trajectories have to contain total energies and forces besides the coordinates
initial = io.read('initial.traj')
final = io.read('final.traj')

#determine appropriate kpts and nbands (CHECK with unit cell sizes
TM = ['Ir','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Pt','Ag','Au','Rh','Pd']
numTM = 0
for atom in initial:
	for m in TM:
		if atom.symbol == m:
			numTM = numTM + 1
print 'number of transition metals'
print str(numTM)

"""
kpts_dict = {18: [4,6,1],  #3x2 cell
	27:[4,4,1], #3x3 cell
	36:[4,4,1], #3x4 cell
	54:[4,2,1], #3x6 cell
	72:[2,3,1]} #6x4 cell
kpts = kpts_dict[numTM]

nbands = {18:-56,
        27:-78,
        36:-92,
        54:-108,
        72:-142}
nbands = nbands[numTM]
"""
# if this is the first neb, and the neb traj's are not available, the default images are created; otherwise counts the images available:
nebfiles = glob.glob('neb*.traj')
intermediate_images = len(nebfiles)
print 'number of images read'
print str(intermediate_images)


# check for restart
for j in range(1,intermediate_images+1):
   if os.path.exists('neb%d.traj' % j) and os.path.getsize('neb%d.traj' % j)!=0:
     restart=True
     print '  neb%d.traj not empty' %j
   else:
     restart=False
     print '  neb%d.traj empty' %j
if restart==True:
   print 'restarting:'
else:
   print 'initial run:'


if climb==True:
   print '  climb mode on'
else:
   print '  no climbing'


#set up multiple espresso calculators
#example shown here: slabs: dipole correction is turned on
m = multiespresso(ncalc=intermediate_images,outdirprefix='ekspressen',
		pw=pw_cutoff,
		dw=dw_cutoff,
		xc=xc,kpts=kpts,
		spinpol=spinpol,
		smearing=smearing,
		nbands=nbands,
		dipole=dipole,
		output = output,
		convergence=convergence
)


mag = [0]*len(initial.get_positions())
initial.set_initial_magnetic_moments(mag)
final.set_initial_magnetic_moments(mag)

images = [initial]
if not restart:
  for i in range(intermediate_images):
     images.append(initial.copy())

if restart:
   for j in range(1,intermediate_images+1):
     nebimage=io.read('neb%d.traj' % j)
     nebimage.set_initial_magnetic_moments(mag)
     images.append(nebimage)

images.append(final)

if not climb:
  neb = NEB(images,k=k)
if climb:
  neb = NEB(images,climb=True,k=k)

m.set_neb(neb)

#only for first step, linear interpolation here.
if not restart:
  neb.interpolate()

if not climb:
   qn = FIRE(neb, logfile='qn.log')
if climb:
   qn = FIRE(neb, logfile='qn.log')

for j in range(1,intermediate_images+1):
    traj = io.PickleTrajectory('neb%d.traj' % j, 'a', images[j])
    qn.attach(traj)

qn.run(fmax=0.05)