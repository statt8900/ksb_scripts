#!/usr/bin/env python
#SBATCH -p iric
#################
#set a job name
#SBATCH --job-name=myjob
#################
#a file for job output, you can check job progress
#SBATCH --output=myjob.out
#################
# a file for errors from the job
#SBATCH --error=myjob.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=50:05:00
#################
#number of nodes you are requesting
#SBATCH --nodes=5
#################
#SBATCH --mem-per-cpu=4000
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END,FAIL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=ksb@stanford.edu
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks-per-node=16
#################

from ase import io
import os
import numpy as np
import glob
from ase.optimize import FIRE
from ase.neb import NEB
from espresso.multiespresso import multiespresso

#####################################
# INPUTS
#####################################
# neb options

spring = [0.3]  #remember to tighten spring constants as you go
images = 5
climb=False

#calculator settings
xc          = 'BEEF-vdW'
kpts        = (4,4,1)
pw_cutoff   = 500   # eV
dw_cutoff   = 5000  # eV
spinpol     = True  #Spin-polarization

#Less-Important Stuff
dipole      = {'status':False}
nbands      = -10
sigma       = 0.1   #Fermi Electron Temperature
magmom      = 2.0   #Initial guess
AFM         = False
output      = {'avoidio':False,
               'removewf':True,
               'wf_collect':False}

#Convergence
smearing='fd'
convergence = {'energy':1e-5,
              'mixing':0.1,
              'nmix':20,
              'mix':4,
              'mixing_mode':'local-TF',
              'maxsteps':500}
fmax       = .05

#######################################################
#####################################

k=spring*(images+1)

#these trajectories have to contain total energies and forces besides the coordinates
initial = io.read('initial.traj')
final = io.read('final.traj')



# if this is the first neb, and the neb traj's are not available, the default images are created; otherwise counts the images available:
if os.path.exists('neb1.traj') and os.path.getsize('neb1.traj')!=0:
	nebfiles = glob.glob('neb*.traj')
	intermediate_images = len(nebfiles)
	print 'number of images read'
	print str(intermediate_images)
else: 
#	print 'number of images interpolated'
#	print str(intermediate_images)
	intermediate_images = 7
	
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
m = multiespresso(ncalc=intermediate_images,outdirprefix='ekspressen',
		pw=pw_cutoff,
		dw=dw_cutoff,
		xc=xc,kpts=kpts,
		spinpol=spinpol,
    sigma=sigma,
		smearing=smearing,
		nbands=nbands,
        #psppath='/nfs/slac/g/suncatfs/sw/external/esp-psp/gbrv1.5pbe',
		dipole=dipole, 
				output = output,
		convergence=convergence
)

############################
#Initialize Magnetic Moments
############################
TM = ['Fe','Co','Ni'] #Naturally ferromagnetic elements

magmoms = [0]*len(initial.get_positions())

if spinpol==True:
    mag_index = []
    for atom in initial:
        for met in TM:
            if atom.symbol == met:
                mag_index.append(atom.index)        

    mag_index_1 = mag_index[::2]
    mag_index_dum = mag_index
    del mag_index_dum[::2]
    mag_index_2 = mag_index_dum
    for j in mag_index_1:     magmoms[j] = magmom
    if AFM:
        for k in mag_index_2: magmoms[k] = -magmom
    else:    # FM case
        for k in mag_index_2: magmoms[k] = magmom

initial.set_initial_magnetic_moments(magmoms)
final.set_initial_magnetic_moments(magmoms)

images = [initial]
if not restart:
  for i in range(intermediate_images):
     images.append(initial.copy())

if restart:
   for j in range(1,intermediate_images+1):
     nebimage=io.read('neb%d.traj' % j)	
     nebimage.set_initial_magnetic_moments(magmoms)
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
    print 'I am here! neb%d.traj'
    qn.attach(traj)

qn.run(fmax=fmax)
