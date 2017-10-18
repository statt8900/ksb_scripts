#!/usr/bin/env python
#above line selects special python interpreter needed to run espresso
#SBATCH -p iric,owners,normal
#SBATCH -x gpu-14-1,sh-20-35
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
#SBATCH --time=10:05:00
#################
#number of nodes you are requesting
#SBATCH --nodes=1
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

#######################################################################################
#######################################################################################
#Imports
#######################################################################################
import numpy as np
import os,sys
from os import listdir
from os.path import isfile, join
import shutil
import glob
from numpy import sqrt
from numpy import pi

from ase import io
from ase import Atom, Atoms
from ase import optimize
from ase.data.molecules import molecule
from ase import units
from ase.io.trajectory import PickleTrajectory
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo
from ase.constraints import FixAtoms
from ase.visualize import view

from ase.lattice.surface import *
from ase.constraints import FixAtoms
from ase.io import read, write

from espresso import espresso
#######################################################################################
#######################################################################################
#Inputs
#######################################################################################
#Important Stuff
xc          = 'BEEF'
kpt         = (4,4,1)
pw_cutoff   = 500   # eV
dw_cutoff   = 5000  # eV
spinpol     = False  #Spin-polarization
parflags    = '-npool 1' # Number of nodes to parallalize
#Less-Important Stuff
dipole      = {'status':True}
nbands      = -10
sigma       = 0.1   #Fermi Electron Temperature
magmom      = 2.0   #Initial guess
AFM         = False #Antiferromagnetic
output      = {'avoidio':False,
               'removewf':True,
               'wf_collect':False}
latticeopt  = 'co.py'

#Convergence
convergence = {'mixing': 0.1,
               'maxsteps': 500,
               'energy':5e-4,
               'diag':'david',
               'mixing_mode':'plain'}
fmax       = .05
calcdir    = 'calcdir'


####################################################################################
####################################################################################
#Main Script
###################################################################################

####################################
# File Management
###################################

path   = os.getcwd() + '/'
files  = [f for f in listdir(path) if isfile(join(path, f))] ; files.sort()

if os.path.exists(path+calcdir): shutil.rmtree(path+calcdir)

if os.path.exists(latticeopt):
    for f in files:
        if f[-5:]=='.traj' or f[-4:] == ".cif": 
            traj=f
            initial=io.read(traj)
else:
    inputs = glob.glob('input_*')

    jobid  = max([int(x.split('.')[1]) for x in glob.glob('nodefile.*')]) 

    if len(inputs)==0: n = 0
    else: 
        n = 1 + max([int(x.split('_')[1]) for x in inputs]) #current iteration number

    os.makedirs('input_'+str(n))
    inputdir = path + 'input_' + str(n)+'/'

    for f in files:
        if f[-5:]=='.traj': 
            traj=f
            initial=io.read(traj)
            shutil.copy(path+f,inputdir)
        if f=='RELAX.py': 
            shutil.copy(path+f,inputdir)
        if f[:8]=='nodefile' or f[:12]=='uniqnodefile':
            if int(f.split('.')[1]) != jobid: os.remove(f)
        if f=='qn.log': os.remove(f)

############################
#Initialize Magnetic Moments
############################
TM = ['Fe','Co','Ni'] #Naturally ferromagnetic elements

magmoms = [0]*len(initial.get_positions())

if spinpol==True:
    mag_index = []
    for atom in initial:
        for m in TM:
            if atom.symbol == m:
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

####################################
#pbc=(1,1,1)
#initial.set_pbc(pbc)

###################################
#Create Coarse and Fine Calculators
###################################
calc = espresso(
        pw          = .8*pw_cutoff,
        dw          = .8*dw_cutoff,
        xc          = xc,
        kpts        = kpt,
        nbands      = nbands,
        sigma       = sigma,
        spinpol     = spinpol,
        dipole      = dipole,
        convergence = convergence,
        outdir      = calcdir,
        output      = output,
        parflags    = parflags)

calc2 = espresso(
        pw          = pw_cutoff,
        dw          = dw_cutoff,
        xc          = xc,
        kpts        = kpt,
        nbands      = nbands,
        sigma       = sigma,
        spinpol     = spinpol,
        dipole      = dipole,
        convergence = convergence,
        outdir      = calcdir,
        output      = output,
        parflags    = parflags)


############################
# Coarse Optimization
############################
initial.set_calculator(calc)

Traj = PickleTrajectory(traj,'a',initial)

initial.rattle(stdev= 0.01)

print 'pw = %d, fmax = %f'%(.8*pw_cutoff,2*fmax)
dyn = optimize.BFGS(initial, logfile='qn.log', trajectory=Traj)
dyn.run(fmax=2*fmax)

print 'pw %d to f %f'%(pw_cutoff,fmax)


##############################
# Fine Optimization
##############################

initial.set_initial_magnetic_moments(magmoms)
initial.set_calculator(calc2)
dyn = optimize.BFGS(initial, logfile='qn.log', trajectory=Traj)
dyn.run(fmax=fmax)

######################
# If job is successful
######################

with open('converged.txt','w') as f: f.write(str(initial.get_potential_energy()))

if not os.path.exists(latticeopt):
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))] ; files.sort()
    for f in files:
        if str(f)[:8]=='nodefile' or str(f)[:12]=='uniqnodefile' or str(f)=='qn.log':
            if os.path.exists(path+calcdir+'/'+f): os.remove(path+calcdir+'/'+f)

            shutil.move(path+f,path + calcdir)

    calcdirs = [f for f in os.listdir(path+calcdir) if os.path.isdir(os.path.join(path+calcdir, f))]
    for f in calcdirs:
        if str(f)[:2]=='qe': shutil.rmdir(path+calcdir+'/'+f)