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
#SBATCH --time=00:30:00
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
kpt         = (8,8,8)
pw_cutoff   = 800   # eV
dw_cutoff   = 8000  # eV
spinpol     = True  #Spin-polarization
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

calcdir    = 'singlepoint_calcdir'


####################################################################################
####################################################################################
#Main Script
###################################################################################

####################################
# File Management
###################################

path   = os.getcwd() + '/'
files  = [f for f in listdir(path) if isfile(join(path, f))] ; files.sort()

for f in files:
    if f[-5:]=='.traj': 
        traj=f
        atoms=io.read(traj)


calc = espresso(
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


atoms.set_calculator(calc)

energy=atoms.get_potential_energy()

with open('energy.txt','w') as f: f.write(str(atoms.get_potential_energy()))