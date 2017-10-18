#!/usr/bin/env python

#LSF -q suncat2
#LSF -n 12
#LSF -W 20:00
#LSF -o opt.log
#LSF -e err.log
#LSF -sp 90
#LSF -N


import os,sys
from os import listdir
import shutil
import glob

from os.path import isfile, join
import numpy as np
from espresso import espresso
from ase import io
from ase.all import *

#######################################################################################
#######################################################################################
#Inputs
#######################################################################################
#Important Stuff
xc          = 'pbesol'
kpt         = (8,8,8)
pw_cutoff   = 800   # eV
dw_cutoff   = 8000  # eV
spinpol     = False  #Spin-polarization
dipole      = {'status':False}

#Less-Important Stuff
nbands      = -10
sigma       = 0.1   #Fermi Electron Temperature
magmom      = 2.0   #Initial guess

output      = {'avoidio':False,
                'removewf':True,
                'wf_collect':False}
calcdir  =   'cellrelax'

#Convergence
convergence={'mixing': 0.1,
            'maxsteps': 500,
            'energy':5e-4,
            'diag':'david',
            'mixing_mode':'plain'}


####################################################
####################################################################################
####################################################################################
#Main Script
###################################################################################

atoms=io.read('init.traj')

############################
#Initialize Magnetic Moments
############################
TM = ['Fe','Co','Ni'] #Naturally ferromagnetic elements

magmoms = [0]*len(atoms.get_positions())
if spinpol==True:
    mag_index = []
    for atom in atoms:
        for m in TM:
            if atom.symbol == m:
                mag_index.append(atom.index)

    for i in mag_index:       magmoms[i] = magmom

atoms.set_initial_magnetic_moments(magmoms)

###################################
#Create Calculator
###################################
calc = espresso(pw     = pw_cutoff,
                dw     = dw_cutoff,
                xc     = xc,
                kpts   = kpt,
                nbands = nbands,
                sigma  = sigma,
                mode   = 'vc-relax',
                cell_dynamics = 'bfgs',
                opt_algorithm = 'bfgs',
                cell_factor   = 5.,
                spinpol = spinpol,
                outdir=calcdir,
                output= output,
                convergence=convergence
                )

atoms.set_calculator(calc)
energy = atoms.get_potential_energy() #trigger espresso to be launched
print 'Optimized unit cell:'
print calc.get_final_structure().cell
print 'Residual stress:'
print calc.get_final_stress()
print 'Total energy (eV):'
print energy

io.write('final.traj',calc.get_final_structure())


######################
# If job is successful
######################

with open('converged.txt','w') as f:
    f.write('Optimized unit cell: %s \nResidual stress: %s\nTotal Energy (ev): %s ' %(str(calc.get_final_structure().cell),calc.get_final_stress(),energy))