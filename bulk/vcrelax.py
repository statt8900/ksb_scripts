#!/home/vossj/suncat/bin/python
#above line selects special python interpreter needed to run espresso
#SBATCH -p iric
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
#SBATCH --time=21:20:00
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#SBATCH --mem-per-cpu=4000
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=FAIL,END
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=ksb@stanford.edu
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks-per-node=16
#################

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
xc          = 'BEEF'
kpt         = (8,8,8)
pw_cutoff   = 800   # eV
dw_cutoff   = 8000  # eV
spinpol     = True  #Spin-polarization

#Less-Important Stuff
dipole      = {'status':True}
nbands      = -10
sigma       = 0.1   #Fermi Electron Temperature
magmom      = 2.0   #Initial guess
AFM         = False #Antiferromagnetism
output      = {'avoidio':False,
                'removewf':True,
                'wf_collect':False}
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
path   = os.getcwd() + '/'
inputs = glob.glob('input*')
jobid  = max([int(x.split('.')[1]) for x in glob.glob('nodefile.*')]) 
files = [f for f in listdir(path) if isfile(join(path, f))]
files.sort()

if len(inputs)==0: n = 0
else: 
    n = 1 + max([int(x.split('_')[1]) for x in inputs]) #current iteration number

os.makedirs('input_'+str(n))
inputdir = path + 'input_' + str(n)+'/'
calcdir  = path + 'calcdir_'  + str(n)+'/'
oldcalcdir  = path + 'calcdir_'  + str(n-1)+'/'

for f in files:
    if f[-5:]=='.traj': 
        traj=f
        atoms=io.read(traj)
        shutil.copy(path+f,inputdir)
    if f=='vcrelax.py': 
        shutil.copy(path+f,inputdir)
    if f[:8]=='nodefile' or f[:12]=='uniqnodefile':
        if int(f.split('.')[1]) != jobid:
            shutil.move(path + f,oldcalcdir)
    if f=='qn.log':
        shutil.move(path + f,oldcalcdir)

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

    mag_index_1 = mag_index[::2]
    mag_index_dum = mag_index
    del mag_index_dum[::2]
    mag_index_2 = mag_index_dum
    for j in mag_index_1:     magmoms[j] = magmom
    if AFM:
        for k in mag_index_2: magmoms[k] = -magmom
    else:    # FM case
        for k in mag_index_2: magmoms[k] = magmom

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
                output= {'avoidio':False,
                        'removewf':True,
                        'wf_collect':False},
                convergence={'mixing': 0.1,
                            'maxsteps': 500,
                            'energy':5e-4,
                            'diag':'david',
                            'mixing_mode':'plain'}
                )

atoms.set_calculator(calc)
energy = atoms.get_potential_energy() #trigger espresso to be launched
print 'Optimized unit cell:'
print calc.get_final_structure().cell
print 'Residual stress:'
print calc.get_final_stress()
print 'Total energy (eV):'
print energy

io.write(traj,calc.get_final_structure())


######################
# If job is successful
######################

with open('converged.txt','w') as f: f.write(str(energy))

files = [f for f in listdir(path) if isfile(join(path, f))]

files.sort()
for f in files:
    if str(f)[:8]=='nodefile' or str(f)[:12]=='uniqnodefile' or str(f)=='qn.log':   
        shutil.move(path+f,calcdir)