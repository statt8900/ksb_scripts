#!/usr/bin/env python

#LSF -q suncat2
#LSF -n 12
#LSF -W 30:00
#LSF -o opt.log
#LSF -e err.log
#LSF -sp 100
#LSF -N

import numpy as np    #vectors, matrices, lin. alg., etc.
import matplotlib
matplotlib.use('Agg') #turn off screen output so we can plot from the cluster
from ase.utils.eos import *  # Equation of state: fit equilibrium latt. const
from ase.units import kJ
from ase.lattice import bulk
from ase import *
from espresso import espresso

#######################################################################################
#######################################################################################
#Inputs
#######################################################################################

#Lattice Opt stuff

metal = 'Pd'
# if you have a metal alloy, specify the second metal as well
metal2 = None
a=3.89     #initial guess for lattice constant
N_steps = 11 # ODD number of steps tested for latt. consts.; assume spacing .01 A

#Important Stuff
xc          = 'pbesol'
kpt         = (11,11,11)
pw_cutoff   = 800   # eV
dw_cutoff   = 8000  # eV
spinpol     = False  #Spin-polarization
parflags    = '-npool 1' # Number of nodes to parallalize
#Less-Important Stuff
dipole      = {'status':False}
nbands      = -10
sigma       = 0.1   #Fermi Electron Temperature
magmom      = 2.0   #Initial guess
AFM         = False #Antiferromagnetic
output      = {'avoidio':False,
               'removewf':True,
               'wf_collect':False}

#Convergence
convergence = {'mixing': 0.1,
               'maxsteps': 500,
               'energy':5e-4,
               'diag':'david',
               'mixing_mode':'plain'}
fmax       = .05
calcdir    = 'calcdir'


#######################################################################################
#######################################################################################
#Main Script
#######################################################################################

strains = np.linspace(-0.01*(N_steps-1)/2,0.01*(N_steps-1)/2,N_steps)

# if Mo then use bcc crystal, otherwise fcc
if metal == 'Mo':
  crystal = 'bcc'
else:
  crystal = 'fcc'

volumes = []  #we'll store unit cell volumes and total energies in these lists
energies = []

#setup up Quantum Espresso calculator
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

for i in strains: #loop over scaling factors
    #build Pt unit cell
    if metal2:
      atoms = bulk(metal, crystal, a=a+i, cubic=True)
      atoms.set_chemical_symbols(metal+'3'+metal2)
    else:
      atoms = bulk(metal, crystal, a+i)
    atoms.set_pbc((1,1,1))                #periodic boundary conditions about x,y & z
    atoms.set_calculator(calc)            #connect espresso to Pt unit cell
    volumes.append(atoms.get_volume())    #append the current unit cell volume
                                          #to list of volumes
    energy=atoms.get_potential_energy()   #append total energy to list of
    energies.append(energy)               #energies

eos = EquationOfState(volumes, energies) #Fit calculated energies at different
v0, e0, B = eos.fit()                    #lattice constants to an
                                         #equation of state


# setup bulk using optimized lattice and save it

if metal2:
  best_a = (v0)**(1./3.) # Angstroms
  atoms = bulk(metal, crystal, a=best_a, cubic=True)
  atoms.set_chemical_symbols(metal+'3'+metal2)
else:
  best_a = (4.*v0)**(1./3.) # Angstroms
  atoms = bulk(metal, crystal, best_a)
atoms.write('bulk.traj')

#output of lattice constant = cubic root of volume of conventional unit cell
#fcc primitive cell volume = 1/4 * conventional cell volume
print 'Lattice constant:', best_a, 'AA'
print 'Bulk modulus:', B / kJ * 1e24, 'GPa'
print '(Fitted) total energy at equilibrium latt. const.:', e0, 'eV'
eos.plot(atoms.get_name()+'-eos.png')    #create a png plot of eos fit
