#!/usr/bin/env python

#LSF -n 12
#LSF -o opt.log
#LSF -e err.log
#LSF -N
#LSF -W 30:00 
#LSF -q suncat2
#LSF -sp 90

import numpy as np
from numpy import *
from numpy import pi
from ase import Atoms,Atom
from ase import io
from ase.lattice.spacegroup.cell import cellpar_to_cell
from ase.structure import bulk
from espresso import espresso
from ase.optimize import BFGS
from ase.visualize import view
from scipy.optimize import fmin

lcfix = 3.95333607
lc = [3.95333607]

xc = 'pbesol'
pw_cutoff   = 800
dw_cutoff   = 8000
kpts        = (8,8,8)
spinpol     = False
dipole      = {'status':False}
nbands      = -20
convergence = {'mixing': 0.1,
               'maxsteps': 500,
               'energy':5e-4,
               'diag':'david',
               'mixing_mode':'plain'}

output      = {'avoidio':False,
               'removewf':True,
               'wf_collect':False}
######################################
######################################
######################################
######################################
def get_energy(x):
    global iteration
    iteration += 1

    cellpar = [lcfix,x[0],lcfix, 90,90,90]
    ucell = cellpar_to_cell(cellpar)

    # setting up the atomic configuration
    atoms = io.read('init.traj')
    #Set the unit cell
    atoms.set_cell(ucell,scale_atoms=True)

    atoms.write('input.traj')

    calc = espresso(pw=pw_cutoff,
                    dw=dw_cutoff,
                    xc=xc,
                    kpts =kpts,
                    nbands=nbands,
                    spinpol=spinpol,
                    dipole = dipole,
                    convergence=convergence,
                    output = output,
                    outdir='calcdir'+str(np.round(x[0],5)))

    atoms.set_calculator(calc)
    trajfile = 'out'+str(np.round(x[0],5))+'.traj'
    logfile = 'out'+str(np.round(x[0],5))+'.log'
    dyn = BFGS(atoms, logfile=logfile, trajectory=trajfile)
    dyn.run(fmax=0.01)
    energy = atoms.get_potential_energy()
    print('%20.8f%20.8f' % (x[0],energy))
    f = open('out'+str(np.round(x[0],5))+'.energy', 'w')
    f.write(repr(x) + ' ' + str(energy))
    f.close()
    calc.stop()
    del calc
    return energy

print('%20s%20s' % ('lc[0]', 'lc[1]''energy'))
iteration = 0.

x = fmin(get_energy, x0=lc, xtol=0.001, ftol=0.00001)
print('Best lc')
print x