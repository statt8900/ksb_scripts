#!/usr/bin/env python

#LSF -n 8
#LSF -o opt.log
#LSF -e err.log
#LSF -N
#LSF -W 20:00 -q suncat
#LSF -sp 90
#!/home/vossj/suncat/bin/python
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
#SBATCH --time=01:20:00
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

from ase import *
from ase import io
from ase.units import Bohr
from ase.io import read
from ase.io.trajectory import PickleTrajectory
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.visualize import *
from ase.io.bader import attach_charges
import numpy as np
import os
import pickle
from espresso import espresso
from os import listdir
from os.path import isfile, join

#######################################################################################
#######################################################################################
#Inputs
#######################################################################################
#Important Stuff
xc          = 'BEEF'
kpt         = (8,8,8)
pw_cutoff   = 800   # eV
dw_cutoff   = 8000  # eV
spinpol     = False  #Spin-polarization
parflags    = '-npool 1' # Number of nodes to parallalize
#Less-Important Stuff
dipole      = {'status':False}
nbands      = -80
sigma       = 0.1   #Fermi Electron Temperature
smearing    = 'mv'
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
path   = os.getcwd() + '/'
files  = [f for f in listdir(path) if isfile(join(path, f))] ; files.sort()


for f in files:
    if f[-5:]=='.traj': 
        traj=f
        atoms=io.read(traj)

atoms.set_pbc((True,True,True))




calc = espresso(pw         = pw_cutoff,
                dw         = dw_cutoff,
                xc         = xc,
                kpts       = kpt,
                nbands     = nbands,
                smearing   = smearing,
                sigma      = sigma,
                spinpol    = spinpol,
                dipole     = dipole,
                convergence= convergence,
                outdir     = calcdir,
		)

print "single point calculation of charge density with cutoff"
atoms.set_calculator(calc)

## charge density - need higher cutoffs for accurate charge densities.
e = atoms.get_potential_energy()
print "e:"
print e

origin,cell,density = calc.extract_charge_density()
f=open('density.pickle','w')
pickle.dump(density,f)
f.close()
f2=open('origin.pickle','w')
pickle.dump(origin,f2)
f2.close()
f3=open('cell.pickle','w')
pickle.dump(cell,f3)
f3.close()

## work function
print "single point calculation of charge density with normal 500eV cutoff and full kpts"


wf = calc.get_work_function(pot_filename="pot.xsf", edir=3)
print 'wf: '
print wf
print 'e :'
print e
fw = open('out.WF','w')
fw.write(str(wf))
fw.close()
os.system('rm pot.xsf')

## bader
cd = pickle.load(open('density.pickle','r'))
nx,ny,nz = np.shape(cd)

#cut away periodic image planes
u=nx-1
v=ny-1
w=nz-1
cd2 = np.empty((u,v,w), np.float)
for i in range(u):
    for j in range(v):
        cd2[i][j][:] = cd[i][j][:w]

bohr = 0.52917721092
integral = np.sum(cd2.flat)*atoms.get_volume()/(u*v*w)/bohr**3
print "total electrons"
print integral


io.write('density.cube',atoms,data=cd2)
xyzname = 'qn.xyz' #name.split('.')[0]+'.xyz'
io.write(xyzname, atoms)
os.system('/nfs/slac/g/suncatfs/sw/external/bader-0.27c/bader  density.cube')  # see http://theory.cm.utexas.edu/henkelman/code/bader/ for print options
atoms = io.read('density.cube')
attach_charges(atoms)
os.system('rm Bv*.cube AtIndex.cube')


dos = calc.calc_pdos(nscf=True, kpts=(8,8,1), Emin=-50.0,Emax=50.0, tetrahedra=False, sigma=0.08)

#save dos and pdos into pickle file
import pickle
f = open('dos.pickle', 'w')
pickle.dump(dos, f)
f.close()

dos = calc.calc_pdos(nscf=True, kpts=(8,8,1), Emin=-50.0,Emax=50.0, tetrahedra=False, sigma=0.2)

f = open('dos_smear.pickle', 'w')
pickle.dump(dos, f)
f.close()
