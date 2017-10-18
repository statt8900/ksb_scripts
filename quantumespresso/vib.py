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
#SBATCH --time=00:31:00
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
#task to run per node; each node has 16 s
#SBATCH --ntasks-per-node=16
#################

import sys,os,ase
from ase.optimize import QuasiNewton
from ase.thermochemistry import HarmonicThermo

from espresso import espresso
from espresso.vibespresso import vibespresso
from ase.vibrations import Vibrations

##################
# Inputs #########
##################
pw       = 800                  # planewave cutoff
dw       = 8000                 # density cutoff
kpts     = (8, 8, 8)            # k points
xc       = 'BEEF'               # exchange correlation method
sigma    = 0.1                  # Fermi temoperature
spinpol  = True                 # Spin polarization
dipole   = {'status': False}    # Dipole (surface) correction
output   = {'avoidio':False,
            'removewf':True,
            'wf_collect':False}

vib_atoms   = [4]            # List of atoms allowed to vibrate
temperature = 300            # K

#Convergence
maxsteps=500    
mixing  = 0.1   

##################
# Slab ###########
##################
path = os.getcwd()+'/'
files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
files.sort()
for f in files:
  if str(f)[-5:]=='.traj': traj=str(f)

atoms = ase.io.read(traj)

atoms.rattle()
atoms.set_masses()


calc = vibespresso(
              pw=pw,          # planewave cutoff
              dw=dw,          # density cutoff
              nbands=-10,     # number of bands
              kpts=kpts,      # k points
              xc=xc,          # exchange correlation method
              sigma=sigma,    # Fermi temperature
              dipole= dipole,
	            spinpol = spinpol,
              convergence={
                'energy': 0.0005,
                'mixing': mixing,
		            'nmix':10,
                'maxsteps': maxsteps,
                'diag': 'david'},
              #mode = 'scf',
              output = {'avoidio':False,
                        'removewf':True,
                        'wf_collect':False},
	             outdirprefix='vibdir'
             )

atoms.set_calculator(calc)

vib = Vibrations(atoms, delta=0.04, indices=vib_atoms)
vib.run()
vib.summary(log='vibrations.txt')
vib.write_jmol()

#dyn = QuasiNewton(atoms, logfile='qn.log', trajectory='qn.traj')
#dyn.run(fmax=0.05)

energy = atoms.get_potential_energy()

vibs = vib.get_energies()

with open('vibrations.txt','r') as f: ZPE = f.readlines()[-1]
gibbs   = HarmonicThermo(vib_energies=vibs)
entropy = gibbs.get_entropy(temperature)

with open('converged.log', 'w') as f: f.write("Energy: %f eV\nEntropy: %f eV/K \n%s"%(energy,entropy,ZPE))

