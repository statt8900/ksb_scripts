#!/home/vossj/suncat/bin/python
#above line selects special python interpreter needed to run espresso
#SBATCH -p iric,owners,normal
#################
#set a job name
#you can also use --job-name=$PWD when submitting
#SBATCH --job-name=myjob
#################
#a file for job output, you can check job progress
#SBATCH --output=myjob.out
#################
# a file for errors from the job
#SBATCH --error=myjob.err
#################
#time you think you need; default is 20 hours
#SBATCH --time=01:00:00
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


from ase.data.molecules import molecule
from espresso import espresso
from ase.optimize import QuasiNewton
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.all import view,read
import os


#######################################################################################
#######################################################################################
#Inputs
#######################################################################################
#Important Stuff
xc          = 'BEEF'
kpt         = (1, 1, 1)   
pw_cutoff   = 800        # eV
dw_cutoff   = 8000       # eV
spinpol     = False      # Spin-polarization

vib_atoms   = [1]        # WHICH ATOMS ALLOWED TO VIBRATE

#Less-Important Stuff
dipole      = {'status': False}
nbands      = -60
sigma       = 0.1   #Fermi Electron Temperature
output      = {'avoidio':False,
               'removewf':True,
               'wf_collect':False}

#Convergence
convergence = {'energy':0.00005,
               'mixing':0.1,
               'nmix':10,
               'mix':4,
               'maxsteps':500,
               'diag':'david',
               'mixing_mode':'local-TF'}
fmax       = .05
calcdir    = 'calcdir'

#######################################################################################
#######################################################################################

path = os.getcwd()+'/'

files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
for f in files:
    if str(f)[-5:]=='.traj': 
        traj=str(f)[:-5]
        atoms=read(traj+'.traj')


calc0 = espresso(pw=pw_cutoff,
				 dw=dw_cutoff,
				 xc=xc,
				 kpts=kpt,
                 nbands=nbands,
                 smearing = 'gaussian',
                 sigma=sigma,
                 convergence=convergence,
				 spinpol=spinpol,
                 outdir=calcdir) 

atoms.set_calculator(calc0)

dyn = QuasiNewton(atoms,logfile='out.log',trajectory='out.traj')
dyn.run(fmax=0.01)
electronicenergy = atoms.get_potential_energy()
calc0.stop()

from espresso.vibespresso import vibespresso
#it's important to use the vibrational calculator only for one vib.run() and no relaxations or any other
#calculation, because the first calculation performed by this calculator must be for the equilibrium structure
#this is why there are two calculators here:
    #calc0 for the relaxation (in case you hadn't relaxed your structure already)
    #and calc for the vibrational analysis

calc = espresso(pw=pw_cutoff,
				 dw=dw_cutoff,
				 xc=xc,
				 kpts=kpt,
                 nbands=-60,
                 smearing = 'gaussian',
                 sigma=sigma,
                 convergence=convergence,
				 spinpol=spinpol,
                 outdir=calcdir) 
atoms.set_calculator(calc)

#indices=[1] only allows the second N atom to be displaced
#in general, you can always fix at least one atom in vibrational calculations,
#which corresponds to not calculating the three zero frequency or k=0 acoustic modes
#due to translational invariance (you'll often fix more atoms, because you're only
#interested in a small subset of the modes)
vib = Vibrations(atoms,indices=vib_atoms)
vib.run()
#this will dump the frequencies to the standard output of the job
vib.summary()
vib_energies = vib.get_energies()

thermo = IdealGasThermo(vib_energies=vib_energies,
                        electronicenergy=electronicenergy,
                        atoms=atoms,
                        geometry='linear',
                        symmetrynumber=2, spin=0)
G = thermo.get_free_energy(temperature=298.15, pressure=101325.)

e = open('e_energy.out','w')
g = open('g_energy.out','w')
e.write(str(electronicenergy))
g.write(str(G))
