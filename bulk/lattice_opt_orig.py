#!/home/vossj/suncat/bin/python
#above line selects special python interpreter needed to run espresso
#SBATCH -p iric
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
#SBATCH --time=72:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=iric
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#SBATCH --mem-per-cpu=4000
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=mstatt@stanford.edu
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks-per-node=16
#################

from ase.constraints import *
from ase.utils.geometry import *
from ase.lattice.spacegroup import *
from ase.structure import bulk
from math import *
from ase.lattice.surface import *
from ase import *
from ase.io import read
from ase.optimize import QuasiNewton
from espresso import espresso
from ase.dft.bee import BEEF_Ensemble
import cPickle as pickle


##########################
#INPUTS
######################
x0 = [3.03]
symbols = ['V']
kpts = (11,11,11)

######################

# read in trajectory file
# it can be the clean surface or
# with atoms adsorbed
try:
	atoms = read('qn.traj')
except:
	try:
		atoms = read('qn.traj.bak')
     	except:
        	try:
	  		atoms = read('init_wN.traj')
		except:
			atoms = read('init.traj')


from scipy.optimize import *
from espresso import espresso
import os
from ase.lattice import bulk



iter = 0


atoms_list = []
calc = espresso(pw=800,
        dw = 8000,
        kpts = kpts,
        nbands = -20,
        xc = 'BEEF',
        convergence = {'energy':1e-5,
                       'mixing':0.1,
                       'maxsteps':200,
                       'diag':'david'},
        dipole = {'status':False},
        spinpol = False,
        outdir = 'outdir',
        output={'removesave':True})


def get_energy(x):
    "Function to create an atoms object with lattice parameters specified by x and obtain its energy using calc."
    global iter
    iter +=1
    a = x[0]
    atoms = bulk(symbols[0], crystalstructure='fcc', a=x)
    atoms_list.append(atoms)
    atoms.set_calculator(calc)
    write('out.traj',atoms_list)
    energy = atoms.get_potential_energy()
    with paropen('lattice_opt.log','a') as logfile: logfile.write('%s\\t%s\\n' %(energy,x))
    return energy

x = fmin(get_energy,x0,ftol=1e-3) #this can be changed to minimize(...method = method...) in scipy > v0.1
print x