#!/nfs/slac/g/suncatfs/sw/py2.7.13/bin/pythonwrapper

from espresso import espresso
from ase import io

#######################################################################################
#######################################################################################
#Inputs
#######################################################################################
#Important Stuff
xc          = 'rpbe'
kpt         = (8,8,8)
pw_cutoff   = 800   # eV
spinpol     = False  #Spin-polarization
dipole      = {'status':False}
spinpol     = False

#Less-Important Stuff
nbands      = -10
sigma       = 0.1   #Fermi Electron Temperature

output      = {'avoidio':False,
                'removewf':True,
                'wf_collect':False}

calcdir  =   'cellrelax'

#Convergence
convergence={'mixing':   0.1,
             'maxsteps': 500,
             'energy':   5e-4,
             'diag':     'david',
             'mixing_mode':'plain'}


####################################################################################
####################################################################################
#Main Script
###################################################################################

atoms=io.read('init.traj')

###################################
#Create Calculator
###################################
calc = espresso(pw            = pw_cutoff,
                dw            = pw_cutoff*10,
                xc            = xc,
                kpts          = kpt,
                nbands        = nbands,
                sigma         = sigma,
                mode          = 'vc-relax',
                cell_dynamics = 'bfgs',
                opt_algorithm = 'bfgs',
                cell_factor   = 5.,
                spinpol       = spinpol,
                outdir        = calcdir,
                output        = output,
                psppath       = '/nfs/slac/g/suncatfs/sw/external/esp-psp/v2',
                convergence   = convergence
                )

atoms.set_calculator(calc)
energy = atoms.get_potential_energy() #trigger espresso to be launched
io.write('intermediate.traj',calc.get_final_structure())

atoms = io.read('intermediate.traj')
atoms.set_calculator(calc)
energy = atoms.get_potential_energy() #trigger espresso to be launched
print 'Optimized unit cell:\n', calc.get_final_structure().cell
print 'Residual stress:\n', calc.get_final_stress()
print 'Total energy (eV):\n', energy

io.write('final.traj',calc.get_final_structure())

######################
# If job is successful
######################

with open('converged.txt','w') as f:
    f.write('Optimized unit cell: %s \nResidual stress: %s\nTotal Energy (ev): %s ' %(str(calc.get_final_structure().cell),calc.get_final_stress(),energy))
