#!/home/vossj/suncat/bin/python
#SBATCH -p iric,owners
#SBATCH --exclusive
#SBATCH --job-name=name
#SBATCH --output=myjob.out
#SBATCH --error=myjob.err
#SBATCH --time=00:00:10
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=FAIL
##SBATCH  --mail-user=ksb@stanford.edu
#SBATCH --ntasks-per-node=16

from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator

from ase.io.trajectory import PickleTrajectory
from ase.io import read,write

atoms = read('POSCAR')

calc = vasp_calculator.Vasp(encut=500,
                  xc = 'PBE',  gga ='BF',luse_vdw=True,zab_vdw=-1.8867,
                  kpts  = (4,3,1),
                  kpar  = 7,           # use this if you run on one node (most calculations).  see suncat confluence page for optimal setting
                  npar  = 1,           # use this if you run on one node (most calculations).  see suncat confluence page for optimal setting
                  gamma = True,     # Gamma-centered (defaults to Monkhorst-Pack)
                  ismear= 0,
                  algo  = 'all',
                  nelm  = 250,
                  sigma = 0.05,
                  ibrion= 2,
                  nelmdl= -9,
                  #isif  = 3, # 3 for vc-relax, otherwise delete the flag
                  ediffg= -0.05,     # forces
                  ediff = 1e-4,       # energy conv. both of these are for the internal relaxation, ie nsw
                  prec  = 'Accurate',
                  nsw   = 0,            # don't use the VASP internal relaxation, only use ASE
                  ispin = 1,
                  lreal = 'auto',  # automatically decide to do real vs recip space calc
                  ldipol= True,
                  lvhar = True,
                  dipol = (0.0,0.0,0.0),
                  idipol= 3,
                  icharg= 1)         # start from CHGCAR if icharg = 1

atoms.set_calculator(calc)
atoms.get_potential_energy()

write('opt.traj',atoms)