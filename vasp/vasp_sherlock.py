#!/home/vossj/suncat/bin/python
#SBATCH -p iric,owners
#SBATCH --exclusive
#SBATCH --job-name=name
#SBATCH --output=myjob.out
#SBATCH --error=myjob.err
#SBATCH --time=00:00:10
#SBATCH --qos=normal
#SBATCH --nodes=7
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=FAIL
#SBATCH  --mail-user=ksb@stanford.edu
#SBATCH --ntasks-per-node=16

from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator

from ase.io.trajectory import PickleTrajectory
from ase             import Atoms
from ase             import Atom
from ase.io          import read,write
from ase.constraints import FixAtoms
from ase.optimize    import QuasiNewton
from ase.optimize    import BFGS

name = 'Cu-fcc'

try:
    # restarting from an old calculation
    atoms=read('POSCAR')
except:
    # new calculation
    atoms=read('init.traj')

calc = vasp_calculator.Vasp(encut=500,
                        xc='BEEF',luse_vdw=True,zab_vdw=True,zab_vdw=-1.8867,gga='BF',
                        kpts  = (4,3,1),
                        npar=2, 
                        gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                        ismear=0,
                        algo = 'fast',
                        nelm=250,
                        sigma = 0.05,
                        ibrion=2,
                        ediffg=-0.05,  # forces
                        ediff=1e-4,  #energy conv. both of these are for the internal relaxation, ie nsw
                        prec='Accurate',
                        nsw=0, # don't use the VASP internal relaxation, only use ASE
                        ispin=1,
                        ldipol=True,
                        lvhar=True,
                       	dipol=(0.5,0.5,0.5),
                        idipol=3,
                        icharg=1) #start from CHGCAR if icharg = 1

atoms.set_calculator(calc)
atoms.get_potential_energy()

write(name+'_opt.traj',atoms)