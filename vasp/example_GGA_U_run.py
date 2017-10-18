from ase.io import *
from ase.io.trajectory import PickleTrajectory
from ase.visualize import view
from ase.io import write
from ase import Atoms
from ase.io import read
from ase import Atoms 
from ase import Atom 
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.optimize import BFGS
from ase.visualize import *
from ase.io import *
import numpy as np

import os
import subprocess

name = 'vasp_run2'
#atoms=read('../../../qn.traj')
atoms=read('restart.traj')

for a in atoms:
	b=a.magmom
	a.magmom=round(b, 2)
	print round(b, 2)
#	if a.symbol=='Ni':
#		if abs(a.magmom) <1.5:
#			a.magmom=1.0
#		else:
#			if a.magmom >0.0:
#				a.magmom=1.8
#			else:
#				a.magmom=-1.8
#	else:
#		a.magmom=0.0

#atoms.center(vacuum=5.5, axis=2)

#del atoms[38]
#del atoms[37]
#del atoms[36]

atoms.write(name+'_init'+'.traj')
atoms.write(name+'_init'+'.cif')

#print cell

#vasp_calculator.int_keys.append("nedos")
#vasp_calculator.int_keys.append("ICHARG")
#vasp_calculator.float_keys.append("AMIX")
#vasp_calculator.float_keys.append("AMIX_MAG")
#vasp_calculator.float_keys.append("BMIX")
#vasp_calculator.float_keys.append("BMIX_MAG")
#vasp_calculator.float_keys.append("POTIM")
#vasp_calculator.int_keys.append("INIMIX")
#vasp_calculator.bool_keys.append("LASPH")
#vasp_calculator.bool_keys.append("LWAVE")
#vasp_calculator.string_keys.append("LREAL")
#vasp_calculator.int_keys.append("kpar")
#vasp_calculator.int_keys.append("ncore")
#vasp_calculator.int_keys.append("NUPDOWN")

from ase.calculators.vasp import Vasp
#import ase.calculators.vasp as vasp_calculator
calc = Vasp(encut=500,
			xc='PBE',
            gga='PE',
            setups={'K':'_pv', 'Nb':'_pv', 'Zr':'_sv'}, 
			kpts  = (4,4,1),
            kpar = 10,
			npar = 8,
			gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
			ismear=0,
			inimix = 0,
			amix   = 0.2,
			bmix   = 0.00001,
			amix_mag= 0.1,
			bmix_mag= 0.00001,
            #NUPDOWN= 0,
            nelm= 250 ,	
            sigma = 0.05,
            algo = 'Normal',
			ibrion=2,
            #POTIM=1.0,
			ediffg=-0.02,  # forces
			ediff=1e-4,  #energy conv.
            nedos=2001,
			prec='Normal',
			#nsw=0, # don't use the VASP internal relaxation, only use ASE
			nsw=300, # don't use the VASP internal relaxation, only use ASE
			lvtot=False,
			ispin=2,
            #ldau=False, 
            ldau=True, 
            ldautype=2,
            lreal  = 'auto',
			lasph  = True, 
			lwave  = False, 
            #ldau_luj={'Ni':{'L':2,  'U':9.0, 'J':0.0},
            ldau_luj={'Ni':{'L':2,  'U':6.45, 'J':0.0},
				'Co':{'L':2,  'U':3.32, 'J':0.0},
                'Cr':{'L':2,  'U':3.5,  'J':0.0},   
				'Fe':{'L':2,  'U':5.3, 'J':0.0},
				'V':{'L':2,  'U':3.10, 'J':0.0},
				'Mn':{'L':2,  'U':3.75, 'J':0.0},
				'Ti':{'L':2,  'U':3.00, 'J':0.0},
				'W':{'L':-1,  'U':0.0, 'J':0.0},
                 'O':{'L':-1, 'U':0.0, 'J':0.0},
                 'C':{'L':-1, 'U':0.0, 'J':0.0},
				'Au':{'L':-1, 'U':0.0, 'J':0.0},
				'Cu':{'L':-1, 'U':0.0, 'J':0.0},
                'Ce':{'L':3,  'U':5.5, 'J':1.0},
				'H':{'L':-1, 'U':0.0, 'J':0.0}},
             ldauprint=2,
             #idipol=3,
              #dipol=(0, 0, 0.5),
              #ldipol=True,
              lmaxmix=4,
              lorbit=11)

atoms.set_calculator(calc)
atoms.get_potential_energy()

#dyn = QuasiNewton(atoms, logfile=name+'.log', trajectory=name+'.traj')
#dyn.run(fmax=0.05)

write('final.traj',atoms)

import os
import subprocess
#subprocess.call('rm -f WAVECAR', shell=True)
#write(name+'opt.traj',atoms)
