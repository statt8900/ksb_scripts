#!/nfs/slac/g/suncatfs/sw/py2.7.13/bin/pythonwrapper

import ase,espresso
print "ASE FILE IS ",ase.__file__
print "ESPRESSO FILE IS ",espresso.__file__

from espresso import espresso
from ase 	  import optimize,io



pw      = 500
xc      = 'BEEF'
kpt     = (4,4,1)
spinpol = False
dipole  = True

eConv    = 5e-5
mixing   = 0.1
nmix     = 10
maxsteps = 500
nbands   = -12
sigma    = 0.2

atoms = io.read('init.traj')

calc = espresso(pw      =  pw
				,dw      = pw*10
				,xc      = xc
				,kpts    = kpt
				,spinpol = spinpol
				,convergence = {'energy':eConv
								,'mixing':mixing
								,'nmix': nmix
								,'maxsteps':maxsteps
								,'diag':'david'}
				,nbands = nbands
				,sigma  = sigma
				,dipole = {'status':dipole}
				,outdir = 'calcdir'	 #output directory
            	#,psppath= 
				,output = {'removesave':True})
atoms.set_calculator(calc)
dyn = optimize.BFGS(atoms=atoms, logfile='qn.log', trajectory='qn.traj',restart='qn.pckl')
dyn.run(fmax=0.05)
io.write('out.traj',calc.get_final_structure())
with open('converged.log','w') as f: f.write(str(atoms.get_potential_energy()))
