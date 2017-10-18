#!/nfs/slac/g/suncatfs/sw/py2.7.13/bin/pythonwrapper

from ase import io
from espresso.vibespresso import vibespresso
from ase.vibrations       import Vibrations
from ase.thermochemistry  import HarmonicThermo


##################
# Inputs #########
##################
vib_atoms = [0] # List of atoms allowed to vibrate

#Convergence
##################
# Slab ###########
##################

atoms = io.read('qn.traj') # optimized

atoms.set_masses()

from espresso import espresso

calc = vibespresso(
              pw=800,        # planewave cutoff
              dw=8000,        # density cutoff
              nbands=-10,   # number of bands
              kpts=(8,8,8),    # k points
              xc='rpbe',        # exchange correlation method
              sigma=0.2,    # Fermi temperature
              dipole={'status': False},
	            spinpol = False,
              convergence={
                'energy': 0.0005,
                'mixing': 0.1,
		            'nmix':10,
                'maxsteps': 500,
                'diag': 'david'},
	             outdirprefix='vibdir'
             )

atoms.set_calculator(calc)

vib = Vibrations(atoms, delta=0.04, indices=vib_atoms)
vib.run()
vib.summary(log='vibrations.txt')
vib.write_jmol()

realFreq = [x*0.00012 for x in vib.get_frequencies() if not isinstance(x,complex)] #eV

S = HarmonicThermo(realFreq).get_entropy(300)
with open('vibrations.txt','a') as f: f.write('\nTS at 300K: '+str(S))

