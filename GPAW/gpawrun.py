from ase import io
atoms = io.read('POSCAR')
from gpaw import GPAW, PW, FermiDirac, restart, Davidson, Mixer, MixerSum, MixerDif
from ase.optimize import QuasiNewton
from ase.parallel import paropen

kpts=(13,13,13)

atoms.calc = GPAW(setups='sg15',mode=PW(800),xc='PBE',eigensolver=Davidson(5), mixer=Mixer(0.1, 5, 100), kpts=kpts,
    occupations=FermiDirac(0.2),nbands=-16,txt='pbe.txt')

atoms.get_potential_energy()
atoms.calc.write('inp.gpw', mode='all')

#This uses libxc for MBEEF, you can use for relaxations, etc... (might be
#a little bit faster)
#atoms,calc = restart('inp.gpw', setups='sg15', xc='MGGA_X_MBEEF+GGA_C_PBE_SOL', convergence={'energy': 5e-6}, txt='mbeef.txt')
#But for xc energy contribution output (i.e. single-point calculations)
#only the python version works:
atoms,calc = restart('inp.gpw', setups='sg15', xc='mBEEF', convergence={'energy': 5e-6}, txt='mbeef.txt')

f = paropen('total_energy_and_xc_coefficients.txt', 'w')
print >>f, atoms.get_potential_energy()

from gpaw.xc.bee import BEEFEnsemble

beef = BEEFEnsemble(calc)

print >>f, beef.mbeef_exchange_energy_contribs()