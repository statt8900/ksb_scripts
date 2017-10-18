from ase.build import mx2
from gpaw import GPAW, FermiDirac
from ase.dft.kpoints import bandpath, special_points
from gpaw.unfold import Unfold, find_K_from_k,plot_spectral_function
from ase.units import Hartree
import pickle


structure = mx2(formula='MoS2', kind='2H', a=3.184, thickness=3.127,
                size=(3, 3, 1), vacuum=7.5)
structure.pbc = (1, 1, 1)

# Create vacancy
del structure[2]

calc = GPAW(mode='pw', #changed from LCAO with dbz basis set
            xc='PBE',
            kpts=(4, 4, 1),
            occupations=FermiDirac(0.01),
            txt='gs_3x3_defect.txt')

structure.set_calculator(calc)
structure.get_potential_energy()
calc.write('gs_3x3_defect.gpw', 'all')


a = 3.184
PC = mx2(a=a).get_cell(complete=True)
path = [special_points['hexagonal'][k] for k in 'MKG']
kpts, x, X = bandpath(path, PC, 48)
    
M = [[3, 0, 0], [0, 3, 0], [0, 0, 1]]

Kpts = []
for k in kpts:
    K = find_K_from_k(k, M)[0]
    Kpts.append(K)

calc_bands = GPAW('gs_3x3_defect.gpw',
                  fixdensity=True,
                  kpts=Kpts,
                  symmetry='off',
                  nbands=220,
                  convergence={'bands': 200})

calc_bands.get_potential_energy()
calc_bands.write('bands_3x3_defect.gpw', 'all')

unfold = Unfold(name='3x3_defect',
                calc='bands_3x3_defect.gpw',
                M=M,
                spinorbit=False)


unfold.spectral_function(kpts=kpts, x=x, X=X,
                         points_name=['M', 'K', 'G'])



calc = GPAW('gs_3x3_defect.gpw', txt=None)
ef = calc.get_fermi_level()

plot_spectral_function(filename='sf_3x3_defect',
                       color='blue',
                       eref=ef,
                       emin=-3,
                       emax=3)









