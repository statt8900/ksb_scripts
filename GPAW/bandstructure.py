

from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac

"""In this tutorial we calculate the electronic band structure of Si along high symmetry directions in the Brillouin zone.

First, a standard ground state calculations is performed and the results are saved to a .gpw file. As we are dealing with small bulk system, plane wave mode is the most appropriate here."""

# Perform standard ground state calculation (with plane wave basis)
si = bulk('Si', 'diamond', 5.43)
calc = GPAW(mode=PW(200),
            xc='PBE',
            kpts=(8, 8, 8),
            random=True,  # random guess (needed if many empty bands required)
            occupations=FermiDirac(0.01),
            txt='Si_gs.txt')
si.calc = calc
si.get_potential_energy()
calc.write('Si_gs.gpw')


"""
Next, we calculate eigenvalues along a high symmetry path in the Brillouin zone kpts={'path': 'GXWKL', 'npoints': 60}. See ase.dft.kpoints.special_points for the definition of the special points for an FCC lattice.

For the band structure calculation, density is fixed to the previously calculated ground state density (fixdensity=True), and as we want to calculate all k-points, symmetry is not used (symmetry='off'). The unoccupied states can be sometimes converged faster with the conjugate gradient eigensolver.
"""

# Restart from ground state and fix potential:
calc = GPAW('Si_gs.gpw',
            nbands=16,
            fixdensity=True,
            symmetry='off',
            kpts={'path': 'GXWKL', 'npoints': 60},
            convergence={'bands': 8})

calc.get_potential_energy()

# Plot
bs = calc.band_structure()
bs.plot(filename='bandstructure.png', show=True, emax=10.0)
