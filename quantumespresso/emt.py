"""Effective medium theory potential."""

from math import sqrt, exp, log

import numpy as np

from ase.data import chemical_symbols
from ase.units import Bohr
from ase.calculators.neighborlist import NeighborList
from ase.calculators.calculator import Calculator


parameters = {
    #      E0     s0    V0     eta2    kappa   lambda  n0
    #      eV     bohr  eV     bohr^-1 bohr^-1 bohr^-1 bohr^-3
    'H':  (-3.21, 1.31, 0.132, 2.652,  2.790,  3.892,  0.00547),
    'Li': (-3, 1.1, 0.1,3, 3.1, 4.3, 0.1 ),
    'Be': (-4, 1.2, 0.1,3, 3.2, 4.2, 0.01 ),
    'B':  (-5, 1.3, 0.1,3, 3.3, 4.1, 0.1 ),
    'C':  (-3.50, 1.81, 0.332, 1.652,2.790,1.892,0.01322),
    'N':  (-5.10, 1.88, 0.132, 1.652,2.790,1.892,0.01222),
    'O':  (-4.60, 1.95, 0.332, 1.652,2.790,1.892,0.00850),
    'F':  (-6, 1.4, 0.1,3, 2.34, 4.4, 0.1 ),
    'Na': (-7, 1.5, 0.1,3, 2.34, 4.3, 0.01 ),
    'Mg': (-8, 1.6, 0.1,3, 2.33, 4.2, 0.1 ),
    'Al': (-3.28, 3.00, 1.493, 1.240,2.000,1.169,0.00700),
    'Si': (-9, 1.7, 0.1,3, 3, 4, 0.1 ),
    'P':  (-3, 1.1, 1.1,4, 4, 5, 0.01 ),
    'S':  (-4, 1.2, 0.4,4, 4, 5, 0.1 ),
    'Cl': (-5, 1.3, 0.1,2,4.9, 5, 0.01 ),
    'K':  (-6, 1.4, 2.6,4.3,4.8, 5, 0.1 ),
    'Ca': (-7, 1.5, 0.1,4.2,4.7, 5, 0.01 ),
    'Sc': (-8, 1.6, 3.1,1,4.6, 5.6, 0.1 ),
    'Ti': (-9, 1.7, 0.1,4,4.5, 5.2, 0.01 ),
    'V':  (-3, 1.8, 2.1,2,2.4, 3.3, 0.1 ),
    'Cr': (-4, 1.9, 2.1,2,2.3, 3.4, 0.01 ),
    'Mn': (-5, 2.0, 0.1,1,2.24, 3, 0.1 ),
    'Fe': (-6, 2.1, 3.1,2,2.2, 3, 0.01 ),
    'Co': (-7, 2.2, 1.1,2,2.1, 3, 0.1 ),
    'Ni': (-4.44, 2.60, 3.673, 1.669,2.757,1.948,0.01030),
    'Cu': (-3.51, 2.67, 2.476, 1.652,2.740,1.906,0.00910),
    'Zn': (-8, 2.3, 0.1,2, 2, 1.32, 0.01 ),
    'Ga': (-9, 2.4, 1.1,2, 2.2, 2.34, 0.1 ),
    'Ge': (-3, 2.53,0.2, 2.2,2, 0.3, 0.01 ),
    'As': (-4, 2.6, 0.3, 2.4,2.3, 3.3, 0.1 ),
    'Se': (-5, 2.7, 0.6,2, 2, 1.3, 0.01 ),
    'Br': (-6, 2.75,0.1,2,2.4, 4.3, 0.1 ),
    'Rb': (-7, 2.8, 0.7,1.5, 2, 3.3, 0.01 ),
    'Sr': (-8, 2.9, 0.1,1.53, 2.2, 3.1, 0.1 ),
    'Y':  (-9, 1.0, 0.7,1.3,2.3, 3.2, 0.01 ),
    'Zr': (-8, 1.2, 0.6,1.4,2.4, 3.3, 0.01 ),
    'Nb': (-7, 1.4, 0.5,1.6, 2.5, 3.4, 0.1 ),
    'Mo': (-6, 2.65,0.4,1.7,2.9, 3.5, 0.01 ),
    'Tc': (-5, 2.8, 0.3,1.3, 4.2, 3.7, 0.01 ),
    'Ru': (-3.3, 3.1, 1.21, 1.17,1.87,3.2,0.07), ### SAME AS Au
    'Rh': (-3, 3.2, 0.2,2.2, 1.2, 3.1, 0.01 ),
    'Pd': (-3.90, 2.87, 2.773, 1.818,3.107,2.155,0.00688),
    'Ag': (-2.96, 3.01, 2.132, 1.652,2.790,1.892,0.00547),
    'Cd': (-2, 3.4, 0.4,3.3, 1.5, 2.3, 0.03 ),
    'In': (-1, 2.5, 1.5,1.1, 1.7, 1.3, 0.01 ),
    'Sn': (-9, 2.2, 2.2,1.2, 1.19, 2.31, 0.001 ),
    'Sb': (-1, 2.3, 3.2,1.3, 1.29, 1.32, 0.003 ),
    'Te': (-2, 2.4, 1.2,1.4, 1.39, 2.33, 0.005 ),
    'I':  (-3, 2.5, 2.2,1.1, 1.49, 1.34, 0.007 ),
    'Cs': (-4, 2.2, 3.2,1.2, 1.59, 2.35, 0.002 ),
    'Ba': (-5.1, 2.3, 3.2,1.3, 2.19, 1.63, 0.034 ),
    'La': (-6, 2.35,1.1,0.3, 1.99, 4.63, 0.024 ),
    'Hf': (-6, 2.4, 2.2,1.4, 1.29, 2.36, 0.006 ),
    'Ta': (-7, 2.5, 3.2,1.1, 1.39, 1.27, 0.008 ),
    'W':  (-1, 2.2, 1.2,1.2, 1.49, 1.38, 0.003 ),
    'Re': (-2, 2.3, 2.2,1.34,1.59, 2.39, 0.004 ),
    'Os': (-3, 2.4, 3.2,1.4, 1.19, 1.31, 0.005 ),
    'Ir': (-4, 2.5, 1.2,1.3, 1.29, 2.32, 0.006 ),
    'Pt': (-5.85, 2.90, 4.067, 1.812,3.145,2.192,0.00802),
    'Au': (-3.80, 3.00, 2.321, 1.674,2.873,2.182,0.00703),
    'Pb': (-5.80, 2.30, 2.721, 1.174,5.823,4.182,0.01703)}

beta = 1.809     # (16 * pi / 3)**(1.0 / 3) / 2**0.5,
                 # but preserve historical rounding


class EMT(Calculator):
    implemented_properties = ['energy', 'forces']

    nolabel = True

    def __init__(self):
        Calculator.__init__(self)

    def initialize(self, atoms):
        self.par = {}
        self.rc = 0.0
        self.numbers = atoms.get_atomic_numbers()
        maxseq = max(par[1] for par in parameters.values()) * Bohr
        rc = self.rc = beta * maxseq * 0.5 * (sqrt(3) + sqrt(4))
        rr = rc * 2 * sqrt(4) / (sqrt(3) + sqrt(4))
        self.acut = np.log(9999.0) / (rr - rc)
        for Z in self.numbers:
            if Z not in self.par:
                p = parameters[chemical_symbols[Z]]
                s0 = p[1] * Bohr
                eta2 = p[3] / Bohr
                kappa = p[4] / Bohr
                x = eta2 * beta * s0
                gamma1 = 0.0
                gamma2 = 0.0
                for i, n in enumerate([12, 6, 24]):
                    r = s0 * beta * sqrt(i + 1)
                    x = n / (12 * (1.0 + exp(self.acut * (r - rc))))
                    gamma1 += x * exp(-eta2 * (r - beta * s0))
                    gamma2 += x * exp(-kappa / beta * (r - beta * s0))

                self.par[Z] = {'E0': p[0],
                               's0': s0,
                               'V0': p[2],
                               'eta2': eta2,
                               'kappa': kappa,
                               'lambda': p[5] / Bohr,
                               'n0': p[6] / Bohr**3,
                               'rc': rc,
                               'gamma1': gamma1,
                               'gamma2': gamma2}
                #if rc + 0.5 > self.rc:
                #    self.rc = rc + 0.5

        self.ksi = {}
        for s1, p1 in self.par.items():
            self.ksi[s1] = {}
            for s2, p2 in self.par.items():
                #self.ksi[s1][s2] = (p2['n0'] / p1['n0'] *
                #                    exp(eta1 * (p1['s0'] - p2['s0'])))
                self.ksi[s1][s2] = p2['n0'] / p1['n0']

        self.forces = np.empty((len(atoms), 3))
        self.sigma1 = np.empty(len(atoms))
        self.deds = np.empty(len(atoms))

        self.nl = NeighborList([0.5 * self.rc + 0.25] * len(atoms),
                               self_interaction=False)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=['positions', 'numbers', 'cell',
                                  'pbc', 'charges','magmoms']):
        Calculator.calculate(self, atoms, properties, system_changes)

        if 'numbers' in system_changes:
            self.initialize(self.atoms)

        positions = self.atoms.positions
        numbers = self.atoms.numbers
        cell = self.atoms.cell

        self.nl.update(self.atoms)

        self.energy = 0.0
        self.sigma1[:] = 0.0
        self.forces[:] = 0.0

        natoms = len(self.atoms)

        for a1 in range(natoms):
            Z1 = numbers[a1]
            p1 = self.par[Z1]
            ksi = self.ksi[Z1]
            neighbors, offsets = self.nl.get_neighbors(a1)
            offsets = np.dot(offsets, cell)
            for a2, offset in zip(neighbors, offsets):
                d = positions[a2] + offset - positions[a1]
                r = sqrt(np.dot(d, d))
                if r < self.rc + 0.5:
                    Z2 = numbers[a2]
                    p2 = self.par[Z2]
                    self.interact1(a1, a2, d, r, p1, p2, ksi[Z2])

        for a in range(natoms):
            Z = numbers[a]
            p = self.par[Z]
            try:
                ds = -log(self.sigma1[a] / 12) / (beta * p['eta2'])
            except (OverflowError, ValueError):
                self.deds[a] = 0.0
                self.energy -= p['E0']
                continue
            x = p['lambda'] * ds
            y = exp(-x)
            z = 6 * p['V0'] * exp(-p['kappa'] * ds)
            self.deds[a] = ((x * y * p['E0'] * p['lambda'] + p['kappa'] * z) /
                            (self.sigma1[a] * beta * p['eta2']))
            self.energy += p['E0'] * ((1 + x) * y - 1) + z

        for a1 in range(natoms):
            Z1 = numbers[a1]
            p1 = self.par[Z1]
            ksi = self.ksi[Z1]
            neighbors, offsets = self.nl.get_neighbors(a1)
            offsets = np.dot(offsets, cell)
            for a2, offset in zip(neighbors, offsets):
                d = positions[a2] + offset - positions[a1]
                r = sqrt(np.dot(d, d))
                if r < self.rc + 0.5:
                    Z2 = numbers[a2]
                    p2 = self.par[Z2]
                    self.interact2(a1, a2, d, r, p1, p2, ksi[Z2])

        self.results['energy'] = self.energy
        self.results['forces'] = self.forces

    def interact1(self, a1, a2, d, r, p1, p2, ksi):
        x = exp(self.acut * (r - self.rc))
        theta = 1.0 / (1.0 + x)
        y1 = (0.5 * p1['V0'] * exp(-p2['kappa'] * (r / beta - p2['s0'])) *
              ksi / p1['gamma2'] * theta)
        y2 = (0.5 * p2['V0'] * exp(-p1['kappa'] * (r / beta - p1['s0'])) /
              ksi / p2['gamma2'] * theta)
        self.energy -= y1 + y2
        f = ((y1 * p2['kappa'] + y2 * p1['kappa']) / beta +
             (y1 + y2) * self.acut * theta * x) * d / r
        self.forces[a1] += f
        self.forces[a2] -= f
        self.sigma1[a1] += (exp(-p2['eta2'] * (r - beta * p2['s0'])) *
                            ksi * theta / p1['gamma1'])
        self.sigma1[a2] += (exp(-p1['eta2'] * (r - beta * p1['s0'])) /
                            ksi * theta / p2['gamma1'])

    def interact2(self, a1, a2, d, r, p1, p2, ksi):
        x = exp(self.acut * (r - self.rc))
        theta = 1.0 / (1.0 + x)
        y1 = (exp(-p2['eta2'] * (r - beta * p2['s0'])) *
              ksi / p1['gamma1'] * theta * self.deds[a1])
        y2 = (exp(-p1['eta2'] * (r - beta * p1['s0'])) /
              ksi / p2['gamma1'] * theta * self.deds[a2])
        f = ((y1 * p2['eta2'] + y2 * p1['eta2']) +
             (y1 + y2) * self.acut * theta * x) * d / r
        self.forces[a1] -= f
        self.forces[a2] += f