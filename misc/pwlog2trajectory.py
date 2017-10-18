

#****************************************************************************
# Copyright (C) 2013 SUNCAT
# This file is distributed under the terms of the
# GNU General Public License. See the file `COPYING'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#****************************************************************************


from sys import argv, exit, stderr
if len(argv) != 3:
    print >>stderr, 'usage: ' + argv[0] + ' pw-log output.traj'
    exit(1)

import numpy as np
from ase.io import PickleTrajectory
from ase import Atoms
from os import popen



import ase.units
hartree = ase.units.Hartree
rydberg = ase.units.Rydberg
bohr = ase.units.Bohr


def get_total_energy(s):
    a = s.readline()
    while a != '' and a[:17] != '!    total energy':
        a = s.readline()
    if a == '':
        return None
    energy = float(a.split()[-2]) * rydberg
    for i in range(11):
        a = s.readline()
        if '     smearing contrib.' in a:
            break
    # correct for finite temperature entropy term
    # in case of finite temp. smearing
    if a[:22] == '     smearing contrib.':
        energy -= 0.5 * float(a.split()[-2]) * rydberg
    return energy


# dummy calculator
class calculator:

    def set_energy(self, e):
        self.energy = e
    def get_potential_energy(self):
        return self.energy
    def get_calculator(self):
        return self
    def notimpl(self, apply_constraint=False):
        raise NotImplementedError

p = popen('grep -n Giannozzi ' + argv[1] + ' 2>/dev/null | tail -1', 'r')
try:
    n = int(p.readline().split()[0].strip(':'))
except:
    print >>stderr, 'No valid pw-log at ' + argv[1] + ' found.'
    p.close()
    exit(2)
p.close()

calc = calculator()

s = open(argv[1], 'r')
# skip over previous runs in log in case the current log has been
# appended to old ones
for i in range(n):
    s.readline()

a = s.readline()
while a[:11] != '     celldm':
    a = s.readline()
alat = float(a.split()[1]) / 1.889726
a = s.readline()
while a[:12] != '     crystal':
    a = s.readline()
cell = []
for i in range(3):
    cell.append([float(x) for x in s.readline().split()[3:6]])
cell = np.array(cell)
a = s.readline()
while a[:12] != '     site n.':
    a = s.readline()
pos = []
syms = ''
y = s.readline().split()
while len(y) > 0:
    nf = len(y)
    pos.append([float(x) for x in y[nf - 4:nf - 1]])
    syms += y[1].strip('0123456789')
    y = s.readline().split()
pos = np.array(pos) * alat
natoms = len(pos)

# create atoms object with coordinates and unit cell
# as specified in the initial ionic step in log
atoms = Atoms(syms, pos, cell=cell * alat, pbc=(1, 1, 1))

atoms.get_calculator = calc.get_calculator
atoms.get_potential_energy = calc.get_potential_energy
atoms.get_forces = calc.notimpl
atoms.get_stress = calc.notimpl
atoms.get_charges = calc.notimpl

# get total energy at first ionic step
en = get_total_energy(s)
if en is not None:
    calc.set_energy(en)
else:
    print >>stderr, 'no total energy found'
    exit(3)

print(atoms)

traj = PickleTrajectory(argv[2], 'w')


traj.write(atoms)

a = s.readline()
while a != '':
    while a[:7] != 'CELL_PA' and a[:7] != 'ATOMIC_' and a != '':
        a = s.readline()
    if a == '':
        break
    if a[0] == 'A':
        coord = a.split('(')[-1]
        for i in range(natoms):
            pos[i][:] = s.readline().split()[1:4]

        if coord == 'alat)':
            atoms.set_positions(pos * alat)
        elif coord == 'bohr)':
            atoms.set_positions(pos * bohr)
        elif coord == 'angstrom)':
            atoms.set_positions(pos)
        else:
            atoms.set_scaled_positions(pos)
        # get total energy at 2nd,3rd,...nth ionic step
        en = get_total_energy(s)
        calc.set_energy(en)

        if en is not None:
            traj.write(atoms)
        else:
            break
    else:
        for i in range(3):
            cell[i][:] = s.readline().split()
        atoms.set_cell(cell * alat, scale_atoms=False)
    a = s.readline()
