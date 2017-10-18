import numpy as np    #vectors, matrices, lin. alg., etc.
import matplotlib
matplotlib.use('Agg') #turn off screen output so we can plot from the cluster
from ase.utils.eos import *  # Equation of state: fit equilibrium latt. const
from ase.units import kJ
from ase.lattice import bulk
from ase import *
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac, MethfesselPaxton


#######################################################################################
#######################################################################################
#Inputs
#######################################################################################
a0 = 3.52

metal   = 'Ni'
crystal = 'fcc'

al = bulk(metal, crystal, a=a0)
cell0 = al.cell

#Important Stuff
xc          = 'PBE'
kpt         = (8,8,8)
pw_cutoff   = 800   # eV
dw_cutoff   = 8000  # eV
N_steps = 11 # ODD number of steps tested for latt. consts.; assume spacing .01 A

spinpol     = False  #Spin-polarization
#dipole      = {'status':False}
#Less-Important Stuff
#nbands      = -10
#sigma       = 0.1   #Fermi Electron Temperature
#magmom      = 2.0   #Initial guess

#######################################################################################
#######################################################################################


calc = GPAW(mode=PW(pw_cutoff),
               xc=xc,
               kpts=kpt,
               parallel={'band': 1},    #WHAT DOES THIS MEAN 
               spinpol = spinpol,
               txt='Al-%d.txt' % pw_cutoff)




strains = np.linspace(-0.01*(N_steps-1)/2,0.01*(N_steps-1)/2,N_steps)
volumes = []  #we'll store unit cell volumes and total energies in these lists
energies = []

for i in strains: #loop over scaling factors
    #build Pt unit cell
    atoms = bulk(metal, crystal, a0+i)
    atoms.set_pbc((1,1,1))                #periodic boundary conditions about x,y & z
    atoms.set_calculator(calc)            #connect espresso to Pt unit cell
    volumes.append(atoms.get_volume())    #append the current unit cell volume
                                          #to list of volumes
    energy=atoms.get_potential_energy()   #append total energy to list of
    energies.append(energy)               #energies

eos = EquationOfState(volumes, energies) #Fit calculated energies at different
v0, e0, B = eos.fit()                    #lattice constants to an
                                         #equation of state


# setup bulk using optimized lattice and save it

best_a = (4.*v0)**(1./3.) # Angstroms
atoms = bulk(metal, crystal, best_a)
atoms.write('bulk.traj')

#output of lattice constant = cubic root of volume of conventional unit cell
#fcc primitive cell volume = 1/4 * conventional cell volume
print 'Lattice constant:', best_a, 'AA'
print 'Bulk modulus:', B / kJ * 1e24, 'GPa'
print '(Fitted) total energy at equilibrium latt. const.:', e0, 'eV'
eos.plot(metal+'-bulk-eos.png')    #create a png plot of eos fit
