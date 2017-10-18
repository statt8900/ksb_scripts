from ase.io import read,write
import ase.calculators.vasp as vasp_calculator

atoms = read('init.traj')

#vasp_calculator.int_keys.append('kpar')
calc = vasp_calculator.Vasp(encut=500
                  ,xc = 'PBE',  gga ='BF',luse_vdw=True,zab_vdw=-1.8867
                  ,kpts  = (4,3,1)
                  ,kpar  = 1               # use this if you run on one node (most calculations).  see suncat confluence page for optimal setting
                  ,npar  = 4              # use this if you run on one node (most calculations).  see suncat confluence page for optimal setting
                  ,lhfcalc=True            # HSE parameter
                  ,hfscreen=0.3            # HSE parameter
                  ,setups={'K':'_pv', 'Nb':'_pv', 'Zr':'_sv'} # don't use default for these, _sv = include semi-core electrons too, _h and _s (hard vs soft)?
                  ,gamma = True           # Gamma-centered (defaults to Monkhorst-Pack)
                  ,ismear=0
                  #,inimix = 0            # density mixing
                  ,amix   = 0.2           # etc.
                  ,bmix   = 0.00001       #
                  ,amix_mag= 0.1          # magnetic mixing
                  ,bmix_mag= 0.00001
                  #NUPDOWN= 0,
                  ,nelm= 250              # max electronic steps
                  ,sigma = 0.05           # fermi temp
                  ,algo = 'Normal'        # diagonalization
                  ,ibrion=2               # conjugate gradient descent
                  #POTIM=1.0,
                  ,ediffg=-0.02           # force cutoff (IF NEGATIVE)
                  ,ediff=1e-4             # energy conv.
                  ,nedos=2001             # # gridpoints in dos
                  ,prec='Normal'          # precision
                  #,nsw=0,                # don't use the VASP internal relaxation, only use ASE
                  ,nsw=300                # don't use the VASP internal relaxation, only use ASE
                  ,lvtot=False            # never want this to be true (if you want work function, use lvhar=True to write locpot)
                  ,ispin=2                # turns on spin polarization (1 = false 2 = true)

                  ,ldau=True              # turns on + U
                  ,ldautype=2             # parameter ?
                  ,lreal  = 'auto'        # automatically decide to do real vs recip space calc
                  ,lasph  = True          # "also includes non-spherical contributions from the gradient corrections inside the PAW spheres. "
                  ,lwave  = False         # don't write out wf

                  ,ldau_luj={'Ni':{'L':2,  'U':6.45, 'J':0.0}, # how to apply U correction (what shell, how much, etc.)
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
                        'H':{'L':-1, 'U':0.0, 'J':0.0}}
                  ,ldauprint=2                 # prints u values to outcar
                  ,idipol=3                   # what direction (z)
                  ,dipol=(0, 0, 0.5)         # location
                  ,ldipol=True               # turns it on
                  ,lmaxmix=4
                  ,lorbit=11) # prints magnetic moments in outcar

atoms.set_calculator(calc)
atoms.get_potential_energy()

write('opt.traj',atoms)

