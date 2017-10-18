import os
import sqlite3 as sqlite
from db import insertObject
from itertools import product

from bulk import BulkTraj
from calc import Calc
from psp  import PSP
from xc   import XC

commands = []

sherlockDB = '/scratch/users/ksb/db/data.db'

###########################
# Add PSP Paths to Database
###########################
# Define domain of psp's to consider
gbrv15pbe = PSP('gbrv15pbe','/home/vossj/suncat/psp/gbrv1.5pbe',[1,0,3,4,3,4,5,6,7,0,9,10,3,4,5,6,7,0,9,10,11,12,13,14,15,16,17,18,19,20,19,12,5,6,7,0,9,10,11,12,13,14,15,16,15,16,19,12,13,14,15,6,7,0,9,10,11,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,13,11,15,16,15,16,11,12,13,14,8,0,0,0,0,0,0,0,0,0,0,8,8])
sgs15     = PSP('sgs_15','/scratch/users/ksb/gpaw/gpaw_sg15/norm_conserving_setups',[1,0,3,4,3,4,5,6,7,0,9,10,3,4,5,6,7,0,9,10,11,12,13,14,15,16,17,18,19,20,19,12,5,6,7,0,9,10,11,12,13,14,15,16,15,16,19,12,13,14,15,6,7,0,9,10,11,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,13,11,15,16,15,16,11,12,13,14,8,0,0,0,0,0,0,0,0,0,0,8,8])

pspDomain = [gbrv15pbe,sgs15]

for p in pspDomain: insertObject(p)

################################
# Add XC Functionals to Database
################################
# Define domain of XC's to consider
PBE   = XC('PBE',  [[]])
BEEF  = XC('BEEF', [[]])
mBEEF = XC('mBEEF',[[]])

xcDomain = [PBE,BEEF,mBEEF]

for x in xcDomain:  insertObject(x)

################################
# Add Calculators to Database
################################
# Define domain of calcs to consider
xcNameDomain  = [x.name for x in xcDomain]
pspNameDomain = [p.name for p in pspDomain]
pwDomain      = [500,700,900,1100,1300,1500,1700,1900] # eV
kptDomain     = [6] # converted to varchar (e.g. '5-5-5'), 6 -> 15x15x15 for 1 atom unit cell
eConvDomain   = [0.0005]  # eV
fMaxDomain    = [0.05]     # eV/A
xTolDomain    = [0.005]    # Angstrom
strainDomain  = [0.02]     # Fraction, +/- 1 over which bulk modulus fit is calculated
mixDomain     = [0.1]      # 0.05, #0.1 ? 
nMixDomain    = [5]        # Int
maxStepDomain = [500]      # Int
nBandDomain   = [12]       # Int
sigmaDomain   = [0.1]      # 0.1
magDomain     = [0,3]      # Bohr magneton

combs = product(xcNameDomain,pspNameDomain,pwDomain,kptDomain,eConvDomain
				,fMaxDomain,xTolDomain,mixDomain,nMixDomain
				,maxStepDomain,nBandDomain,sigmaDomain,magDomain)

for c in combs: 
	insertObject(Calc(*((None,)+c)))


#########################
trajRoot = '/scratch/users/ksb/db/traj'
trajDomain = [BulkTraj(None,trajRoot+'/'+x) for x in os.listdir(trajRoot)]

for t in trajDomain: 	insertObject(t)


