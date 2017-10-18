from db import sqlexecute
from ase.db 	import connect

"""
Common elements from InitializeJobs and InitializeSurfaceJobs
"""

"""
###########
# Domains #
###########
"""
asedb = connect('/scratch/users/ksb/db/ase.db')
######
# Bulk
######

trajDomain 		= [x.id for x in asedb.select()] 
dftDomain    	= ['gpaw','quantumespresso']
xcDomain     	= ['PBE','BEEF','mBEEF','RPBE']
preXcDomain 	= ['PBE','None']
pspDomain   	= ['gbrv15pbe','sg15']
pwDomain    	= [400,500,700,900,1100,1300,1500,1700]
kptDomain   	= [2,4] 	# pts / A^-1

##########
# Surfaces
##########

facetDomain 		= [111,100]
xyDomain 			= [11,22]
layerDomain 		= [3,4]
constrainedDomain 	= [2]
symmetricDomain 	= [True,False]
vacuumDomain 		= [10]
vacancyDomain 		= [[],[1]] #indices of Atoms objects. Will I ever need to consider multiple vacancies?

adsDomain 			= [{},{'H':['B1']}]

