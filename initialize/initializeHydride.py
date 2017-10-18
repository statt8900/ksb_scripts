import warnings
warnings.filterwarnings("ignore", message="Moved to ase.neighborlist")
from emt  					import EMT

from contextlib import contextmanager
import sys, os

from constraint 		import *
from initializeCommon 	import *
from db 				import plotQuery,query1,insertObject
from job 				import db2object,Job
from ase.data 			import covalent_radii,chemical_symbols
from ase.db				import connect
from ase.visualize 		import view

from pymatgen.analysis.defects.point_defects 		import Interstitial
from pymatgen.transformations.site_transformations 	import InsertSitesTransformation
from pymatgen.io.ase 								import AseAtomsAdaptor

from pymatgen.symmetry.analyzer 					import SpacegroupAnalyzer
asedb = connect('/scratch/users/ksb/db/ase.db')
covalent_radii = {symb:rad for symb,rad in zip(chemical_symbols,covalent_radii)}

"""
#############
# CONSTRAINTS
#############
"""

CONSTRAINTS = [COMPLETED,PD,QE] # some constraint that detects whether past the minmum
inter = 'H'

def checkForDuplicates(ase_new,structure,emt):
	for row in  asedb.select(relaxed=False):
		if list(row.numbers) == list(ase_new.get_atomic_numbers()):
			if row.structure == structure:
				if row.emt - emt < 0.001:
					print 'Possible duplicate found (same stoich, similar EMT)'
					return False
	return True

def main():
	jIDs = [x[0] for x in plotQuery(['jobid'],CONSTRAINTS)]

	question = "Going to add an %s interstitial to %d bulk objects.\n(y/n)--> "%(inter,len(jIDs))

	if raw_input(question).lower() in ['y','yes']:
		for j in jIDs:
			aseinitID 	= query1('aseid','jobid',j)
			aseinit 	= asedb.get_atoms(id=aseinitID)

			# Access information about previous Job / Atoms object
			jobinit  	= db2object(j)
			xc,pw,kptden 	= jobinit.xc,		jobinit.pw,		jobinit.kptden
			psp,xtol,strain = jobinit.psp,		jobinit.xtol,	jobinit.strain
			precalc,dftcode = jobinit.precalc,	jobinit.dftcode
			nameinit = jobinit.name()
			structure = jobinit.structure()
			vacancies = jobinit.vacancies()

			# Create conventional PyMatGen Object
			pmg_init	= AseAtomsAdaptor.get_structure(aseinit)
			pmg_init2 	= SpacegroupAnalyzer(pmg_init).get_conventional_standard_structure()
			
			interstitial = Interstitial(pmg_init2,None,covalent_radii) #accuracy=high breaks...
			os.system('cls' if os.name == 'nt' else 'clear')
			
			for i,site in enumerate(interstitial.enumerate_defectsites()):
				coordination =  int(round(interstitial.get_defectsite_coordination_number(i)))
				mult 	= 0 # interstitial.get_defectsite_multiplicity(i) -- broken ???
				insert  = InsertSitesTransformation([inter],[site.coords],coords_are_cartesian=True)
				pmg_new = insert.apply_transformation(pmg_init2.copy())
				ase_new = AseAtomsAdaptor.get_atoms(pmg_new)
				
				ase_new.set_calculator(EMT())
				emt = ase_new.get_potential_energy()


				if   coordination == 4: siteName='T'
				elif coordination == 6: siteName='O'
				else: siteName = '%d-fold'%coordination

				question = ('site: %s\ncoordination: %s\nmultiplicity: %s'%(site.coords,coordination,mult)
							+'\ninitial ase id: %d \nxc: %s \npw: %d '%(aseinitID,xc,pw)
							+'\nkptden: %f \npsp: %s\nxtol: %f\nstrain: %f'%(kptden,psp,xtol,strain)
							+'\nprecalc: %s \ndftcode: %s'%(precalc,dftcode)
							+'\n\nDoes this structure look good?\n(y/n)--> ')

				view(ase_new)
				if raw_input(question).lower() in ['y','yes']:
					if checkForDuplicates(ase_new,structure,emt):

						newquestion 	= 'What structure does this have?\n(leave blank for general triclinic case)\n--> '
						structure 		= raw_input(newquestion)
						if structure is '': structure = 'triclinic'

						info  = {'name': 		nameinit+'_%s-%s'%(inter,siteName)
								,'emt': 		emt 				# EMT for relaxed structures useless, only relevant for deciding when to relax something
								,'relaxed': 	False 				# Could always doing this be a problem?
								,'comments':	'Generated from initializeHydride.py'
								,'parent': 		aseinitID
								,'kind': 		'bulk'
								,'structure': 	structure
								,'interstitial':siteName} 		
						if vacancies is not None: info['vacancies'] = vacancies 

						newaseid = asedb.write(ase_new,key_value_pairs=info)

						newjob = Job(None,'bulkrelax',newaseid,None,None,xc,pw,kptden,psp,xtol,strain
									,2 if xc=='mBEEF' else 1,precalc,dftcode,None,None,'initialized')
						
						insertObject(newjob)


if __name__ == '__main__':
	main()
