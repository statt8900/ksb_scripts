import warnings
warnings.filterwarnings("ignore", message="Moved to ase.neighborlist")

import sys,base64
import numpy 				as np
import matplotlib.pyplot 	as plt

from itertools 			import product 	# For cartesian product
from operator 			import mul 		# For product of list

from ase 				import io
from ase.db 			import connect
from ase.visualize 		import view
from ase.constraints 	import FixAtoms
from emt  				import EMT
from constraint 		import *
from filters 			import *
from db 				import plotQuery

from pymatgen 						import Structure,Molecule
from pymatgen.symmetry.analyzer 	import SpacegroupAnalyzer
from pymatgen.io.ase 				import AseAtomsAdaptor
from pymatgen.core.surface 			import SlabGenerator
from pymatgen.analysis.adsorption 	import get_rot,reorient_z,cart_to_frac,frac_to_cart,patches,color_dict,AdsorbateSiteFinder,Path


"""
Functions to allow a user to create surfaces from completed bulk jobs

Constraints are used to select subsets of completed bulk jobs.

Filters are used to select subsets of the parameter space 
	for creating a surface (e.g. #layers, vacuum, facet)

Like initializing BulkJobs, there are domains for the parameter space 
	that need to be defined.

PyMatGen functions are used when possible. 

Output: addition to ase.db


Instead of using strings to store all special info, use data?
>>> db.write(atoms, ..., data={'parents': [7, 34, 14], 'stuff': ...})


"""

USER_APPROVAL 	= True
magmomInit 		= 3
magElems 		= ['Fe','Mn','Cr','Co','Ni']

#############################
# Constraints on Job that
# Generated Initial Structure
##############################
essentialConstraints 	= [completed] #don't change
constraintsDuJour 	 	= [al,rpbe,qe]

initialConstraints = essentialConstraints + constraintsDuJour

#################
# Surface domains
#################
bulkDomain 			= [x[0] for x in plotQuery('bulkresult',['aseid'],initialConstraints)]
facetDomain 		= [[1,1,1],[1,0,0]]
xyDomain 			= [[1,1],[2,2]]
layerDomain 		= [3,4]
constrainedDomain 	= [2]
symmetricDomain 	= [True,False]
vacuumDomain 		= [10]
vacancyDomain 		= [[],[1]] #indices of Atoms objects. Will I ever need to consider multiple vacancies?

###########
# Filters #
###########

filtsDuJour = 	[facetFilter([1,1,1])
				,xyFilter([2,2])
				,layerFilter(4)
				,slabSymFilter(False)
				,vacuumFilter(10)
				,vacancyFilter([])]

filts = andFilters(essentialBareSurfFilters + filtsDuJour)

filtBulk 		= [x for x in bulkDomain 		if filts['aseid'](x)] 
filtFacet	 	= [x for x in facetDomain 		if filts['facet'](x)]
filtXY 			= [x for x in xyDomain 			if filts['xy'](x)]
filtLayer		= [x for x in layerDomain 		if filts['layers'](x)]
filtConstrained = [x for x in constrainedDomain if filts['constrained'](x)]
filtSym 		= [x for x in symmetricDomain 	if filts['symmetric'](x)]
filtVacuum 		= [x for x in vacuumDomain 		if filts['vacuum'](x)]
filtVacancy 	= [x for x in vacancyDomain 	if filts['vacancy'](x)]

domains = [filtBulk,filtFacet,filtXY,filtLayer,filtConstrained,filtSym,filtVacuum,filtVacancy]

domainProduct = product(*domains)

ncombo 		= reduce(mul,[len(x) for x in domains])
nBulk 		= len(filtBulk)
##############
# Adsorbates #
##############
class Adsorbate(object):
	def __init__(self,pmg,vector): 
		self.pmg  	= pmg #should be centered around 0,0,0
		self.vector = vector
		self.offset = offset
		
CO = Molecule(["C", "O"], [[0, 0, 0],[0, 0, 1.3]])
H  = Molecule(["H"],[[0,0,0]])

molDict = {'CO':CO,'H':H}

#####################################
# Functions for constructing surfaces
#####################################
def findHeight(unit,facet,layers,sym):
	"""
	Finds the slab thickness input necessary for SlabGenerator to produce slab with correct number of layers
	"""
	x0,dx,n = 1,1,None # Angstrom
	getN 	= lambda x: len(SlabGenerator(unit,facet,x,4,primitive=True).get_slabs(symmetrize=sym)[0])
	while n != layers:
		x0+=dx
		n = getN(x0)
		if n > 20: raise ValueError, ('over ten atoms in primitive surface column'
										+' ... primitive cell may contain multiple atoms ('
										+str(getN(0.1))+') not compatible with requested # of layers ('+str(layers)+')')
	return x0

def plot(mgStruct,pth='/scratch/users/ksb/img/'): 
	plot_slab(mgStruct,plt.gca())
	plt.gcf().set_size_inches(9,9)
	plt.savefig(pth+'/sites.png')
	with open(pth+'/sites.png', 'rb') as f: return base64.b64encode(f.read())

def plot_slab(slab, ax, scale=0.8, repeat=3, window=1, decay = 0.2):
	"""
	Function that helps visualize the slab in a 2-D plot, for
	convenient viewing of output of AdsorbateSiteFinder.
	Args:
		slab (slab): 	Slab object to be visualized
		ax (axes): 		matplotlib axes with which to visualize
		scale (float): 	radius scaling for sites
		repeat (int): 	number of repeating unit cells to visualize
		window (float): window for setting the axes limits, is essentially a fraction of the unit cell limits
		decay (float): 	how the alpha-value decays along the z-axis
	"""

	orig_slab 	= slab.copy()
	slab 		= reorient_z(slab)
	orig_cell 	= slab.lattice.matrix.copy()
	if repeat: 	slab.make_supercell([repeat, repeat, 1])
	coords 	= np.array(sorted(slab.cart_coords, key=lambda x: x[2]))
	sites 	= sorted(slab.sites, key = lambda x: x.coords[2])
	alphas 	= 1 - decay*(np.max(coords[:, 2]) - coords[:, 2])
	alphas 	= alphas.clip(min=0)
	corner 	= [0, 0, cart_to_frac(slab.lattice, coords[-1])[-1]]
	corner 	= frac_to_cart(slab.lattice, corner)[:2]
	verts 	= orig_cell[:2, :2]
	lattsum = verts[0]+verts[1]
	# Draw circles at sites and stack them accordingly
	for n, coord in enumerate(coords):
		r = sites[n].specie.atomic_radius*scale
		ax.add_patch(patches.Circle(coord[:2]-lattsum*(repeat//2), r, color='w', zorder=2*n))
		color = color_dict[sites[n].species_string]
		ax.add_patch(patches.Circle(coord[:2]-lattsum*(repeat//2), r,facecolor=color, alpha=alphas[n], edgecolor='k', lw=0.3, zorder=2*n+1))
	# Adsorption sites
	asf 		= AdsorbateSiteFinder(orig_slab)
	ads_sites 	= asf.find_adsorption_sites(symm_reduce=0)['all']
	sop 		= get_rot(orig_slab)
	ads_sites 	= [sop.operate(ads_site)[:2].tolist() for ads_site in ads_sites]

	b_sites = asf.find_adsorption_sites(symm_reduce=0)['bridge']
	b_sites = [(sop.operate(ads_site)[:2].tolist(),'B'+str(i)) for i,ads_site in enumerate(b_sites)]
	o_sites = asf.find_adsorption_sites(symm_reduce=0)['ontop']
	o_sites = [(sop.operate(ads_site)[:2].tolist(),'O'+str(i)) for i,ads_site in enumerate(o_sites)]
	h_sites = asf.find_adsorption_sites(symm_reduce=0)['hollow']
	h_sites = [(sop.operate(ads_site)[:2].tolist(),'H'+str(i)) for i,ads_site in enumerate(h_sites)]

	for site in b_sites+o_sites+h_sites:
		ax.text(site[0][0],site[0][1],site[1], zorder=10000,ha='center',va='center')

	# Draw unit cell
	verts = np.insert(verts, 1, lattsum, axis=0).tolist()
	verts += [[0., 0.]]
	verts = [[0., 0.]] + verts
	codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
	verts = [(np.array(vert) + corner).tolist() for vert in verts]
	path  = Path(verts, codes)
	patch = patches.PathPatch(path, facecolor='none',lw=2,alpha = 0.5, zorder=2*n+2)
	ax.add_patch(patch)

	ax.set_aspect("equal")
	center 		= corner + lattsum / 2.
	extent 		= np.max(lattsum)
	lim_array 	= [center-extent*window, center+extent*window]
	x_lim 		= [ele[0] for ele in lim_array]
	y_lim 		= [ele[1] for ele in lim_array]
	ax.set_xlim(x_lim);	ax.set_ylim(y_lim)
	return ax

def reindex_by_height(atoms):
	h 			= [atom.position[2] for atom in atoms]
	atoms_new 	= atoms.copy()
	ind 		= np.argsort(h)
	del atoms_new[:]
	for i in ind: atoms_new.append(atoms[i])
	return atoms_new

def tag_atoms(atoms, num_bins=20):
	# reindex the input Atoms object by height and tag using evenly placed bins
	# that divide the layers into 20
	atoms 	= reindex_by_height(atoms)
	height 	= [round(atom.position[2],4) for atom in atoms]
	binH 	= np.digitize(height,np.linspace(min(height)*0.99,max(height)*1.01,num_bins))
	binind 	= list(np.sort(list(set(binH))))
	binind.reverse()
	for i in range(len(binH)):
		atoms[i].tag = binind.index(binH[i])+1
	return atoms

def constrainAtoms(atoms,nConstrained,sym):
	constraints = []
	def modTag(tag):
		if sym: raise NotImplementedError
		else: return tag

	maxTag = max([a.tag for a in atoms])
	constrainedLayers = range(maxTag-nConstrained+1,maxTag+1)
	constrainedAtoms = []
	for i,a in enumerate(atoms):
		if a.tag in constrainedLayers: constrainedAtoms.append(i)
	return constrainedAtoms

#############
# Main Script
#############

def main():
	asedb = connect('/scratch/users/ksb/db/ase.db')
	initialQuestion = ('\n%d relaxed bulk structure(s) passed filters+constraints.'%(nBulk)
						+'\nDo you want to create %d slabs?\n(y/n)--> '%(ncombo))

	if raw_input(initialQuestion).lower() in ['y','yes']:

		for combo in domainProduct:
			ind,fac,xy,l,c,sym,vac,vacancy = combo

			a 		= asedb.get_atoms(id=ind)
			pmg_a	= AseAtomsAdaptor.get_structure(a)
			unit 	= SpacegroupAnalyzer(pmg_a,symprec=0.1, angle_tolerance=5).get_conventional_standard_structure()
			
			gen = SlabGenerator(unit,fac,findHeight(unit,fac,l,sym)
								,vac,center_slab=sym,primitive=True,lll_reduce=True)

			slabs = gen.get_slabs()

			if len(slabs)>1: print "Warning, multiple slabs generated..."
			
			slab = slabs[0]
			slab.make_supercell(list(xy)+[1])
			bare = reorient_z(slab)
			img  = plot(bare)

			aseAtoms 	= AseAtomsAdaptor.get_atoms(bare)
			taggedAtoms = tag_atoms(aseAtoms)
			taggedAtoms.set_constraint(FixAtoms(indices=constrainAtoms(taggedAtoms,c,sym)))

			magmoms = [magmomInit if (magmomInit and e in magElems) else 0 for e in taggedAtoms.get_chemical_symbols()] 

			taggedAtoms.set_initial_magnetic_moments(magmoms)
			taggedAtoms.set_calculator(EMT())
			taggedAtoms.wrap()
			
			for i in vacancy: del taggedAtoms[i]

			info = 	{'name': 		'%s_%s_%s'%(asedb.get(ind).get('name'),'-'.join([str(x) for x in fac]),'x'.join([str(x) for x in xy]))
					,'kind': 		'surface' 
					,'structure':	asedb.get(ind).get('structure')
					,'sites': 		img
					,'bare': 		True
					,'parent': 		ind 		#ase object surface was generated from
					,'facet': 		str(fac)
					,'xy': 			str(xy)
					,'layers': 		l
					,'constrained': c
					,'symmetric':	sym
					,'vacuum':		vac
					,'vacancies':	str(vacancy)
					,'surfcomments': 	'Autogenerated by initializeBareSurf'
					,'emt': 		taggedAtoms.get_potential_energy()
					,'relaxed': 	False}


			if USER_APPROVAL:
				view(taggedAtoms)
				question = ('Does this structure look right?\nfacet= %s\nlayers = %d\nfixed = %d\n'%('-'.join([str(x) for x in fac]),l,c)
							+'vacuum = %d\nxy = %s\n(y/n) --> '%(vac,'x'.join([str(x) for x in xy])))

				if raw_input(question).lower() in ['y','yes']:
					asedb.write(taggedAtoms,key_value_pairs=info)
			else: asedb.write(taggedAtoms,key_value_pairs=info)

if __name__=='__main__':
	main()

