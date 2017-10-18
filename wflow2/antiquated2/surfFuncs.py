import base64,pickle
import numpy 				as np
import matplotlib.pyplot 	as plt

from ase.constraints 				import FixAtoms
from pymatgen.symmetry.analyzer 	import SpacegroupAnalyzer
from pymatgen.io.ase 				import AseAtomsAdaptor
from pymatgen.core.surface 			import SlabGenerator,Slab
from pymatgen.analysis.adsorption 	import get_rot,reorient_z,cart_to_frac,frac_to_cart,patches,color_dict,AdsorbateSiteFinder,Path

# Internal modules
import printParse,gas

"""
Functions to allow a user to create surfaces from completed bulk jobs
"""

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
	plt.gcf().set_size_inches(8,8)
	plt.savefig(pth+'/sites.png')
	with open(pth+'/sites.png', 'rb') as f: return base64.b64encode(f.read())

def plot_slab(slab, ax, scale=0.8, repeat=3, window=1, decay = 0.2):
	"""(ADAPTED FROM PYMATGEN) Function that helps visualize the slab in a 2-D plot, for convenient viewing of output of AdsorbateSiteFinder."""

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

	for site in b_sites+o_sites+h_sites: ax.text(site[0][0],site[0][1],site[1], zorder=10000,ha='center',va='center')

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
	"""H/T MSTATT"""
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
	maxTag = max([a.tag for a in atoms])
	if sym:
		if maxTag%2 != nConstrained%2: raise ValueError, 'Impossible to have symmetric slab with %n layers and %n fixed'%(maxTag,nConstrained)
		else: 	skip = (maxTag - nConstrained) / 2 #should be integer
	else: skip = 0
	constrainedLayers = range(maxTag-nConstrained+1-skip,maxTag+1-skip)
	constrainedAtoms = []
	for i,a in enumerate(atoms):
		if a.tag in constrainedLayers: constrainedAtoms.append(i)
	return constrainedAtoms

##############
# Main Scripts
##############
def makeSlab(a,facList,lay,sym,xy,vac):
	pmg_a	= AseAtomsAdaptor.get_structure(a)
	unit 	= SpacegroupAnalyzer(pmg_a,symprec=0.1, angle_tolerance=5).get_conventional_standard_structure()
	gen = SlabGenerator(unit,facList,findHeight(unit,facList,lay,sym),vac,center_slab=sym,primitive=True,lll_reduce=True)
	slabs = gen.get_slabs()
	if len(slabs)>1: print "Warning, multiple slabs generated..."
	slab = slabs[0]
	slab.make_supercell([xy[0],xy[1],1])
	return reorient_z(slab)


def bulk2surf(bulkpkl,facet,xy,layers,constrained,symmetric,vacuum,vacancies):
	magmomInit 	= 3
	magElems 	= ['Fe','Mn','Cr','Co','Ni']
	bare 		= makeSlab(pickle.loads(bulkpkl),facet,layers,symmetric,xy,vacuum)
	img  		= plot(bare)
	aseAtoms 	= AseAtomsAdaptor.get_atoms(bare)
	taggedAtoms	= tag_atoms(aseAtoms)
	magmoms 	= [magmomInit if (magmomInit and e in magElems) else 0 for e in taggedAtoms.get_chemical_symbols()] 
	taggedAtoms.set_constraint(FixAtoms(indices=constrainAtoms(taggedAtoms,constrained,symmetric)))
	taggedAtoms.set_initial_magnetic_moments(magmoms)

	taggedAtoms.wrap()

	for i in vacancies: del taggedAtoms[i]
	return taggedAtoms,img

def makePMGSlabFromASE(a,facet):
	lattice 			= a.get_cell()
	species 			= a.get_chemical_symbols()
	coords 				= a.get_positions()
	miller_index 		= facet
	oriented_unit_cell	= AseAtomsAdaptor.get_structure(a)
	shift 				= 0
	scale_factor 		= None #???????
	return Slab(lattice, species, coords, miller_index,oriented_unit_cell, shift, scale_factor,coords_are_cartesian=True)

def adsorbedSurface(surfpckl,facet,adsorbates): 
	"""
	Adds adsorbates to bare surface
	"""
	magmomInit 	= 3
	magElems 	= ['Fe','Mn','Cr','Co','Ni']
	initASE 	= pickle.loads(surfpckl)
	constrnts 	= initASE.constraints
	tags 		= initASE.get_tags().tolist()
	slab 		= makePMGSlabFromASE(initASE,facet)
	asf 		= AdsorbateSiteFinder(slab)

	b_sites = asf.find_adsorption_sites(distance=1.1,symm_reduce=0)['bridge']
	o_sites = asf.find_adsorption_sites(distance=1.1,symm_reduce=0)['ontop']
	h_sites = asf.find_adsorption_sites(distance=1.1,symm_reduce=0)['hollow']

	for ads,sites in adsorbates.items():
		a = gas.molDict[ads]
		for (kind,num) in [printParse.alphaNumSplit(x) for x in sites]:
			asf = AdsorbateSiteFinder(slab)
			if   kind == 'B': slab=asf.add_adsorbate(a,b_sites[int(num)-1])
			elif kind == 'O': slab=asf.add_adsorbate(a,o_sites[int(num)-1])
			elif kind == 'H': slab=asf.add_adsorbate(a,h_sites[int(num)-1])
			else: raise ValueError, "Bad site character in "+str(sites)

	aseSlab 	= AseAtomsAdaptor.get_atoms(slab)
	magmoms 	= [magmomInit if (magmomInit and e in magElems) else 0 for e in aseSlab.get_chemical_symbols()] 

	newtags = np.zeros(len(aseSlab))

	for i,t in enumerate(tags): newtags[i]=t
	aseSlab.set_tags(newtags)
	aseSlab.set_constraint(constrnts)
	aseSlab.set_pbc([1,1,1])
	aseSlab.set_initial_magnetic_moments(magmoms)
	aseSlab.wrap()

	return aseSlab







