
from pymatgen.core.lattice        import Lattice
from pymatgen.core.structure      import Structure
from pymatgen.core.surface        import SlabGenerator
from pymatgen.io.ase              import AseAtomsAdaptor
from pymatgen.symmetry.analyzer   import SpacegroupAnalyzer
from pymatgen.analysis.adsorption import *
from pymatgen import Molecule
from pymatgen.transformations.advanced_transformations import SlabTransformation

from ase import io
import  matplotlib.pyplot as plt
import  numpy as np

#latRockSalt  = Lattice.from_parameters(3,3,3,90,90,90)
#primRockSalt = Structure(latRockSalt,['Na','Cl'],[[0,0,0],[0.5,0.5,0.5]])
#latBCC       = Lattice.from_parameters(3,3,3,90,90,90)
#primBCC      = Structure(latBCC,['Li'],[[0,0,0]])
def write(filename,x): io.write(filename,AseAtomsAdaptor.get_atoms(x))
def plot(mgStruct): 
	plot_slab2(mgStruct,plt.gca())
	plt.gcf().set_size_inches(9,9)
	plt.show()
def plot_slab2(slab, ax, scale=0.8, repeat=3, window=0.7
			  , decay = 0.2):

	"""
	Function that helps visualize the slab in a 2-D plot, for
	convenient viewing of output of AdsorbateSiteFinder.

	Args:
		slab (slab): Slab object to be visualized
		ax (axes): matplotlib axes with which to visualize
		scale (float): radius scaling for sites
		repeat (int): number of repeating unit cells to visualize
		window (float): window for setting the axes limits, is essentially
			a fraction of the unit cell limits
		decay (float): how the alpha-value decays along the z-axis
	"""
	orig_slab = slab.copy()
	slab = reorient_z(slab)
	orig_cell = slab.lattice.matrix.copy()
	if repeat: 	slab.make_supercell([repeat, repeat, 1])
	coords = np.array(sorted(slab.cart_coords, key=lambda x: x[2]))
	sites = sorted(slab.sites, key = lambda x: x.coords[2])
	alphas = 1 - decay*(np.max(coords[:, 2]) - coords[:, 2])
	alphas = alphas.clip(min=0)
	corner = [0, 0, cart_to_frac(slab.lattice, coords[-1])[-1]]
	corner = frac_to_cart(slab.lattice, corner)[:2]
	verts =  orig_cell[:2, :2]
	lattsum = verts[0]+verts[1]
	# Draw circles at sites and stack them accordingly
	for n, coord in enumerate(coords):
		r = sites[n].specie.atomic_radius*scale
		ax.add_patch(patches.Circle(coord[:2]-lattsum*(repeat//2), 
									r, color='w', zorder=2*n))
		color = color_dict[sites[n].species_string]
		ax.add_patch(patches.Circle(coord[:2]-lattsum*(repeat//2), r,
									facecolor=color, alpha=alphas[n], 
									edgecolor='k', lw=0.3, zorder=2*n+1))
	# Adsorption sites
	asf = AdsorbateSiteFinder(orig_slab)
	ads_sites = asf.find_adsorption_sites(symm_reduce=0)['all']
	sop = get_rot(orig_slab)
	ads_sites = [sop.operate(ads_site)[:2].tolist()
				 for ads_site in ads_sites]
	#ax.plot(*zip(*ads_sites), color='k', marker='$B2$',
	#		markersize=10, mew=1, linestyle='', zorder=10000)



	b_sites = asf.find_adsorption_sites(symm_reduce=0)['bridge']
	b_sites = [(sop.operate(ads_site)[:2].tolist(),'B'+str(i))
				 for i,ads_site in enumerate(b_sites)]
	o_sites = asf.find_adsorption_sites(symm_reduce=0)['ontop']
	o_sites = [(sop.operate(ads_site)[:2].tolist(),'O'+str(i))
				 for i,ads_site in enumerate(o_sites)]

	h_sites = asf.find_adsorption_sites(symm_reduce=0)['hollow']
	h_sites = [(sop.operate(ads_site)[:2].tolist(),'H'+str(i))
				 for i,ads_site in enumerate(h_sites)]

	for site in b_sites+o_sites+h_sites:
		ax.text(site[0][0]-0.3,site[0][1]-0.1,site[1], zorder=10000)


	# Draw unit cell
	verts = np.insert(verts, 1, lattsum, axis=0).tolist()
	verts += [[0., 0.]]
	verts = [[0., 0.]] + verts
	codes = [Path.MOVETO, Path.LINETO, Path.LINETO, 
			 Path.LINETO, Path.CLOSEPOLY]
	verts = [(np.array(vert) + corner).tolist() for vert in verts]
	path = Path(verts, codes)
	patch = patches.PathPatch(path, facecolor='none',lw=2,
			alpha = 0.5, zorder=2*n+2)
	ax.add_patch(patch)

	ax.set_aspect("equal")
	center = corner + lattsum / 2.
	extent = np.max(lattsum)
	lim_array = [center-extent*window, center+extent*window]
	x_lim = [ele[0] for ele in lim_array]
	y_lim = [ele[1] for ele in lim_array]
	ax.set_xlim(x_lim)
	ax.set_ylim(y_lim)
	return ax




ff = Structure(Lattice.from_parameters(2.45,2.45,2.45,60,60,60)
			,['Cu']
			,[[0,0,0]])

CO = Molecule(["C", "O"], [[0, 0, 0],[0, 0, 1.5]])
H  = Molecule(["H"],[[0,0,0]])


placement = [(CO,'ontop',0),(CO,'hollow',2),(H,'hollow',2)]

unit = SpacegroupAnalyzer(ff).get_conventional_standard_structure()

surfU = SlabGenerator(unit,[1,1,1],6,10)

slab= (surfU.get_slabs()[0])

slab=reorient_z(slab)
slab.make_supercell([2,2,1])


for p in placement:
	osites=len(AdsorbateSiteFinder(slab).find_adsorption_sites(symm_reduce=0)['ontop'])
	hsites=len(AdsorbateSiteFinder(slab).find_adsorption_sites(symm_reduce=0)['hollow'])
	bsites=len(AdsorbateSiteFinder(slab).find_adsorption_sites(symm_reduce=0)['bridge'])
	print osites,hsites,bsites
	sites=AdsorbateSiteFinder(slab).find_adsorption_sites()[p[1]]
	slab = AdsorbateSiteFinder(slab).add_adsorbate(p[0],sites[p[2]])
	

plot(slab)
write('out.traj',slab)


#symAtoms=AseAtomsAdaptor.get_atoms(ss.get_slabs(symmetrize=True)[0])
#io.write('sym.traj',symAtoms)

