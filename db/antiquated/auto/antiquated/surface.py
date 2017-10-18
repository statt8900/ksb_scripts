import copy
import numpy as np

from objectClass      import Object
from printParse       import underscoreSep,derscoreSep,printList,parseList,printAds,parseAds,parseSym,alpha,digit
from molecules        import molDict
from miniDatabase     import trajDict
from pymatgen.analysis.adsorption import get_rot,reorient_z,cart_to_frac,frac_to_cart,patches,color_dict,AdsorbateSiteFinder,Path

class SurfacePure(Object):
	def __init__(self,bulkJob,facet,scale,symmetric,constrained,vacuum,adsorbates):
		self.trajName    = None
		self.bulkJob     = bulkJob     # Job
		self.facet       = facet       # (Int,Int,Int)
		self.scale       = scale       # (Int,Int,Int)
		self.symmetric   = symmetric   # Bool
		self.constrained = constrained # Int
		self.vacuum      = vacuum      # Int
		self.adsorbates  = adsorbates  # {String : [String]} 

	def __str__(self):  return ('PureSurface_'
								+str(self.bulkJob)+'---'
								+underscoreSep(2,	[printList(self.facet)
													,printList(self.scale)
													,'S' if self.symmetric else 'A'
													,str(self.constrained)
													,str(self.vacuum)
													,printAds(self.adsorbates)]))

	def bulkResult(self): return self.bulkJob.getBulkResult()
	def bulkUnit(self):   
		return self.bulkResult().makeUnitCell() 
	def atomList(self):   return self.makeAtoms().get_chemical_symbols()

	def makeSurfs(self):
		from pymatgen.core.surface import SlabGenerator
		return SlabGenerator(	self.bulkUnit()
								,self.facet
								,self.findHeight()
								,self.vacuum,center_slab=self.symmetric,primitive=True,lll_reduce=True)#,,,,) # maxnormalsearch makes fcc go in multiples of 3...why?

	def findHeight(self):
		"""
		Finds the slab thickness input necessary for SlabGenerator to produce slab with correct number of layers
		"""
		from pymatgen.core.surface import SlabGenerator
		x0,dx,n = 1,1,None # Angstrom
		getN = lambda x: len(SlabGenerator(self.bulkUnit(),self.facet,x,4,primitive=True).get_slabs(symmetrize=self.symmetric)[0])
		while n != self.scale[2]:
			#print x0,n
			x0+=dx
			n = getN(x0)
			if n > 20: raise ValueError, 'over ten atoms in primitive surface column ... primitive cell may contain multiple atoms ('+str(getN(0.1))+') not compatible with requested # of layers ('+str(self.scale[2])+')'
		return x0

	def bareSurface(self):
		"""
		Creates bare pymatgen slab object and scales to correct size
		"""

		slab = self.makeSurfs().get_slabs()[0]
		slab.make_supercell([self.scale[0],self.scale[1],1])
		return reorient_z(slab)

	def adsorbedSurface(self): 
		"""
		Adds adsorbates to bare surface
		"""
		from pymatgen.analysis.adsorption import AdsorbateSiteFinder,get_rot
		slab = copy.deepcopy(self.bareSurface())
		asf = AdsorbateSiteFinder(slab)

		b_sites = asf.find_adsorption_sites(distance=1.1,symm_reduce=0)['bridge']
		o_sites = asf.find_adsorption_sites(distance=1.1,symm_reduce=0)['ontop']
		h_sites = asf.find_adsorption_sites(distance=1.1,symm_reduce=0)['hollow']

		for ads,sites in self.adsorbates.items():
			a = molDict[ads]
			for (kind,num) in map(alphaNumSplit,sites):
				asf = AdsorbateSiteFinder(slab)
				if   kind == 'B': slab=asf.add_adsorbate(a,b_sites[int(num)])
				elif kind == 'O': slab=asf.add_adsorbate(a,o_sites[int(num)])
				elif kind == 'H': slab=asf.add_adsorbate(a,h_sites[int(num)])
				else: raise ValueError, "Bad site character in "+str(sites)
		return slab
	
	def makeAtoms(self): 
		"""
		Creates ASE atoms object based on input specifications
		"""
		from pymatgen.io.ase import AseAtomsAdaptor
		from ase.constraints import FixAtoms

		atoms = AseAtomsAdaptor.get_atoms(self.adsorbedSurface())
		atoms.set_tags(self.assignTags(atoms))
		atoms.set_constraint(FixAtoms(indices=self.constrainAtoms(atoms)))
		return atoms

	def guessTraj(self,pth): 
		from ase import io
		io.write(pth+'/init.traj',self.makeAtoms())

	def constrainAtoms(self,atoms):
		constraints = []
		def modTag(tag):
			if self.symmetric:
				mid = self.scale[2]/2.0+0.5
				diff = tag - mid
				return abs(floor(diff) if diff > 0 else ceil(diff) )+1
			else: return tag
		
		for i, atom in enumerate(atoms):	
			if modTag(atom.tag) <= self.constrained: constraints.append(i)

		return constraints

	def assignTags(self,atoms):
		decimalPointTol = 1; shift = 0
		zs = map(lambda x: round(x[2]+shift,decimalPointTol),atoms.get_positions())
		zs = sorted(list(set(zs)))
		tagDict = dict(zip(zs,range(len(zs))))
		return   map(lambda x: tagDict[round(x[2],1)]+1,atoms.get_positions())
		

	def makeAdsorptionPlot(self,pth): plot(self.bareSurface(),pth)

	def elementList(self):   return self.makeAtoms().get_chemical_symbols()
	def getNelectrons(self): return sum(map(lambda x: nElecDict[x],self.elementList()))

class SurfaceTraj(Object):
	def __init__(self,trajName,scale,adsorbates):
		self.trajName    = trajName    # String
		self.scale       = scale       # (Int,Int,Int)
		self.adsorbates  = adsorbates  # {String : [String]} 
		self.bulkJob     = None        # Job
		self.facet       = None        # (Int,Int,Int)
		self.symmetric   = None        # Bool
		self.constrained = None        # Int
		self.vacuum      = None        # Int

	def __str__(self):  return  'TrajSurface_'+underscoreSep(2,[self.trajName,printList(self.scale),printAds(self.adsorbates)])

	def atomList(self): return self.makeAtoms().get_chemical_symbols()

	def traj(self): return trajDict[self.trajName]
	def trajInput(self): 
		from ase import io
		return io.read(self.traj().path)

	def constrainAtoms(self): return self.trajInput().constraints

	def trajToMG(self):
		from pymatgen.io.ase import AseAtomsAdaptor
		from pymatgen.core.surface import Slab
		atoms = self.trajInput()
		struct= AseAtomsAdaptor.get_structure(atoms)
		return Slab(struct.lattice,struct.species,struct.cart_coords,(1,1,0),struct,3,(1,1,1),coords_are_cartesian=True,to_unit_cell=True)
		
	def bareSurface(self):
		"""
		Creates bare pymatgen slab object and scales to correct size
		"""
		slab = self.trajToMG()
		slab.make_supercell([self.scale[0],self.scale[1],1])
		return reorient_z(slab)

	def adsorbedSurface(self): 
		"""
		Adds adsorbates to bare surface
		"""
		from pymatgen.analysis.adsorption import AdsorbateSiteFinder,get_rot
		slab = copy.deepcopy(self.bareSurface())
		asf = AdsorbateSiteFinder(slab)

		b_sites = asf.find_adsorption_sites(distance=1.1,symm_reduce=0)['bridge']
		o_sites = asf.find_adsorption_sites(distance=1.1,symm_reduce=0)['ontop']
		h_sites = asf.find_adsorption_sites(distance=1.1,symm_reduce=0)['hollow']

		for ads,sites in self.adsorbates.items():
			a = molDict[ads]
			for (kind,num) in map(alphaNumSplit,sites):
				asf = AdsorbateSiteFinder(slab)
				if   kind == 'B': slab=asf.add_adsorbate(a,b_sites[int(num)])
				elif kind == 'O': slab=asf.add_adsorbate(a,o_sites[int(num)])
				elif kind == 'H': slab=asf.add_adsorbate(a,h_sites[int(num)])
				else: raise ValueError, "Bad site character in "+str(sites)
		return slab
	
	def makeAtoms(self): 
		"""
		Creates ASE atoms object based on input specifications
		"""
		from pymatgen.io.ase import AseAtomsAdaptor
		atoms = AseAtomsAdaptor.get_atoms(self.adsorbedSurface())
		atoms.set_constraint(self.constrainAtoms())
		return atoms

	def guessTraj(self,pth): 
		from ase import io
		io.write(pth+'/init.traj',self.makeAtoms())
		

	def makeAdsorptionPlot(self,pth): plot(self.bareSurface(),pth)

####################################################
####################################################
####################################################
####################################################
def parseSurface(x):
	if x[0]=='P': return parsePureSurface(x)
	elif x[0]=='T': return parseTrajSurface(x)
	else: raise ValueError, 'Bad str for surface '+x

def parsePureSurface(x):
	from job              import parseBulkJob
	x = x[12:] #chop off 'PureSurface_'
	bulkJob,surf  = x.split('---')[0],x.split('---')[1]
	facet,scale,sym,constrained,vacuum,adsorbates = derscoreSep(2,surf)
	return SurfacePure(parseBulkJob(bulkJob)
				,parseList(facet)
				,parseList(scale)
				,parseSym(sym)
				,int(constrained)
				,int(vacuum)
				,parseAds(adsorbates))

def parseTrajSurface(x):
	x = x[12:] #chop off 'TrajSurface_'
	trajName,scale,ads = derscoreSep(2,x)
	return SurfaceTraj(trajName,parseList(scale),parseAds(ads))


def alphaNumSplit(x): return (alpha(x),digit(x))

def plot(mgStruct,pth): 
	import  matplotlib.pyplot as plt
	plot_slab(mgStruct,plt.gca())
	plt.gcf().set_size_inches(9,9)
	plt.savefig(pth+'/sites.png')

def plot_slab(slab, ax, scale=0.8, repeat=3, window=1
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
