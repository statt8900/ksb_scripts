import numpy as np
from math 			import sqrt
from db  			import getCluster
from ase.parallel 	import parprint
from ase  			import optimize


class Job(object):

	def gpawRestart(self):
		from gpaw import restart,PW,Davidson,Mixer,MixerSum,FermiDirac

		spinpol = any([x>0 for x in self.magmoms])

		return restart('preCalc_inp.gpw',mode    = PW(self.pw)
										,xc      = self.xc
										,kpts    = self.kpt
										,spinpol = spinpol
										,convergence = {'energy':self.econv} #eV/electron
										,mixer =  ((MixerSum(beta=self.mixing,nmaxold=self.nmix,weight=100)) 
												if spinpol else (Mixer(beta=self.mixing,nmaxold=self.nmix,weight=100)))
										,maxiter       = self.maxstep
										,nbands        =  -1*self.nbands
										,occupations   = FermiDirac(self.sigma)
										,setups        = self.psp #(pspDict[calc.psp]).pthFunc[cluster]
										,eigensolver   = Davidson(5)
										,poissonsolver = None # {'dipolelayer': 'xy'} if isinstance(self,SurfaceJob) else None
										,txt='%d_%s.txt'%(self.pw,self.xc)
										,symmetry={'do_not_symmetrize_the_density': True}) #ERROR IN LI bcc 111 surface
										

	def gpawCalc(self,xc,spinpol):
		from gpaw import GPAW,PW,Davidson,Mixer,MixerSum,FermiDirac
		return GPAW(mode         = PW(self.pw)                         #eV
					,xc          = xc
					,kpts        = self.kpt
					,spinpol     = spinpol
					,convergence = {'energy':self.econv} #eV/electron
					,mixer       = ((MixerSum(beta=self.mixing,nmaxold=self.nmix,weight=100)) 
									if spinpol else (Mixer(beta=self.mixing,nmaxold=self.nmix,weight=100)))
					,maxiter       = self.maxstep
					,nbands        =  -1*self.nbands
					,occupations   = FermiDirac(self.sigma)
					,setups        = 'sg15' 
					,eigensolver   = Davidson(5)
					,poissonsolver = None # {'dipolelayer': 'xy'} if isinstance(self,SurfaceJob) else None
					,txt='%d_%s.txt'%(self.pw,xc)
					,symmetry={'do_not_symmetrize_the_density': True}) #ERROR IN LI bcc 111 surface
					
	def qeCalc(self,xc,spinpol,restart):
		from espresso import espresso	
				
		pspDict = 	{'sherlock':
						{'gbrv15pbe':'/home/vossj/suncat/psp/gbrv1.5pbe'}
					,'suncat':
						{'gbrv15pbe':'/nfs/slac/g/suncatfs/sw/external/esp-psp/gbrv1.5pbe'}}

		cluster = getCluster()
		psppath = pspDict[cluster][self.psp]

		return espresso( pw      = self.pw
						,dw      = self.pw*10
						,xc      = xc
						,kpts    = self.kpt
						,spinpol = spinpol
						,convergence = {'energy':self.econv
										,'mixing':self.mixing
										,'nmix':self.nmix
										,'maxsteps':self.maxstep
										,'diag':'david'}
						,nbands = -1*self.nbands
						,sigma  = self.sigma
						,dipole = {'status': False}
						,outdir = 'calcdir'	 # output directory
						,startingwfc='file' if restart else 'atomic+random'
						,psppath= psppath
						,output = {'removesave':True})

	def calc(self,precalc=False,restart=False):
		""" Creates an ASE calculator object """

		xc = self.precalcxc if precalc else self.xc
		spinpol = any([x>0 for x in self.magmoms])

		if   self.dftcode == 'gpaw':            return self.gpawCalc(xc,spinpol)
		elif self.dftcode == 'quantumespresso': return self.qeCalc(xc,spinpol,restart)
		else: raise ValueError, 'Unknown dftcode '+self.dftcode

	def optimizePos(self,atoms,calc,saveWF=False):
		atoms.set_calculator(calc)

		maxForce = np.amax(abs(atoms.get_forces()))
		if maxForce > self.fmax:
			parprint("max force = %f, optimizing positions"%(maxForce))
			dyn = optimize.BFGS(atoms=atoms, logfile='qn.log', trajectory='qn.traj',restart='qn.pckl')
			dyn.run(fmax=self.fmax)

		if saveWF and self.dftcode=='gpaw' and self.kind=='bulk': 
			atoms.calc.write('preCalc_inp.gpw', mode='all') #for use in getXCContribs
			atoms,calc = self.gpawRestart()





def kptdensity2monkhorstpack(atoms, kptdensity, even=True):
    """Convert k-point density to Monkhorst-Pack grid size.
    atoms: Atoms object  ---  Contains unit cell and information about boundary conditions.
    kptdensity: float    ---  Required k-point density.  Default value is 3.5 point per Ang^-1.
    even: bool           ---  Round up to even numbers.
    """
    recipcell = atoms.get_reciprocal_cell()
    kpts = []
    for i in range(3):
        if atoms.pbc[i]:
            k = 2 * 3.14159 * sqrt((recipcell[i]**2).sum()) * kptdensity
            if even:
                kpts.append(2 * int(np.ceil(k / 2)))
            else:
                kpts.append(int(np.ceil(k)))
        else:
            kpts.append(1)
    return kpts


