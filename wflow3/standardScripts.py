
#External Modules
import os,fireworks,json,glob
#Adding comments in kris
#I want this only to show up in kris_edit
#################################
# Auxillary functions for scripts
#################################
def mergeDicts(listDicts): 
    import itertools
    return dict(itertools.chain.from_iterable([x.items() for x in listDicts])) #presumes no overlap in keys

def getCluster():
    import os
    hostname = os.environ['HOSTNAME'].lower()
    if      'sh'    in hostname: return 'sherlock'
    elif   'gpu-15' in hostname: return 'sherlock'
    elif    'su'    in hostname: return 'suncat' #important to distinguish suncat2 and 3?
    elif    'kris'  in hostname: return 'kris'
    else: raise ValueError, "getCluster did not detect SH or SU in %s"%hostname

def rank():
    import sys
    # Check for special MPI-enabled Python interpreters:
    if '_gpaw' in sys.builtin_module_names:
        import _gpaw        # http://wiki.fysik.dtu.dk/gpaw
        world = _gpaw.Communicator()
    elif '_asap' in sys.builtin_module_names:
        import _asap # http://wiki.fysik.dtu.dk/asap, can't import asap3.mpi here (import deadlock)
        world = _asap.Communicator()
    elif 'asapparallel3' in sys.modules: # Older version of Asap
        import asapparallel3
        world = asapparallel3.Communicator()
    elif 'Scientific_mpi' in sys.modules:
        from Scientific.MPI import world
    elif 'mpi4py' in sys.modules:
        world = MPI4PY()
    else:
        from ase.parallel import DummyMPI
        world = DummyMPI()# This is a standard Python interpreter:
    rank = world.rank
    size = world.size
    return rank

#############################################
# Universally-applicable functions for scripts
#############################################

def writeScript(listOfFuncs,mainFunc):
    import inspect

    with open('script.py','w') as f:
        for fnc in listOfFuncs:
             f.write('\n\n'+inspect.getsource(fnc))
        f.write("\n\nif __name__ == '__main__': %s()\n\n"%mainFunc)

    os.system('chmod 755 script.py')
     

def writeMetaScript(dftcode):
    if dftcode == 'gpaw':
        script = """#!/bin/bash
# Slurm parameters
NTASKS=`echo $SLURM_TASKS_PER_NODE|tr '(' ' '|awk '{print $1}'`
NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`
NCPU=`echo " $NTASKS * $NNODES " | bc`
# load gpaw-specific paths
source /scratch/users/ksb/gpaw/paths.bash
# run parallel gpaw
mpirun -n $NCPU gpaw-python script.py"""

    elif dftcode == 'quantumespresso':
        script = """#!/bin/bash\npython script.py"""

    with open('metascript.sh','w') as f: f.write(script)

    os.system('chmod 755 metascript.sh') 
    

def initialize():
    """
    Does three things to initialize (any) job:
        1. Delete old error and output files for restarted jobs
        2. Load job input parameters
        3. Write initial Atoms object to init.traj
    """
    import json,os,glob
    import ase.io as aseio
    import ase.parallel as asepar

    def keepNewest(string):                 # Step 1
        listPth = glob.glob(string)
        ordered = sorted([(os.path.getmtime(pth),pth) for pth in listPth])
        for t,pth in ordered[:-1]: 
            os.remove(pth)

    if rank()==0:
        keepNewest('*.error')
        keepNewest('*.out')

    try:   os.remove('result.json')
    except OSError: pass
    try:   os.remove('runtime.json')
    except OSError: pass

    with asepar.paropen('params.json','r') as f: prms  = json.loads(f.read())    # Step 2 
    atoms   = makeAtoms(prms)

    aseio.write('init.traj',atoms)         # Step 3

    return prms,atoms

def makeAtoms(params): 
    import cPickle
    return cPickle.loads(str(params['inittraj_pckl']))

def makeCalc(params):
    dftcode = params['dftcode']
    jobkind = params['jobkind']
    relax   = jobkind in ['latticeopt','relax','bulkmod']
    kpt     = makeKPT(params)
    spinpol = makeSpinpol(params)

    def makeGPAWcalc(p):
        from gpaw import GPAW,PW,Davidson,Mixer,MixerSum,FermiDirac
        return GPAW(mode         = PW(p['pw'])                        
                    ,xc          = p['xc']
                    ,kpts        = kpt
                    ,spinpol     = spinpol
                    ,convergence = {'energy':p['econv']} #eV/electron
                    ,mixer       = ((MixerSum(beta=p['mixing'],nmaxold=p['nmix'],weight=100)) 
                                    if spinpol else (Mixer(beta=p['mixing'],nmaxold=p['nmix'],weight=100)))
                    ,maxiter       = p['maxstep']
                    ,nbands        = p['nbands']
                    ,occupations   = FermiDirac(p['sigma'])
                    ,setups        = p['psp']
                    ,eigensolver   = Davidson(5)
                    ,poissonsolver = None 
                    ,txt           ='log'
                    ,symmetry={'do_not_symmetrize_the_density': True}) 

    def makeQEcalc(p):
        from espresso import espresso   

        pspDict =   {'sherlock': {'gbrv15pbe':'/home/vossj/suncat/psp/gbrv1.5pbe'}
                    ,'suncat':   {'gbrv15pbe':'/nfs/slac/g/suncatfs/sw/external/esp-psp/gbrv1.5pbe'}}
        pspPath =  pspDict[getCluster()][params['psp']]

        return espresso( pw         = p['pw']
                        ,dw         = p['pw']*p['dwrat']
                        ,xc         = p['xc']
                        ,kpts       = kpt
                        ,spinpol    = spinpol
                        ,convergence=   {'energy':  p['econv']
                                        ,'mixing':  p['mixing']
                                        ,'nmix':    p['nmix']
                                        ,'maxsteps':p['maxstep']
                                        ,'diag':    'david'}
                        ,nbands     = p['nbands']
                        ,sigma      = p['sigma']
                        ,dipole     = {'status': p['kind'] == 'surface'}
                        ,outdir     = 'calcdir'
                        ,startingwfc= 'atomic+random' 
                        ,psppath    = pspPath
                        ,mode       = 'vc-relax' if jobkind=='vcrelax' else 'scf'
                        ,cell_factor= 2 if jobkind == 'vcrelax' else 1
                        ,output     = {'removesave':True})
    def makeQEVCcalc(p):
        from espresso import espresso
        pspDict =   {'sherlock': {'gbrv15pbe':'/home/vossj/suncat/psp/gbrv1.5pbe'}
                    ,'suncat':   {'gbrv15pbe':'/nfs/slac/g/suncatfs/sw/external/esp-psp/gbrv1.5pbe'}}
        pspPath =  pspDict[getCluster()][params['psp']]
        return  espresso( pw            = p['pw']
                        ,dw             = p['pw']*p['dwrat']
                        ,xc             = p['xc']
                        ,kpts           = kpt
                        ,nbands         = p['nbands']
                        ,dipole         = {'status': False}
                        ,sigma          = p['sigma']
                        ,mode           = 'vc-relax'
                        ,cell_dynamics  = 'bfgs'
                        ,opt_algorithm  = 'bfgs'
                        ,cell_factor    = 2.
                        ,spinpol        = spinpol
                        ,outdir         = 'calcdir'  
                        ,output         = {'removesave':True}
                        ,psppath        = pspPath
                        ,convergence=   {'energy':  p['econv']
                                        ,'mixing':  p['mixing']
                                        ,'nmix':    p['nmix']
                                        ,'maxsteps':p['maxstep']
                                        ,'diag':    'david'})

    def makeQEvibcalc(p):
        from espresso.vibespresso import vibespresso

        pspDict =   {'sherlock': {'gbrv15pbe':'/home/vossj/suncat/psp/gbrv1.5pbe'}
                    ,'suncat':   {'gbrv15pbe':'/nfs/slac/g/suncatfs/sw/external/esp-psp/gbrv1.5pbe'}}
        pspPath =  pspDict[getCluster()][p['psp']]

        return vibespresso( pw          = p['pw']
                            ,dw         = p['pw']*p['dwrat']
                            ,xc         = p['xc']
                            ,kpts       = kpt
                            ,spinpol    = spinpol
                            ,convergence=   {'energy':  p['econv']
                                            ,'mixing':  p['mixing']
                                            ,'nmix':    p['nmix']
                                            ,'maxsteps':p['maxstep']
                                            ,'diag':    'david'}
                            ,nbands     = p['nbands']
                            ,sigma      = p['sigma']
                            ,dipole     = {'status': p['kind'] == 'surface'}
                            ,outdir     = 'calcdir'  
                            ,startingwfc= 'atomic+random' 
                            ,psppath    = pspPath
                            ,output     = {'removesave':True})

    if dftcode =='gpaw':    
        if relax: return makeGPAWcalc(params)
        else: raise NotImplementedError, 'no GPAW calculator-maker for this kind of job'
    elif dftcode =='quantumespresso': 
        if relax:                   return makeQEcalc(params)
        elif jobkind=='vcrelax':    return makeQEVCcalc(params)
        elif jobkind == 'vib':      return makeQEvibcalc(params)
        else: raise NotImplementedError, 'no QE calculator-maker for this kind of job'

def makeKPT(params):    
    """
    Convert k-point density to Monkhorst-Pack grid size. Values forced to be even numbers.
    Special considerations if modeling molecule/bulk/surface. 
    """
    import math  as m
    import numpy as np

    recipcell,kpts = makeAtoms(params).get_reciprocal_cell(),[]
    for i in range(3):
        k = 2 * 3.14159 * m.sqrt((recipcell[i]**2).sum()) * params['kptden'] 
        kpts.append(2 * int(np.ceil(k / 2)))

    kind = params['kind']
    if   kind=='surface':   return np.array(kpts[:2]+[1])
    elif kind=='molecule':  return np.array([1,1,1])
    else:                   return np.array(kpts)

def makeSpinpol(params):
    magmomsinit = makeAtoms(params).get_initial_magnetic_moments()
    return any([x>0 for x in magmomsinit])

def optimizePos(atoms,calc,fmax):
    import ase.optimize as aseopt

    atoms.set_calculator(calc)
    dyn = aseopt.BFGS(atoms=atoms, logfile='qn.log', trajectory='qn.traj',restart='qn.pckl')
    dyn.run(fmax=fmax)

def trajDetails(atoms):
    """ Returns dictionary summary of an (optimized) Atoms object """
    import cPickle
    import numpy as np

    try: mag = atoms.get_magnetic_moments()
    except: mag = np.array([0]*len(atoms))
    return {'finaltraj_pckl':cPickle.dumps(atoms)
            ,'finalpos_pckl':cPickle.dumps(atoms.get_positions())
            ,'finalcell_pckl':cPickle.dumps(atoms.get_cell())
            ,'finalmagmom_pckl':cPickle.dumps(mag)}

def log(params,optatoms):
    import datalog
    datalog.log(optatoms,job_name=params['name'])

############################
############################
# JobKind-specific functions
############################
############################

###########################
# OptimizeLattice functions
###########################

def OptimizeLatticeScript():
    import cPickle,json,os,ase,shutil
    import scipy.optimize as opt
    import ase.parallel   as asepar

    global energies,lparams

    #######################
    print "Initializing..."
    #----------------------
    params,initatoms = initialize()      # Remove old .out/.err files, load from fw_spec, and write 'init.traj'
    shutil.copy('init.traj','out.traj')  # Initial guess for atomic coordinates inside getBulkEnergy reads from 'out.traj'
    energies,lparams = [],[]             # Initialize variables

    if rank()==0:
        for d in ['qn.traj','qn.log','lattice_opt.log','energy_forces.pckl']:          # Remove partially completed calculations that may still be held over from failed job
            if os.path.exists(d): os.remove(d)
            print 'Removed existing file ',d

    ################################
    print "Loading initial guess..."
    #-------------------------------
    
    try:
        with asepar.paropen('lastguess.json','r') as f: iGuess = json.loads(f.read())
        print '\tread lastguess from lastguess.json: ',iGuess
    except: 
        iGuess = getInitGuess(params['structure'],initatoms.get_cell())
        print '\tgenerating initial guess from getInitGuess(): ',iGuess
    
    ########################################
    print "Optimizing cell and positions..."
    #---------------------------------------
    optimizedLatticeParams  = opt.fmin(getBulkEnergy,iGuess,args=(params,),ftol=1,xtol=params['xtol'])
    print 'Optimized lattice params: ',optimizedLatticeParams
    
    ################################
    print "Storing Results..."
    #-------------------------------

    optAtoms    = ase.io.read('out.traj') #read optimized cell and positions
    with asepar.paropen('energy_forces.pckl') as f: eng,forces=cPickle.loads(str(f.read()))
    resultDict  = mergeDicts([params,trajDetails(optAtoms)
                            ,{'raw_energy':     eng
                            ,'forces_pckl':     cPickle.dumps(forces)
                            ,'latticeopt_pckl': cPickle.dumps(zip(energies,lparams))}])  
    
    with open('result.json', 'w') as f: f.write(json.dumps(resultDict))
    if rank()==0:
        log(params,optAtoms)
        #with open('result.json', 'r') as outfile: json.load(outfile) #WHY DOES THIS CRASH
    return 0

def getInitGuess(structure,cell):
    """
    Convert cell into parameters a,b,c,alpha,beta,gamma
    Depending on structure, return sufficient information required to reconstruct cell (e.g. only 'a' is needed if the structure is known cubic)
    """
    import numpy as np
    
    def angle(v1,v2): return np.arccos(np.dot(v1,np.transpose(v2))/(np.linalg.norm(v1)*np.linalg.norm(v2)))
    a,b,c            = np.linalg.norm(cell[0]),np.linalg.norm(cell[1]),np.linalg.norm(cell[2])
    alpha,beta,gamma = angle(cell[1],cell[2]), angle(cell[0],cell[2]), angle(cell[0],cell[1])

    if   structure in ['fcc','bcc','rocksalt','diamond','cscl','zincblende']: return [a]
    elif structure in ['hexagonal']: return [a,c]
    elif structure in ['triclinic']: return [a,b,c,alpha,beta,gamma]
    else: raise ValueError, 'Bad entry in "structure" field for Atoms object info dictionary: '+structure

def getBulkEnergy(latticeParams,params):
    #For a given set of bravais lattice parameters, optimize atomic coordinates and return minimum energy
    import cPickle,ase,json
    import ase.parallel as asepar 
    import ase.io as aseio
    global energies,lparams

    with asepar.paropen('lastguess.json','w') as f: f.write(json.dumps(list(latticeParams)))
    atomsInitUnscaled = makeAtoms(params)
    atoms = fromParams(atomsInitUnscaled,latticeParams,params['structure'])
    optimizePos(atoms,makeCalc(params),params['fmax'])
    energy,forces = atoms.get_potential_energy(),atoms.get_forces()
    energies.append(energy);lparams.append(latticeParams)
    with asepar.paropen('lattice_opt.log','a') as logfile: logfile.write('%s\t%s\n' %(energy,latticeParams))
    aseio.write('out.traj',atoms)
    with open('energy_forces.pckl','w') as f: f.write(cPickle.dumps((energy,forces)))
    return energy

def fromParams(atomsInput,cellParams,structure): 
    """
    Params is a list of 1 to 6 numbers (a,b,c,alpha,beta,gamma). 
    ANGLES OF INPUT ARE IN RADIANS, ASE CELLS WANT ANGLES IN DEGREES
    Depending on structure, we can construct the cell from a subset of these parameters
    """
    import math     as m
    
    a = cellParams[0]
    if   structure in ['cubic','cscl']:                          cell = [a,a,a,90,90,90]
    elif structure in ['fcc','diamond','zincblende','rocksalt']: cell = [a,a,a,60,60,60]
    elif structure in ['bcc']:                                   cell = [a,a,a,109.47122,109.47122,109.47122]
    elif structure in ['hcp','hexagonal']:                       cell = [a,a,cellParams[1],90,90,120]                                                                                          # Input is assumed to be two parameters, a and c
    elif structure in ['triclinic']:                             cell = [cellParams[0],cellParams[1],cellParams[2],m.degrees(cellParams[3]),m.degrees(cellParams[4]),m.degrees(cellParams[5])] # You should really be using VC relax for this....
    else: raise NotImplementedError, 'fromParams(atomsInput,cellParams,structure) cannot handle unknown structure = '+structure
    
    atoms = atomsInput.copy() 
    atoms.set_cell(cell,scale_atoms=True)

    return atoms

###########################
# BulkModulus functions
###########################
def BulkModulusScript():
    import ase,json,base64,copy,os,matplotlib
    import numpy        as np
    import ase.eos      as aseeos
    import ase.parallel as asepar 
    matplotlib.use('Agg')
    #######################
    print "Initializing..."
    #----------------------

    params,optAtoms = initialize()  # Remove old .out/.err files, load from fw_spec, and write 'init.traj'

    optCell,vol,eng = optAtoms.get_cell() , [],[]
    strains         = np.linspace(1 - params['strain'],1 + params['strain'],9)      
    
    if rank()==0:
        for d in ['qn.traj','qn.log']:          # Remove partially completed calculations that may still be held over from failed job
            if os.path.exists(d): os.remove(d)
            print 'Removed existing file ',d

    ##################################################
    print "Calculating energies at various strains..."
    #-------------------------------------------------

    for i, strain in enumerate(strains):
        atoms = optAtoms.copy()
        atoms.set_cell(optCell*strain,scale_atoms=True)
        optimizePos(atoms,makeCalc(params),params['fmax'])

        volume,energy = atoms.get_volume(),atoms.get_potential_energy()
        vol.append(copy.deepcopy(volume));eng.append(copy.deepcopy(energy))

    ################################
    print "Analyzing dE/dV curve..."
    #-------------------------------

    aHat,quadR2 = quadFit(np.array(copy.deepcopy(vol)),np.array(copy.deepcopy(eng)))

    try:        
        eos = aseeos.EquationOfState(vol,eng)
        v0, e0, b = eos.fit()
        eos.plot(filename='bulk-eos.png',show=False)
        b0= b/ase.units.kJ*1e24                                     #GPa: use this value if EOS doesn't fail
        with open('bulk-eos.png', 'rb') as f: img = base64.b64encode(f.read())

    except ValueError:                                  # too bad of a fit for ASE to handle
        b0 = aHat*2*vol[4]*160.2                        # units: eV/A^6 * A^3 * 1, where 1 === 160.2 GPa*A^3/eV
        img = None

    ################################
    print "Storing Results..."
    #-------------------------------

    resultDict  = mergeDicts([params,trajDetails(optAtoms)
                                ,{'bulkmod':b0,'bfit':quadR2,'bulkmodimg_base64':img,'voleng_json':json.dumps(zip(vol,eng))}])

    with open('result.json', 'w') as outfile:   outfile.write(json.dumps(resultDict))
    
    if rank()==0: log(params,optAtoms)
    return 0

def quadFit(xIn,yIn):
    """ Input x vector: units A^3, Input y vector: units eV """
    import numpy            as np
    import scipy.optimize   as opt

    # Center data around 4th data point
    x = xIn-xIn[4]; y = yIn-yIn[4]
    def model(a):  return a*np.square(x)
    def errVec(a): return a*np.square(x) - y   #create fitting function of form mx+b
    aHat, success = opt.leastsq(errVec, [0.1])
    yhat    = model(aHat)
    ybar    = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssresid = np.sum((yhat-y)**2)  
    sstotal = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    r2      = 1 - (ssresid / float(sstotal)) # 
    return float(aHat[0]),float(r2)

###########################
# XC Contribs functions
##########################
def XCcontribsScript():
    import cPickle,json,ase
    from gpaw        import restart
    from gpaw.xc.bee import BEEFEnsemble
    params,optAtoms = initialize()  # Remove old .out/.err files, load from fw_spec, and write 'init.traj'
    pbeParams = copy.deepcopy(params)
    pbeParams['xc']='PBE'

    atoms.set_calculator(makeCalc(pbeParams))
    atoms.get_potential_energy()
    atoms.calc.write('inp.gpw', mode='all') 

    atoms,calc = restart('inp.gpw', setups='sg15', xc='mBEEF', convergence={'energy': 5e-4}, txt='mbeef.txt')
    atoms.get_potential_energy()
    beef        = BEEFEnsemble(calc)
    xcContribs  = beef.mbeef_exchange_energy_contribs()
    ################################
    print "Storing Results..."
    #-------------------------------
    ase.io.write('final.traj',atoms)
    resultDict  = mergeDicts([params,trajDetails(ase.io.read('final.traj'))
                                ,{'xcContribs': cPickle.dumps(xcContribs)}])
    with open('result.json', 'w') as outfile:   outfile.write(json.dumps(resultDict))
    if rank()==0:
        with open('result.json', 'r') as outfile: json.loads(outfile.read()) #test that dictionary isn't 'corrupted'
        log(params,atoms)
    return 0

#################
# Relax functions
#################

def RelaxScript():
    import ase,json,cPickle,os
    import ase.parallel   as asepar

    #######################
    print "Initializing..."
    #----------------------

    params,initatoms = initialize()  # Remove old .out/.err files, load from fw_spec, and write 'init.traj'

    if rank()==0:
        if os.path.exists('qn.traj'): 
            if os.stat('qn.traj').st_size > 100:    
                initatoms = ase.io.read('qn.traj')
                print '\treading from qn.traj...'
            else:                                   
                os.remove('qn.traj')
                print 'removed empty qn.traj'
        if os.path.exists('qn.log') and os.stat('qn.log').st_size < 100: 
            os.remove('qn.log')
            print '\tremoved empty qn.log...'

    #######################
    print "Optimizing positions..."
    #----------------------
    optimizePos(initatoms,makeCalc(params),params['fmax'])

    ############################
    print "Storing Results..."
    #--------------------------
    ase.io.write('final.traj',initatoms)    
    e0 = initatoms.get_potential_energy()
    f0 = initatoms.get_forces()
    optAtoms = ase.io.read('final.traj')
    resultDict  = mergeDicts([params,trajDetails(optAtoms)
                                ,{'raw_energy': e0
                                ,'forces_pckl':cPickle.dumps(f0)} ])

    with open('result.json', 'w') as f: f.write(json.dumps(resultDict))
    if rank()==0:

        # For some unknown reason, when trying to pickle dump initatoms in log(), get error (cannot dump open file)
        log(params,optAtoms)
    return 0


#####################
# Vibration Functions
#####################
def VibScript():
    from ase.vibrations       import Vibrations
    import glob,json,cPickle,ase,os
    
    #######################
    print "Initializing..."
    #----------------------
    
    params,atoms = initialize()  # Remove old .out/.err files, load from fw_spec, and write 'init.traj'
    
    prev = glob.glob('*.pckl') #delete incomplete pckls - facilitating restarted jobs
    
    if rank()==0:
        for p in prev:
            if os.stat(p).st_size < 100: os.remove(p)
    
    atoms.set_calculator(makeCalc(params))
    vib = Vibrations(atoms,delta=params['delta'],indices=json.loads(params['vibids_json']))
    
    vib.run(); vib.write_jmol()
    
    ##########################
    print "Storing Results..."
    #-------------------------
    
    vib.summary(log='vibrations.txt')
    
    with open('vibrations.txt','r') as f: vibsummary = f.read()

    ase.io.write('final.traj',atoms)    
    optAtoms = ase.io.read('final.traj')

    vib_energies,vib_frequencies = vib.get_energies(),vib.get_frequencies()

    resultDict  = mergeDicts([params,  trajDetails(optAtoms),{'vibfreqs_pckl': cPickle.dumps(vib_frequencies)
                                        ,'vibsummary':vibsummary
                                        ,'vibengs_pckl':cPickle.dumps(vib_energies)}])

    with open('result.json', 'w') as outfile:   outfile.write(json.dumps(resultDict))
    with open('result.json', 'r') as outfile: json.loads(outfile.read()) #test that dictionary isn't 'corrupted'
    if rank()==0: log(params,optAtoms)
    return 0


#####################
# VCRelax Functions
#####################
def VCRelaxScript():
    import ase,json,cPickle,os
    #######################
    print "Initializing..."
    #----------------------
    params,atoms = initialize()  # Remove old .out/.err files, load from fw_spec, and write 'init.traj'
    
    if not os.path.exists('intermediate.traj'):
        ##########################################
        print "Running VC Relax for first time..."
        #-----------------------------------------
        atoms.set_calculator(makeCalc(params))
        energy = atoms.get_potential_energy()   # Trigger espresso to be launched
        ase.io.write('intermediate.traj',atoms.calc.get_final_structure())

    ###########################################
    print "Running VC Relax for second time..."
    #------------------------------------------

    atoms = ase.io.read('intermediate.traj')
    atoms.set_calculator(makeCalc(params))
    energy = atoms.get_potential_energy() #trigger espresso to be launched
    ase.io.write('final.traj',atoms.calc.get_final_structure())
 
    ################################
    print "Storing Results..."
    #-------------------------------
    e0,f0   = atoms.get_potential_energy(),atoms.get_forces()
    atoms   = ase.io.read('final.traj')

    resultDict  = mergeDicts([params,trajDetails(atoms),
                                {'raw_energy':  e0
                                ,'forces_pckl': cPickle.dumps(f0)}])

    with open('result.json', 'w') as outfile:   outfile.write(json.dumps(resultDict))
    with open('result.json', 'r') as outfile: json.loads(outfile.read()) #test that dictionary isn't 'corrupted'
    log(params,atoms)
    return 0

#####################
# Not Implemented
#####################

def DosScript(): raise NotImplementedError
def NebScript(): raise NotImplementedError

#############
#############
# FIRETASKS
#############
#############


def executeAndCollect():
    os.system('./metascript.sh')
    try:
        with open('result.json','r') as f: resultDictjson = f.read()
    except IOError:
        errFile = glob.glob('*.error')[0]
        with open(errFile,'r') as f: raise Exception, f.read()
        
    return fireworks.core.firework.FWAction(stored_data={'resultjson':resultDictjson})

@fireworks.utilities.fw_utilities.explicit_serialize
class OptimizeLattice(fireworks.core.firework.FiretaskBase):
    def run_task(self,fw_spec):

        params = fw_spec['params']

        with open('params.json','w') as f: f.write(json.dumps(params))

        writeScript([initialize  # Universal        
                    ,makeAtoms
                    ,makeCalc
                    ,makeKPT
                    ,makeSpinpol
                    ,optimizePos
                    ,trajDetails
                    ,log
                    ,OptimizeLatticeScript # Jobkind-specific
                    ,getInitGuess
                    ,getBulkEnergy
                    ,fromParams
                    ,rank
                    ,mergeDicts # Auxillary
                    ,getCluster
                    ],'OptimizeLatticeScript')

        writeMetaScript(params['dftcode'])

        return executeAndCollect()



@fireworks.utilities.fw_utilities.explicit_serialize
class GetBulkModulus(fireworks.core.firework.FiretaskBase):
    def run_task(self,fw_spec): 
        

        params = fw_spec['params']

        with open('params.json','w') as f: f.write(json.dumps(params))

        writeScript([initialize  # Universal        
                    ,makeAtoms
                    ,makeCalc
                    ,makeKPT
                    ,makeSpinpol
                    ,optimizePos
                    ,trajDetails
                    ,log
                    ,BulkModulusScript # Jobkind-specific
                    ,quadFit
                    ,rank
                    ,mergeDicts # Auxillary
                    ,getCluster
                    ],'BulkModulusScript')

        writeMetaScript(params['dftcode'])

        os.system('./metascript.sh')

        with open('result.json','r') as f: resultDictjson = f.read()

        return fireworks.core.firework.FWAction(stored_data={'resultjson':resultDictjson})

@fireworks.utilities.fw_utilities.explicit_serialize
class GetXCcontribs(fireworks.core.firework.FiretaskBase):
    def run_task(self,fw_spec): 

        params = fw_spec['params']

        with open('params.json','w') as f: f.write(json.dumps(params))

        writeScript([initialize         # Universal     
                    ,makeAtoms
                    ,makeCalc
                    ,makeKPT
                    ,makeSpinpol
                    ,optimizePos
                    ,trajDetails
                    ,log
                    ,XCcontribsScript   # Jobkind-specific
                    ,mergeDicts         # Auxillary
                    ,rank
                    ,getCluster
                    ],'XCcontribsScript')

        writeMetaScript(params['dftcode'])

        return executeAndCollect()

#################################################################
@fireworks.utilities.fw_utilities.explicit_serialize
class Relax(fireworks.core.firework.FiretaskBase):
    def run_task(self,fw_spec): 

        params = fw_spec['params']

        with open('params.json','w') as f: f.write(json.dumps(params))

        writeScript([initialize    # Universal     
                    ,makeAtoms
                    ,makeCalc
                    ,makeKPT
                    ,makeSpinpol
                    ,optimizePos
                    ,log
                    ,trajDetails
                    ,RelaxScript    # Jobkind-specific
                    ,rank           # Auxillary
                    ,mergeDicts         
                    ,getCluster
                    ],'RelaxScript')

        writeMetaScript(params['dftcode'])

        return executeAndCollect()

@fireworks.utilities.fw_utilities.explicit_serialize
class Vibrations(fireworks.core.firework.FiretaskBase):
    def run_task(self,fw_spec): 

        params = fw_spec['params']

        with open('params.json','w') as f: f.write(json.dumps(params))

        writeScript([initialize         # Universal     
                    ,makeAtoms
                    ,makeCalc
                    ,makeKPT
                    ,makeSpinpol
                    ,optimizePos
                    ,trajDetails
                    ,log
                    ,VibScript  # Jobkind-specific
                    ,rank
                    ,mergeDicts         # Auxillary
                    ,getCluster
                    ],'VibScript')

        writeMetaScript(params['dftcode'])
        
        return executeAndCollect()
        
@fireworks.utilities.fw_utilities.explicit_serialize
class DOS(fireworks.core.firework.FiretaskBase):
    def run_task(self,fw_spec):

        params = fw_spec['params']

        with open('params.json','w') as f: f.write(json.dumps(params))

        writeScript([initialize         # Universal     
                    ,makeAtoms
                    ,makeCalc
                    ,makeKPT
                    ,makeSpinpol
                    ,optimizePos
                    ,trajDetails
                    ,log
                    ,DosScript          # Jobkind-specific
                    ,rank
                    ,mergeDicts         # Auxillary
                    ,getCluster
                    ],'DosScript')

        writeMetaScript(params['dftcode'])
        
        return executeAndCollect()

@fireworks.utilities.fw_utilities.explicit_serialize
class NEB(fireworks.core.firework.FiretaskBase):
    def run_task(self,fw_spec):

        params = fw_spec['params']

        with open('params.json','w') as f: f.write(json.dumps(params))

        writeScript([initialize         # Universal     
                    ,makeAtoms
                    ,makeCalc
                    ,makeKPT
                    ,makeSpinpol
                    ,trajDetails
                    ,log
                    ,NebScript          # Jobkind-specific
                    ,rank
                    ,mergeDicts         # Auxillary
                    ,getCluster
                    ],'NebScript')

        writeMetaScript(params['dftcode'])

        return executeAndCollect()

@fireworks.utilities.fw_utilities.explicit_serialize
class VCRelax(fireworks.core.firework.FiretaskBase):
    def run_task(self,fw_spec):

        params = fw_spec['params']

        with open('params.json','w') as f: f.write(json.dumps(params))

        writeScript([initialize         # Universal     
                    ,makeAtoms
                    ,makeCalc
                    ,makeKPT
                    ,makeSpinpol
                    ,optimizePos
                    ,trajDetails
                    ,log
                    ,VCRelaxScript      # Jobkind-specific
                    ,rank
                    ,mergeDicts         # Auxillary
                    ,getCluster
                    ],'VCRelaxScript')

        writeMetaScript(params['dftcode'])

        return executeAndCollect()

