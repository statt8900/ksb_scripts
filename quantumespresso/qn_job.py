#/usr/bin/env python
#LSF -q suncat2
#LSF -n 12
#LSF -W 30:00
#LSF -o opt.log
#LSF -e err.log
#LSF -sp 100
#LSF -N

from espresso import espresso
from ase.parallel import *
from ase import io
from ase.atoms import np
import os
import time
import shutil
from ase import *
from ase.structure import molecule
from ase.dft.bee import BEEF_Ensemble
import numpy as np
import os

###################################
"""INPUTS:"""
###################################

xc = 'BEEF'
kpts = (4, 4, 1)
nbands = -10
pw_cutoff = 500
dw_cutoff = 5000
spinpol = False
outdir = 'esp.log'

dipole = {'status':False}

output      = {'avoidio':False,
           'removewf':True,
           'wf_collect':False}
convergence = {'mixing': 0.1,
               'maxsteps': 500,
               'nmix':10,
               'energy':5e-4,
               'diag':'david'}
###################################
###################################




calc = espresso(pw=pw_cutoff,
    dw = dw_cutoff,
    beefensemble=True, # important!
    kpts = kpts,
    nbands = nbands,
    xc = xc,
    convergence = convergence,
    dipole = dipole,
    spinpol = spinpol,
    outdir = outdir,
    output = output,
#    mode = 'scf',
    )


#!/usr/bin/env python



################################################################################
#      Class definitions. Preferred to import to improve reproducability.      #
################################################################################

class Calculation:
    def __init__(self,path,**kwargs):
        self.path = path
        self.name = 'atoms'
        self.atoms_suffix = 'traj'
        self.converged_file = 'converged.log'
        self.keep_files = [] #files which will not be moved
        self.keep_time = 30 #do not delete any files that have not existed for this many seconds (protects logs)
        self.keep_suffixs = [] #do not remove any files with these suffixes
        self.coarse_grain_dict = {} #Dictionary of arguments to use for coarse-graining calculations

        self.task_prefix = None #task prefix; must be set
        self.calculator = None #calculator to use; must be set
        self.atoms = None #atoms object to optimize; must be set

        for key in kwargs:
            setattr(self,key,kwargs[key])

    def get_last_run_number(self):
        #get the number of the last run. Return -1 if no runs have occured
        last = -1
        files = os.listdir(self.path)
        runs = [f for f in files if f.startswith(self.task_prefix)]
        ns = []
        for r in runs:
            try:
                n = int(r.replace(self.task_prefix,''))
                ns.append(n)
            except:
                pass
        last = max(ns+[last])
        return last

    def get_next_run_dir(self):
        last = self.get_last_run_number()
        if last > -1:
            dirname = self.task_prefix+str(last+1)
        else:
            dirname = self.task_prefix+'0'
        return os.path.join(self.path,dirname)

    def prepare_atoms(self):
        self.atoms.calc = self.calculator
        self.atoms.cell = np.round(self.atoms.cell,12) #Ensure cell is not numerically non-orthogonal

    def prepare_calculator(self):
        if not self.coarse_grain_dict:
            self.coarse_grained = False
            return
        else:
            self.coarse_grained = True
            #if coarse_grain_dict is defined then modify the calculator appropriately.
            self.check_coarse_grain_dict()
            params,cgiter = self.get_coarse_grain_params()
            if params:
                task_params = {}
                for p in params.keys():
                    #If the parameters belong to the task, remove them from the params
                    if p in self.task_parameters:
                        task_params[p] = params[p]
                        del params[p]
                self.calculator.set(**params)
                cg_file = 'coarse_grain_'+str(cgiter)
                if rank == 0: #remake the file even if it exists so that it is protected by "keep_time"
                    f = open(cg_file,'w')
                    f.write('parameters = '+str(params))
                    f.close()
                self.coarse_grain_converged_file = 'coarse_grain_'+str(cgiter)+'_converged.log'
            else:
                self.coarse_grained = False

    def check_coarse_grain_dict(self):
        cgiters = max([len(self.coarse_grain_dict[k]) for k in self.coarse_grain_dict])
        for key in self.coarse_grain_dict:
            plen = len(self.coarse_grain_dict[key])
            if plen < cgiters:
                if plen == 1:
                    self.coarse_grain_dict[key] = self.coarse_grain_dict[key]*cgiters
                else:
                    raise ValueError('Length of coarse_grain_dict values must be the same!')

    def get_coarse_grain_params(self):
        files = os.listdir(self.path)
        cg_iter = 0
        cg_files = [int(fi.replace('coarse_grain_','')) for fi in files if (fi.startswith('coarse_grain') and 'converged' not in fi)]
        cg_converged_files = []
        for path,dirs,files in os.walk(self.path):
            cg_converged_files += [int(fi.replace('coarse_grain_','').replace('_converged.log','')) for fi in files if (fi.startswith('coarse_grain') and 'converged' in fi)]

        if len(cg_files) > 1:
            raise ValueError('Only one coarse_grain file should exist! Progress is ambiguous')
        if cg_files: cg_iter = cg_files[0]
        if cg_converged_files: cg_iter = max(cg_converged_files)+1
        param_dict = {}
        for key in self.coarse_grain_dict:
            if  cg_iter < len(self.coarse_grain_dict):
                param_dict[key] = self.coarse_grain_dict[key][cg_iter]
            else: #all coarse graining is finished
                param_dict = None
        return param_dict, cg_iter

    def make_converged_file(self):
        if rank == 0:
            if self.coarse_grained:
                converged_file = self.coarse_grain_converged_file
            else:
                converged_file = self.converged_file

            if hasattr(self,'converged_text'):
                f = open(converged_file,'w')
                f.write(self.converged_text)
                f.close()
            else:
                raise AttributeError('Calculation is not converged!')

class GeometryOptimization(Calculation):
    def __init__(self,*args,**kwargs):
        self.optimizer = 'QuasiNewton'
        self.task_parameters = {'fmax':0.05}
        self.optimizer_kwargs = {}
        Calculation.__init__(self,*args,**kwargs)
        self.task_prefix = 'qn'
    def get_all_atoms(self,sort_mode = 'energy'):
        #go through previous calculations and find ASE atoms object
        energy_atoms = []
        empty_atoms = []
        for path,dirs,files in os.walk(self.path):
            for f in [fi for fi in files if fi.endswith('.'+self.atoms_suffix)]:
                atompath = os.path.join(path,f)
                name = f.replace('.'+self.atoms_suffix,'')
                ctime = os.path.getctime(atompath) #creation time
                try:
                    atoms = io.read(atompath,':')
                    try:
                        energy = atoms[-1].get_potential_energy()
                        energy_atoms.append([ctime,name,atoms])
                    except:
                        empty_atoms.append([ctime,name,atoms])
                except:
                    pass

        if sort_mode == 'energy':
            if not energy_atoms:
                sort_mode = 'time'
                print 'No atoms had energy. Sorting by time instead'
            else:
                cts,names,all_atoms = zip(*energy_atoms)
                all_atoms = list(all_atoms)
                names = list(names)
                all_atom_obj = []
                all_names = []
                for atm,name in zip(all_atoms,names):
                    all_atom_obj += atm
                    all_names += [name]*len(atm)
                all_atoms = all_atom_obj
                energy_atoms = []
                for a,name in zip(all_atoms,all_names):
                    energy_atoms.append([a.get_potential_energy(),name,a])
                energy_atoms.sort()
                energy_atoms.reverse()
                es,names,atoms = zip(*energy_atoms)
                self.name = names[-1]
                return list(atoms)

        if sort_mode == 'time':
            all_atoms = energy_atoms + empty_atoms
            all_atoms.sort()
            ct,name,atoms = all_atoms[0]
            self.name = name
            return atoms

    def get_atoms(self,*args,**kwargs):
        all_atoms = self.get_all_atoms(*args,**kwargs)
        self.atoms =  all_atoms[-1]
        return self.atoms

    def prepare_directory(self):
        #move all old output to proper subdir; only files should be those used to run the job.

        def keep(filename):
            #function to determine whether or not to keep a file based on its path
            keep_file = False
            fbase = os.path.basename(filename)
            age = time.time() - os.path.getctime(filename)
            if fbase == os.path.basename(__file__):
                keep_file = True
            if filename.rsplit('.',1)[-1] in self.keep_suffixs:
                keep_file = True
            if filename in self.keep_files:
                keep_file = True
            if age < self.keep_time: #HACKish
                keep_file = True
            if fbase.startswith(self.task_prefix): #check to see if its a run directory
                if os.path.isdir(filename):
                    n = fbase.replace(self.task_prefix,'')
                    try:
                        float(n)
                        keep_file = True
                    except:
                        pass
            return keep_file

        if rank == 0:
            next_dir = self.get_next_run_dir()
            assert not os.path.exists(next_dir) #directory should not ever exist
            os.mkdir(next_dir)
            for f in os.listdir(self.path):
                fpath = os.path.join(self.path,f)
                if not keep(fpath):
                    if os.path.isdir(fpath):
                        dest = os.path.join(next_dir,f)
                        shutil.move(fpath,dest)
                    else:
                        shutil.move(fpath,next_dir)

    def run(self):
        barrier()
        import_str = 'from ase.optimize import '+self.optimizer+' as QN'
        exec(import_str)
        kwargs = dict(trajectory=self.name+'.traj',logfile=self.name+'_qn.log',restart=self.name+'_qn.pkl')
        kwargs.update(self.optimizer_kwargs)
        dyn = QN(self.atoms,**kwargs)
        dyn.run(**self.task_parameters)
        barrier()
        opt.converged_text = str(self.atoms.get_potential_energy())
        opt.make_converged_file()


################################################################################
#     End class definitions. Search EOCD to skip to here                       #
################################################################################

opt = GeometryOptimization('.')
opt.calculator = calc
opt.get_atoms()
opt.prepare_atoms()
opt.prepare_calculator()
opt.prepare_directory()
opt.run()

 # ensemble
ens = BEEF_Ensemble(calc)
ens.get_ensemble_energies()
ens.write('qn.bee')

del sys, calc, ens