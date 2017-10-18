import sys
from details_ksb import *
import databaseFuncs as db
###########################################################################################################
###########################################################################################################
###########################################################################################################

class Detail(object):
    def __init__(self,colname,datatype,inputcols,func,author,description,source):
        self.colname        = colname            # Column name in database
        self.datatype       = datatype           # SQL datatype
        self.inputcols      = inputcols          # Dependencies on other columns
        self.func           = func               # Function applied to (pre-existing) columns to get result. Return None if undefined, otherwise should be a total function
        self.author         = author             # Author of the Detail
        self.description    = description        # Description for others to read to know what the column represents / why it is important
        self.source         = source             # Path to source code of function
        self.axislabel      = colname #just for kris, temporarily
    def apply(self,rowid,verbose=True):
        inputvals = self.query(rowid)
        if verbose: print 'in apply %s with inputvals = '%self.colname,inputvals
        try:
            output = self.func(*inputvals)
            if output is None: db.sqlexecute('update derived SET {0} = ?,{0}_undefined=?,{0}_failed = ? where derived_job = ?'.format(self.colname),[None,1,'',rowid]) 
            else: db.sqlexecute('update derived SET {0} = ?,{0}_undefined=?,{0}_failed = ? where derived_job = ?'.format(self.colname),[output,0,'',rowid]) 
        except Exception as e:
            failedOutput = self.colname+' failed with inputvals '+str(inputvals)+' and exception '+str(e)
            if verbose: print failedOutput
            db.sqlexecute('update derived SET {0} = ?,{0}_undefined=?,{0}_failed = ? where derived_job = ?'.format(self.colname),[None,0,failedOutput,rowid]) 
    
    def query(self,rowid): return db.query(self.inputcols,'jobid = %d'%rowid)[0]
    
    def document(self):
        return '\n'.join([self.colname,self.datatype,self.author,self.description,self.source])
    
    def addCol(self):
        """ Add column and supporting info columns """
        try:     db.addCol(self.colname,          self.datatype,'derived')
        except: pass
        try:     db.addCol(self.colname+'_failed',  'varchar', 'derived')
        except: pass
        try:    db.addCol(self.colname+'_undefined','integer', 'derived')
        except: pass

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

                    

def sortDetails(ds):
    """
    Order Detail objects in a way such that one's dependency occurs after itself
    """
    print "Sorting details ..."
    dlist = [(d.colname,i+1) for i,d in enumerate(ds)]
    dmap = dict(dlist)
    maxdep = {d.colname:max([dmap[xx] if xx in dmap.keys() else 0 for xx in d.inputcols]) for d in ds} # max index of dependencies
    triple = [(d,dmap[d.colname],maxdep[d.colname]) for d in ds]
    details = [t for t in triple if t[2] == 0] #initialize output details
    lenIn,lenOut=1,0
    while lenIn!=lenOut:
        currentIDs = [t[1] for t in details]
        lenIn = len(details)
        for t in triple:
            if t[2] in currentIDs and t[1] not in currentIDs: details.append(t)
        lenOut = len(details)
    if lenOut != len(ds): raise ValueError, 'One or more details has a dependency not present in input to sortDetails'
    return [d[0] for d in details]


# {SQLQUERY : VALUE}


###########################################################################################################

allDetails = [
    ## Only applies to KSB jobs
    # input parameters
    Detail('job_type_ksb',         'varchar',['user' ,'storage_directory'],   getJobtypeKSB,    'ksb','Looks for ksb specific tag to detect latticeopt,relax,vcrelax,bulkmod,dos,neb,etc','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('kptden_ksb',           'numeric',['user' ,'storage_directory'],   getKptdenKSB,    'ksb','Gets kpt density from result.json','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('system_type_ksb',      'varchar',['user' ,'storage_directory'],   getKindKSB,      'ksb','Checks input parameters whether bulk,surface, or molecule','/home/ksb/scripts/datalog/details_ksb.py')

    ,Detail('structure_ksb',       'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getStructureKSB, 'ksb','Gets crystal structure (fcc,hcp,etc.) from result.json','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('bulkvacancy_ksb',     'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getBvacKSB, 'ksb','Gets crystal structure (fcc,hcp,etc.) from result.json','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('bulkscale_ksb',       'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getBscaleKSB, 'ksb','Gets crystal structure (fcc,hcp,etc.) from result.json','/home/ksb/scripts/datalog/details_ksb.py')

    ,Detail('facet_ksb',       'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getFacetKSB,      'ksb','Gets facet from result.json','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('xy_ksb',          'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getXyKSB,         'ksb','Gets scale of surface (x by y) from result.json','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('layers_ksb',      'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getLayersKSB,     'ksb','Gets # of layers from result.json','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('constrained_ksb', 'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getConstrainedKSB,'ksb','Gets # of constrained layers from result.json','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('symmetric_ksb',   'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getSymmetricKSB,  'ksb','Gets whether or not surface is symmetric from result.json','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('vacuum_ksb',      'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getVacuumKSB,     'ksb','Gets vacuum from result.json','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('vacancies_ksb',   'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getVacanciesKSB,  'ksb','Gets indices of deleted atoms from result.json','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('adsorbates_ksb',  'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getAdsorbatesKSB, 'ksb','Gets adsorbate json dictionary from result.json','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('sites_ksb',       'varchar',['user' ,'system_type_ksb' ,'storage_directory'],    getSitesKSB, 'ksb','Gets site image (base64) from result.json','/home/ksb/scripts/datalog/details_ksb.py')

    ,Detail('bulkparent_ksb',      'integer',['user' ,'job_type_ksb' 
                                             ,'system_type_ksb' ,'storage_directory'],          getBulkparentKSB, 'ksb','Extracts bulkparent from parameter dictionary','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('parent_ksb',          'integer',['user' ,'job_type_ksb' 
                                             ,'storage_directory'],          getParentKSB, 'ksb','Extracts parent from parameter dictionary','/home/ksb/scripts/datalog/details_ksb.py')

    
    # other
    ,Detail('kptden_x',          'numeric',['storage_directory'
                                           ,'kpts_json'],               kx,             'ksb','kx','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('kptden_y',          'numeric',['storage_directory'
                                           ,'kpts_json'],               ky,             'ksb','ky','/home/ksb/scripts/datalog/details_ksb.py')
    ,Detail('kptden_z',          'numeric',['storage_directory'
                                           ,'kpts_json'],               kz,             'ksb','kz','/home/ksb/scripts/datalog/details_ksb.py')

    ,Detail('raw_energy',          'numeric',['storage_directory'],      getEng,        'ksb','Reads from raw.pckl','/home/ksb/scripts/datalog/details_ksb.py')

    ,Detail('fwid',                'integer',['user'
                                            ,'storage_directory'],      getFWID,         'ksb','Looks for fwid in FW_submit.py','/home/ksb/scripts/datalog/details_ksb.py')

   ,Detail('strjob',              'varchar',['user'
                                            ,'storage_directory'],      makeStrJob,      'ksb','Converts a job into a string that can be compared to other jobs for equivalence','/home/ksb/scripts/datalog/details_ksb.py')

   ,Detail('bulk_modulus_quadfit','numeric',['user'
                                            ,'job_type_ksb'
                                            ,'storage_directory'],      getQuadfitKSB,    'ksb','Looks in result.json and gets quadratic fit','/home/ksb/scripts/datalog/details_ksb.py')
    
   ,Detail('bulk_modulus',         'numeric',['user'
                                            ,'job_type_ksb'
                                            ,'storage_directory'],      getBulkmodKSB,    'ksb','Looks in result.json and gets quadratic fit','/home/ksb/scripts/datalog/details_ksb.py')

    ,Detail('error_lattice_A',     'numeric',['lattice_parameter'
                                             ,'keld_lattice_A'
                                             ,'structure_ksb'],         errAFunc,     'ksb','Error in Lattice Parameter, Angstrom','/home/ksb/scripts/datalog/details_ksb.py')

    ,Detail('error_BM',            'numeric',['bulk_modulus'
                                             ,'keld_lattice_BM'],       errBM,        'ksb','Error in Lattice Parameter, Angstrom','/home/ksb/scripts/datalog/details_ksb.py')


   ,Detail('fwstatus',            'varchar',['fwid'],                   getFWstatus,     'ksb','Current status of firework','/home/ksb/scripts/datalog/details_ksb.py')

   ,Detail('timeperstep',         'numeric',['user','job_type_ksb'
                                            ,'storage_directory'
                                            ,'fwid'],               getTPS,                'ksb','Gets time per step for a latticeopt or relaxation','/home/ksb/scripts/datalog/details_ksb.py')

   ,Detail('pawjob',              'integer',['numbers_str'],        pawElem,               'ksb','1 if there are any alkali/alkali earth metals or zinc present else 0','/home/ksb/scripts/datalog/details_ksb.py')

    ## Applies to all jobs
    ,Detail('number_of_atoms',     'integer',['storage_directory'],    getNatoms,       'ksb','Counts number of atoms','/home/ksb/scripts/datalog/details_ksb.py')

    ,Detail('psppath',             'varchar',['dftcode'  
                                             ,'storage_directory'],    getPSPpath,      'ksb','Extracts pseudopotential path','/home/ksb/scripts/datalog/details_ksb.py')
        
    ,Detail('system_type',         'varchar',['storage_directory'],    systemType,      'ksb','Looks at Atoms, decides bulk/surface/molecule from fill fraction (cutoffs at 15%% and 0.1%%)','/home/ksb/scripts/datalog/details_ksb.py')
    
    ,Detail('crystal_system',      'varchar',['system_type_ksb'
                                             ,'storage_directory'],    getCrystalSystem,'ksb','Looks at Atoms, decides what crystal lattice best describes it','/home/ksb/scripts/datalog/details_ksb.py')
    
    ,Detail('keld_lattice_A',      'numeric',['job_name'],             getExptA,        'ksb','Checks if name is in keld database, returns lattice parameter in A','/home/ksb/scripts/datalog/details_ksb.py')
    
    ,Detail('keld_lattice_BM',     'numeric',['job_name'],             getExptBM,       'ksb','Checks if name is in keld database, returns bulk modulus in GPa','/home/ksb/scripts/datalog/details_ksb.py')

    ,Detail('lattice_parameter',   'numeric',['system_type_ksb'
                                             ,'storage_directory'],    getLatticeA,     'ksb','Looks at final image and calculates lattice parameter (if bulk)','/home/ksb/scripts/datalog/details_ksb.py')

   ,Detail('pointgroup',          'varchar', ['system_type_ksb'
                                            ,'storage_directory'],   getPointgroup,       'ksb','Gets point group if molecule','/home/ksb/scripts/datalog/details_ksb.py')

   ,Detail('symmetry_number',     'integer', ['pointgroup'],         pointGroupToSymnum,  'ksb','Maps point group to symmetry number (for vibrational analysis)','/home/ksb/scripts/datalog/details_ksb.py')

   ,Detail('molecule_shape',      'varchar', ['pointgroup'],         getMolShape,         'ksb','Maps point group to {monatomic,linear,nonlinear}','/home/ksb/scripts/datalog/details_ksb.py')

   ,Detail('spacegroup',          'integer', ['system_type_ksb'
                                            ,'storage_directory'],   getSpacegroup,        'ksb','Gets space group if bulk crystal','/home/ksb/scripts/datalog/details_ksb.py')

   ,Detail('numbers_str',         'varchar',['storage_directory'],  getNumbersStr,         'ksb','Prints string of Atoms.get_atomic_numbers()','/home/ksb/scripts/datalog/details_ksb.py')

   ,Detail('blank',               'varchar',['storage_directory'],  lambda x: '',         'ksb','Prints empty string','/home/ksb/scripts/datalog/details_ksb.py')

   ,Detail('surface_area',        'numeric',['storage_directory'
                                            ,'system_type_ksb'
                                            ,'symmetric_ksb'],        surfaceArea,          'ksb','Prints empty string','/home/ksb/scripts/datalog/details_ksb.py')

    ,Detail('eform',              'numeric',['raw_energy'
,'dftcode'
,'xc'
,'pw'
,'kptden_x'
,'kptden_y'
,'kptden_z'

                                            ,'jobcalc'
                                            ,'raw_energy'
                                            ,'numbers_str'],        getEform, 'ksb','Looks up reference','')

]

dDict = {d.colname:d for d in allDetails}
###########################################################################################################
