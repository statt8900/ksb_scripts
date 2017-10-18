#External Modules
import ast,json,subprocess,os,sys,shutil,itertools,signal,fireworks,pickle
#Internal Modules
import databaseFuncs as db
import detailClass,misc
from misc import readOnSherOrSlac,parseLine
from details_ksb import getEng,getPSPpath,getDFTcode

#####################
# Auxillary Functions
#--------------------
lpad = fireworks.LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')
def ask(x): return raw_input(x+'\n(y/n)---> ').lower() in ['y','yes']

def abbreviateDict(d): return '\n'.join([str(k)+':'+(str(v) if len(str(v))<100 else '...') for k,v in d.items()])

######################
# Database maintanence
# --------------------
def delete(constraint,check=True):
    """
    Delete completed jobs that meet some constraint. Removes storage directory (but not working directory). 
    """
    output = db.query(['jobid','fwid','storage_directory'],constraint,order='jobid')
    for jid,fwid,path in output: 
        lpad.archive_wf(fwid)                                       # archive firework
        db.updateDB('deleted','jobid',jid,1,tableName='completed')  # note deletion in deleted column
        if not check or ask('Do you want to delete %s?'%path):      # delete storage directory  
            if 'scratch' in path: shutil.rmtree(path)
            elif 'nfs' in path: 
                d   = subprocess.Popen(['ssh','ksb@suncatls1.slac.stanford.edu', 'rm -r %s'%path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                dout, err = d.communicate()
            else: raise NotImplementedError
            print 'deleted!'

def duplicates(deleteFlag=False,cnst='1'):
    """
    Search database for duplicate strjobs. 
    """
    output = db.query(['fwid','strjob'],cnst)       # list of (fwid,strjob) pairs
    rptDict={}                                      # dictionary for repeat values (key = first fwid, value = list of duplicates)
    for fwid,strjob in output:
        for f,s in output:                          # double-FOR loop
            if f is None: print 'NONE FWID??? ',f,s # hopefully this isn't flagged
            if strjob == s and f!=fwid:             # condition for duplicate job
                if fwid not in list(itertools.chain.from_iterable(rptDict.values())): 
                    if fwid in rptDict.keys(): rptDict[fwid].append(f)  # add to the list
                    else: rptDict[fwid] = [f]                           # create a key,value pair
    print 'FWIDs with equal strjob entries: \n',abbreviateDict(rptDict) # summarize results
    if deleteFlag:
        delfws = list(itertools.chain.from_iterable(rptDict.values()))
        if ask('Are you sure you want to delete %d duplicates?'%len(delfws)):
            for f in delfws: delete('fwid = %d'%f,False)

def updateDB(cols = None, constraint = '1'
            , retry_completed = False, retry_failed = False, retry_undefined = False
            ,deps = True, verbose = False, new=False):
    """
    Select certain details to update (a space-separated string)
    Optionally give constraint over which rows to update
    Specify if computations should be done when they have previously been completed/failed/marked undefined.
    """
    if new: constraint+=' and planewave_cutoff is null' # SQL constraint that very likely picks out newly-added jobs
    def updateRow(rowid,d,retry_completed,retry_failed,retry_undefined,verbose):
        """
        May update a particular row/column pair, depending on the current state of the row/column pair (completed/failed/undefined/uninitialized)
        """
        if retry_completed or db.query1(d.colname,'jobid',rowid) is None:   # only try to compute if cell is empty OR we want to redo completed cells
            undef = db.query1(d.colname + '_undefined', 'jobid', rowid)
            if retry_undefined or undef is not 1:                           # also if we've previously determined the cell is undefined, only proceed if we want to redo those
                err = db.query1(d.colname + '_failed', 'jobid', rowid)
                if retry_failed or err is None or len(err) == 0:            # also if we've previously tried to calculate and the calculation failed, proceed only if we want to redo those
                    d.apply(rowid,verbose=False)                            # calculate
    
    loadJobs()                                                                             # load recently completed jobs
    if cols is None: ds = detailClass.allDetails                                           # No cols specified = all columns
    elif not deps: ds = filter(lambda x: x.colname in cols.split(),detailClass.allDetails) # only update column names explicitly specified, otherwise include dependencies (see below)
    else: ds = filter(lambda x: x.colname in [q.colname for q in detailClass.allDetails],list(set([item for sublist in [[detailClass.dDict[z] for z in d.inputcols if z in detailClass.dDict.keys()]+[d] for d in detailClass.allDetails if d.colname in cols.split()] for item in sublist]))) #sorry...


    jids = db.queryCol('jobid',constraint,order='jobid')        # list of row id's which we want to update

    for i,d in enumerate(detailClass.sortDetails(ds)):          # now, for every column (ordered such that dependencies come first)
        print 'Updating %s (%d/%d)...'%(d.colname,i+1,len(ds))
        d.addCol()                                              # add column to schema IF we've never added it before
        for j,jobid in enumerate(jids):
            print '%d/%d'%(j+1,len(jids)) ; sys.stdout.write("\033[F") # Cursor up one line
            updateRow(jobid,d,retry_completed,retry_failed,retry_undefined,verbose)


#############################
# Database startup/generation
# ---------------------------

def resetAll():
    #Repopulate completed table using the data stored in .../ksb/shared/jobs/...
    os.remove('/scratch/users/ksb/share/suncatdata.db')
    db.createSuncatDataDB()
    loadJobs()

def loadJobs():
    print "Loading recently completed jobs ..."
    sharedjobsSuncat,sharedjobsSherlock = '/nfs/slac/g/suncatfs/ksb/share/jobs','/scratch/users/ksb/share/jobs'
    suncout, err = subprocess.Popen(['ssh','ksb@suncatls1.slac.stanford.edu', 'cd %s;ls -d $PWD/*/*'%sharedjobsSuncat], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    sherout      = os.popen('cd %s; ls -d $PWD/*/*'%sharedjobsSherlock).read().split('\n')
    allJobs = [x+'/' for x in sherout + suncout.split('\n') if len(x)>0]
    newJobs = list(set(allJobs)-set(db.queryCol('storage_directory'))) # identify jobs not present already in database

    if len(allJobs) == len(newJobs):
        print allJobs
        print newJobs
        sys.exit()

    if len(newJobs) == 0:
        print 'No new jobs to add' ; return None
    else:
        print 'adding %d new jobs'%len(newJobs)
        for path in newJobs:  loadNewFolder(path)                                        # Extract calc / update calc table, update kwargs/completed/refeng


def loadNewFolder(storpath):
    """
    INPUT: a json of runtime data (generated by log()). OUTPUT: populates a new row of completed table
    """
    print 'loading ',storpath
    
    runtime              = json.loads(readOnSherOrSlac(storpath+'runtime.json')) # Load runtime information
    
    calcid               = getCalcID(storpath) 
    runtime['jobcalc']   = calcid                                               # Look at *calc* table and get ID (if no matching calc exists, add a row to table)
    
    kwargs               = json.loads(runtime.pop('kwargs'))                    # Take kwargs and move to seperate object

    listOfColNames,binds = ','.join(runtime.keys()), runtime.values()           # Runtime column names and corresponding values

    command = 'INSERT OR IGNORE into completed (%s) VALUES (%s) '%( listOfColNames, ','.join(['?']*len(binds)) ) # Make SQL command and
    db.sqlexecute(command,binds)                                                                                 # Insert it into *completed* 

    idnum = db.query1('jobid','storage_directory',runtime['storage_directory'],'completed')     # Get jobid of new job in *completed*
    if len(kwargs)>0: insertKWARGS(kwargs,idnum)                                                # Insert kwargs into *keyword* table
    print 'added kwargs, now jobid = %d and calcid = %d'%(idnum,calcid)
    addToRefEng(runtime,kwargs.get('job_name'),idnum,calcid)                                                # Adds row to *refeng* table if name in misc.refnames
    db.sqlexecute("insert or ignore into derived (derived_job) values (%d)"%idnum)              # Initialize row in *derived* table


def getCalcID(storpath):
    params   = getCalcParams(storpath)
    constrnt = ' AND '.join([k+'='+misc.s(v) for k,v in params.items()]).replace('=None',' is null')
    calcID   = db.queryCol('calc_id',constrnt,'calc',order='calc_id')
    if len(calcID)==1: return calcID[0]
    elif len(calcID)>1: raise ValueError, 'there should not be more than 1 calc ID matching input criteria'
    else: 
        command = 'INSERT into calc (%s) VALUES (%s) '%(','.join(params.keys()),','.join(['?']*len(params)))
        db.sqlexecute(command,params.values())  
        return db.queryCol('calc_id',constrnt,'calc',order='calc_id')[0]

def getCalcParams(strpth):
    dftcode     = getDFTcode(strpth)
    xc          = getXC(strpth,dftcode)
    pw          = getPW(strpth,dftcode)
    dw          = getDW(strpth,dftcode)
    kpts_json   = getKPTSJSON(strpth,dftcode)
    fmax        = getFMAX(strpth,dftcode)
    psp         = getPSP(strpth,dftcode)
    econv       = getECONV(strpth,dftcode)
    xtol        = getXTOL(strpth,dftcode)
    strain      = getSTRAIN(strpth,dftcode)
    dw          = getDW(strpth,dftcode)
    sigma       = getSIGMA(strpth,dftcode)
    nbands      = getNBANDS(strpth,dftcode)
    mixing      = getMIXING(strpth,dftcode)
    diag        = getDIAG(strpth,dftcode)
    spinpol     = getSPINPOL(strpth,dftcode)
    dipole      = getDIPOLE(strpth,dftcode)
    maxstep     = getMAXSTEP(strpth,dftcode)
    return {'dftcode':dftcode,'xc':xc,'pw':pw,'kpts_json':kpts_json,'fmax':fmax
            ,'psp':psp,'econv':econv,'xtol':xtol,'strain':strain,'dw':dw
            ,'sigma':sigma,'nbands':nbands,'mixing':mixing,'diag':diag
            ,'spinpol':spinpol,'dipole':dipole,'maxstep':maxstep}


def getXC(storagepath,dftcode):
    def getXC_GPAW(log): return parseLine(log,'xc:').split(':')[-1].strip()
    def getXC_QE(log):   return parseLine(log,'input_dft').split('=')[-1].replace(',','').replace("'",'').strip()
    if   dftcode == 'gpaw':            return getXC_GPAW(readOnSherOrSlac(storagepath+'log'))
    elif dftcode == 'quantumespresso': return getXC_QE(readOnSherOrSlac(storagepath+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'

def getPW(strpth,dftcode):
    def getPW_GPAW(log): return int(round(float((parseLine(log,'ecut').split()[-1].replace(',','')))))
    def getPW_QE(log):   return int(round(13.60569 * float(parseLine(log,'ecutwfc').split('=')[1][:7])))
    if   dftcode == 'gpaw':            return getPW_GPAW(readOnSherOrSlac(strpth+'log'))
    elif dftcode == 'quantumespresso': return getPW_QE(readOnSherOrSlac(strpth+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'

def getKPTSJSON(strpth,dftcode):
    def getKPTJSON_GPAW(log): return json.dumps(ast.literal_eval(','.join(parseLine(log,'kpts').split(': ')[-1].split())))
    def getKPTJSON_QE(log): return json.dumps([int(x) for x in log.split('\n')[-2].split()[:3]])
    if   dftcode == 'gpaw':            return getKPTJSON_GPAW(readOnSherOrSlac(strpth+'log'))
    elif dftcode == 'quantumespresso': return getKPTJSON_QE(readOnSherOrSlac(strpth+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'

def getFMAX(strpth,dftcode): 
    #print 'currently no way to ascertain fmax, assuming 0.05.  Try implementing something that scrapes all files to see if they'' json dictionaries'
    return json.loads(readOnSherOrSlac(strpth+'result.json')).get('fmax')

def getXTOL(strpth,dftcode):
    #print 'currently no way to ascertain xtol, assuming 0.01. Try implementing something that scrapes all files to see if they'' json dictionaries'
    return json.loads(readOnSherOrSlac(strpth+'result.json')).get('xtol')

def getSTRAIN(strpth,dftcode):
    #print 'currently no way to ascertain strain, assuming 0.04. Try implementing something that scrapes all files to see if they'' json dictionaries'
    return json.loads(readOnSherOrSlac(strpth+'result.json')).get('strain')


def getPSP(strpth,dftcode):
    pspth = getPSPpath(strpth,dftcode)
    if   'gbrv'   in pspth: return 'gbrv1.5pbe'
    elif 'dacapo' in pspth: return 'dacapo'
    elif 'sg15'   in pspth: return 'sg15'
    elif 'paw'    in pspth: return 'paw'
    else: raise NotImplementedError, 'New psp? path = ',pspth

def getECONV(strpth,dftcode):
    def getECONV_QE(log):    return round(13.60569 * float(parseLine(log,'conv_thr',-1).split('=')[1][:-1]),6)
    def getECONV_GPAW(log):  return float(parseLine(log,'convergence:').split('energy:')[-1].replace('}',''))
    if   dftcode == 'gpaw':            return getECONV_GPAW(readOnSherOrSlac(strpth+'log'))
    elif dftcode == 'quantumespresso': return getECONV_QE(readOnSherOrSlac(strpth+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'


def getDW(strpth,dftcode):
    def getDW_QE(log):   return round(13.60569 * float(parseLine(log,'ecutrho').split('=')[1][:7]))
    if   dftcode == 'gpaw':            return None
    elif dftcode == 'quantumespresso': return getDW_QE(readOnSherOrSlac(strpth+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'

def getSIGMA(strpth,dftcode):
    def getSIGMA_QE(log):    return round(13.60569 * float(parseLine(log,'degauss').split('=')[1][:7]),3)
    def getSIGMA_GPAW(log):  return float(parseLine(log,'Fermi-Dirac').split('=')[1].split()[0])
    if   dftcode == 'gpaw':            return getSIGMA_GPAW(readOnSherOrSlac(strpth+'log'))
    elif dftcode == 'quantumespresso': return getSIGMA_QE(readOnSherOrSlac(strpth+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'

def getNBANDS(strpth,dftcode): 
    if dftcode=='gpaw':
        return int(parseLine(readOnSherOrSlac(strpth+'log'),'nbands').split(':')[1])
    else:
        #print 'assuming nbands = -12'
        return json.loads(readOnSherOrSlac(strpth+'result.json')).get('nbands')

def getMIXING(strpth,dftcode):
    def getMIXING_QE(log):    return round(float(parseLine(log,'mixing_beta').split("=")[1][:-3]),3)
    def getMIXING_GPAW(log):  return float(parseLine(log,'Linear mixing parameter').split(':')[1])
    if   dftcode == 'gpaw':            return getMIXING_GPAW(readOnSherOrSlac(strpth+'log'))
    elif dftcode == 'quantumespresso': return getMIXING_QE(readOnSherOrSlac(strpth+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'

def getDIAG(strpth,dftcode):
    def getDIAG_QE(log):    return parseLine(log,'diagonalization').split("='")[1][:-2]
    def getDIAG_GPAW(log):  
        if 'Davidson' in log: return 'david'
        else: raise NotImplementedError, 'Different diagonalizer to check for in GPAW ? '+log
    if   dftcode == 'gpaw':            return getDIAG_GPAW(readOnSherOrSlac(strpth+'log'))
    elif dftcode == 'quantumespresso': return getDIAG_QE(readOnSherOrSlac(strpth+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'

def getSPINPOL(strpth,dftcode):
    def getSPINPOL_QE(log):    return int('starting_magnetization' in log)
    def getSPINPOL_GPAW(log):  return int('True' in parseLine(log,'spinpol'))
    if   dftcode == 'gpaw':            return getSPINPOL_GPAW(readOnSherOrSlac(strpth+'log'))
    elif dftcode == 'quantumespresso': return getSPINPOL_QE(readOnSherOrSlac(strpth+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'

def getDIPOLE(strpth,dftcode):
    def getDIPOLE_QE(log):    return int('dipfield=.true.' in log)
    if   dftcode == 'gpaw':            return 0
    elif dftcode == 'quantumespresso': return getDIPOLE_QE(readOnSherOrSlac(strpth+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'

def getMAXSTEP(strpth,dftcode):
    def getMAXSTEP_QE(log):    return int(parseLine(log,'electron_maxstep').split('=')[1][:-1])
    def getMAXSTEP_GPAW(log):  return int(parseLine(log,'maxiter').split(':')[-1])
    if   dftcode == 'gpaw':            return getMAXSTEP_GPAW(readOnSherOrSlac(strpth+'log'))
    elif dftcode == 'quantumespresso': return getMAXSTEP_QE(readOnSherOrSlac(strpth+'pw.inp'))
    else: raise NotImplementedError, 'unknown dftcode?'

def countElem(storpth,elemNum): # Counts how many of element X are in an atoms object
    return list(pickle.loads(readOnSherOrSlac(storpth+'raw.pckl')).get_atomic_numbers()).count(elemNum)

def addToRefEng(runtimedict,name,jobid,calcid):
    num = misc.nameToNum.get(name)
    stordir = runtimedict['storage_directory']
    if num is None: return None
    else:
        eng = getEng(stordir)
        count = countElem(stordir,num)
        if eng is None: return None
        else:
            command = 'INSERT OR IGNORE into refeng (element,refeng_job,refeng_calc,refeng) values (?,?,?,?)'
            db.sqlexecute(command,[num,jobid,calcid,eng/count])

def insertKWARGS(dic,idnum):
    """
    Take the keyword args dictionary in runtime.json and add elements to keyword table
    """
    for c in dic.keys():
        try: db.addCol(c,'varchar','keyword')  # If this is the first time this keyword has appeared, add the column
        except: pass                            
    cols,qMarks = 'keyword_job,'+','.join(dic.keys()),'?'+',?'*len(dic) # Column names in dictionary
    command = 'INSERT into keyword (%s) VALUES (%s) '%(cols,qMarks)     # SQL insert command
    binds = [idnum]+dic.values()                                        # Values in dictionary
    db.sqlexecute(command,binds)                                        

