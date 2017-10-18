# ---------------------------------------------------------------------------------------
# Things absolutely necessary to log during job (Everything else done in post-processing)
# ---------------------------------------------------------------------------------------
# user
# timestamp
# working directory
# pickled version of ASE atoms object
# copies of output/error/log files and job script

# -------------------------------------
# Optional things declarable at runtime
# -------------------------------------
# job_name, job_comments, dftcode, AND MORE?
# WE WANT THIS TO BE COMPLETE, BUT COULD BE EXPANDED IN THE FUTURE

def safeMkdir(newDir):
    import os
    try:            os.makedirs(newDir)
    except OSError: pass # may fail due to multiple threads or already existing 

def safeCopy(path,name):
    import shutil
    try:            shutil.copyfile(path,name)
    except IOError: pass

def clusterRoot():
    import os
    hostname = os.environ['HOSTNAME'].lower()
    if      'sh'    in hostname: return '/scratch/users/ksb/share/'
    elif   'gpu-15' in hostname: return '/scratch/users/ksb/share/'
    elif    'su'    in hostname: return '/nfs/slac/g/suncatfs/ksb/share/' #important to distinguish suncat2 and 3?
    else: raise ValueError, "clusterRoot did not detect SH or SU in %s"%hostname

def log(optimized_atoms, job_name='',job_comments=''):             
    """  
    Execute at the end of all scripts. Copies info to external storage, where it should be added to a database
    """
    import getpass,time,os,sys,pickle,sqlite3,json
        
    ######################################################
    print "Collecting environment data..."
    #-----------------------------------------------------
    public_storage    = clusterRoot()+'jobs/'
    
    user              = getpass.getuser()
    timestamp         = time.time()
    working_directory = os.getcwd()
    script            = sys.argv[0].split('/')[-1] # e.g. opt.py

    ######################################################
    print "Agglomerating data to be logged in database..."
    #-----------------------------------------------------
    newDir      = public_storage+'%s/%s/'%(user,str(timestamp).replace('.','')) #name of new directory
    raw_pckl,logpth,pwinp,qnlogpth,errpth,outpth,scriptpth = [newDir+x for x in ['raw.pckl','log','pw.inp','qn.log','myjob.err','myjob.out','script.py']]
    
    allCols     = ['deleted','user',    'timestamp', 'working_directory','storage_directory','job_name','job_comments'] 
    binds       = [0,         user,  int(timestamp),  working_directory,  newDir,             job_name,  job_comments ]       
    runtime     = dict(zip(allCols,binds))
    
    #################################################
    print "Creating and Populating copy directory..."
    #------------------------------------------------
    safeMkdir(newDir)
    
    with open(newDir+'raw.pckl','w') as f: f.write(pickle.dumps(optimized_atoms)) # Write final Atoms object to file
    with open(newDir+'runtime.json','w') as f: f.write(json.dumps(runtime))       # Write runtime info to file

    for root, dirs, files in os.walk(working_directory, topdown=False):
        for name in files:
            originalPath = os.path.join(root, name)
            if name in ['log','pw.inp','qn.log'] : safeCopy(originalPath,newDir+name) 
            elif name == script: safeCopy(originalPath,newDir+'script.py')             
            elif '.err' in name: safeCopy(originalPath,newDir+'myjob.err')
            elif '.out' in name: safeCopy(originalPath,newDir+'myjob.out')
            elif os.stat(originalPath).st_size < 1e9: safeCopy(originalPath,newDir+name) 

    return 0

#-------------------------------------------------------
# 'Derivable' observables: 
# ------------------------
# Everything above concerned an easy/simple way of storing all 'irreducible' data. 
# Now we need to transform it in a way that's convenient for anyone to query/visualize the data. 
# Ideas for that in /home/ksb/scripts/datalog/datatransform.txt


