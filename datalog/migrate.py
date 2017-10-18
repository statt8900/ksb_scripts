import os,getpass,time,os,pickle,sqlite3,shutil,json,dbase


deleted     = 0  # False
dftcode     = 'gpaw'
user        = 'ksb'
jobsroot    = '/scratch/users/ksb/fireworks/jobs/'
public_database   = '/scratch/users/ksb/share/suncatdata.db'
public_storage    = '/scratch/users/ksb/share/jobs/ksb/'

os.chdir(jobsroot)

def checkRepeat(direct):
    connection  = sqlite3.connect('/scratch/users/ksb/share/suncatdata.db')
    cursor      = connection.cursor()
    cursor.execute("select * from completed where working_directory = '%s'"%direct,[])
    output_array = cursor.fetchall()
    connection.commit(); cursor.close(); connection.close()
    return len(output_array) == 0 


for d in os.listdir(jobsroot):
    working_directory   = jobsroot+d
    if os.path.exists(d+'/result.json') and checkRepeat(working_directory):
        if os.stat(d+'/result.json').st_size == 0 :     
            print 'removing directory ',working_directory
            shutil.rmtree(working_directory)
        else:
            print 'entering ',d,' directory = ',working_directory
            with open(d+'/result.json','r') as f: result = json.loads(f.read())
            if 'name' not in result.keys(): 

            job_name    = result['name']    
            job_comments,script= '',''
            timestamp   = time.time()
            
            newDir  = public_storage+'%s/'%(str(timestamp).replace('.',''))
            os.makedirs(newDir)
            with open(newDir+'raw.pckl','w') as f: f.write(result['finaltraj_pckl'])
               
            def safeCopy(path,name):
                try:            shutil.copyfile(path,newDir+name)
                except IOError: pass

            raw_pckl,logpth,pwinp,qnlogpth,errpth,outpth,scriptpth = [newDir+x for x in ['raw.pckl','log','pw.inp','qn.log','myjob.err','myjob.out','script.py']]
            
            for root, dirs, files in os.walk(working_directory, topdown=False):
                for name in files:
                    originalPath = os.path.join(root, name)
                    if name == 'qn.log': safeCopy(originalPath,name) 
                    elif name == script: safeCopy(originalPath,'script.py')             
                    elif '.err' in name: safeCopy(originalPath,'myjob.err')
                    elif '.out' in name: safeCopy(originalPath,'myjob.out')
                    elif '.txt' in name: safeCopy(originalPath,'log') 


            allCols     = ['deleted','user',    'timestamp', 'working_directory','raw_pckl','log_path','qnlog_path','pwinp_path','err_path','out_path','script_path','dftcode','job_name','job_comments'] 
            binds       = [0,        'ksb',  int(timestamp),  working_directory,  raw_pckl,  logpth,    qnlogpth,    pwinp,       errpth,    outpth,    scriptpth,    dftcode,  job_name,  job_comments ]       
            listOfColNames,questionMarks   = ','.join(allCols)  , '?'+',?'*(len(allCols)-1) 
            command     = 'INSERT into completed (%s) VALUES (%s)'%(listOfColNames,questionMarks)
            connection  = sqlite3.connect(public_database,timeout=60)
            cursor      = connection.cursor()
            cursor.execute(command, binds); connection.commit(); cursor.close(); connection.close()
