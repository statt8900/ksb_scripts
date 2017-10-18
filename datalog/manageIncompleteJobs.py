#External modules
import fireworks,os,subprocess,datetime,manageSharedDatabase
#Internal Modules
import jobs,misc

lpad = fireworks.LaunchPad(host='suncatls2.slac.stanford.edu',name='krisbrown',username='krisbrown',password='krisbrown')

def fizzleLostRuns(expiration_secs=3600*10, max_runtime=None, min_runtime=None):
    lpad.detect_lostruns(expiration_secs, fizzle=True, max_runtime=max_runtime, min_runtime=min_runtime)
##################
# Helper functions
#-----------------
def filterError(keywords): return [i for i in fizzledFWIDS() if getErr(i) is not None and any([k in getErr(i) for k in keywords])]

def delete(fwid): lpad.archive_wf(fwid)

##################
# Lists of fwids
#-----------------
def allFWIDS():         return lpad.get_fw_ids()
def fizzledFWIDS():     return lpad.get_fw_ids({'state':'FIZZLED'})
def completeFWIDS():    return lpad.get_fw_ids({'state':'COMPLETED'})
def archivedFWIDS():    return lpad.get_fw_ids({'state':'ARCHIVED'})
def readyFWIDS():       return lpad.get_fw_ids({'state':'READY'})
def reservedFWIDS():    return lpad.get_fw_ids({'state':'RESERVED'})
def runningFWIDS():     return lpad.get_fw_ids({'state':'RUNNING'})
def incompleteFWIDS():  return fizzledFWIDS()+runningFWIDS()+reservedFWIDS()+readyFWIDS() #include archived? don't want to repeat unconverged. but if it's archived b/c of a mistake...
def runningTimeout():   return [i for i in runningFWIDS() if getTimeout(i)]

def jobsSinceXhours(x): return [i for i in allFWIDS() if timeSince(i) < x and timeSince(i) is not None]

##############
# Known Errors
#-------------

errorDict = {'timeout':     ['TIME','time','User defined signal 2']
            ,'kohnsham':    ['KohnShamConvergenceError']
            ,'stres_vdw':   ['ValueError: Extra data:']
            ,'scf':         ['Error in routine stres_vdW_DF']
            ,'diagonalize': ['RuntimeError: SCF calculation failed']}

def sortErrLog(err):
    for e in errorDict.keys():
        if any([x in err for x in errorDict[e]]): return e
    return 'unknown'

def errFWIDS(err): return filterError(errorDict[err])

def errReport():
    from prettytable import PrettyTable
    fizz = fizzledFWIDS()

    lFizz = len(fizz)
    lReady,lRes,lRun,lArch,lComp,lTot,lrTO = map(len,[readyFWIDS(),reservedFWIDS(),runningFWIDS(),archivedFWIDS(),completeFWIDS(),lpad.get_fw_ids(),runningTimeout()])

    print "\nStatus of all fireworks..."
    x = PrettyTable(['Ready','Reserved','Running','Running (timed out)','Fizzled','Archived','Completed','Total'])
    x.add_row([lReady,lRes,lRun-lrTO,lrTO,lFizz,lArch,lComp,lTot])
    print x
    print '\nDiagnosis of fizzled fireworks'
    
    unknowns,countDict = [],{x:0 for x in errorDict.keys()+['unknown']}
    for i in fizz:
        n = countDict['unknown']
        countDict[sortErrLog(getErr(i))]+=1
        if countDict['unknown']>n: unknowns.append(i)

    cD = {k:v for k,v in countDict.items() if v != 0}
    y = PrettyTable(cD.keys())
    y.add_row(cD.values())
    print y
    if len(unknowns) > 0: print "unknown FWIDS ",' '.join(map(str,unknowns))

##########################
# Things to get from a job
#-------------------------
def timeSince(fwid):
    try:
        fwdict      = lpad.get_fw_dict_by_id(fwid)
        states      = fwdict['launches'][-1]['state_history']

        starttime   = [state['created_on'] for state in states if state['state']=='RUNNING'][0]

        format      = '%Y-%m-%dT%H:%M:%S.%f'
        start       = datetime.datetime.strptime(starttime,format)
        now         = datetime.datetime.utcnow()
        delta       = (now - start).total_seconds()/3600.

        return delta
    except TypeError: return None
    except KeyError: return None
    except IndexError: return None

def getTimeout(fwid):
    try:
        fwdict      = lpad.get_fw_dict_by_id(fwid)
        starttime   = fwdict['launches'][-1]['action']['time_start']
    
        format      = '%Y-%m-%dT%H:%M:%S.%f'
        start       = datetime.datetime.strptime(starttime,format)
        now         = datetime.datetime.utcnow()

        delta       = (now - start).total_seconds()

        walltime    = fwdict['spec']['_queueadapter']['walltime']
        walltimelist= walltime.split(':')
        wallseconds = sum([int(x)*y for x,y in zip(walltimelist,[3600,60,1])])  # HH:MM or HH:MM:SS
        if delta > wallseconds: return True
        else: return False
    except TypeError: return False
def getErr(fwid):
    import glob 
    ld =  lpad.get_fw_dict_by_id(fwid)['spec']['_launch_dir']
    if 'scratch' in ld:
        err = glob.glob(ld+'/*.error') 
        with open(err[0],'r') as f: return f.read()
    elif 'nfs' in ld:
            d   = subprocess.Popen(['ssh','ksb@suncatls1.slac.stanford.edu', 'cat %s/*.error'%ld], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            dout, err = d.communicate()
            return dout
    return None

def fwid2strjob(fwid):  return str(jobs.Job(lpad.get_fw_dict_by_id(fwid)['spec']['params']))
def fwid2params(fwid):  return lpad.get_fw_dict_by_id(fwid)['spec']['params']

#########################
# More complicated things
#------------------------
def param(key,value): return lambda x: fwid2params(x)[key]==value

def relaunch(predicate=lambda x:True):
    """Relaunch timed out jobs and unconverged GPAW jobs"""
    fizzleLostRuns();lpad.detect_unreserved(expiration_secs=3600*24*7, rerun=True)
    timeouts,unconvergeds = filter(predicate,errFWIDS('timeout')),filter(predicate,errFWIDS('kohnsham'))
    tQuestion   = "Do you want to relaunch %d timed out runs?"%len(timeouts)
    uQuestion   = "Do you want to relaunch %d unconverged jobs?"%len(unconvergeds)
    if misc.ask(tQuestion):
        for fwid in timeouts:
            q               = lpad.get_fw_dict_by_id(fwid)['spec']['_queueadapter']
            q['walltime']   =  misc.doubleTime(q['walltime'])
            lpad.update_spec([fwid],{'_queueadapter':q})
            lpad.rerun_fw(fwid)
            if q['walltime'][:2] == '40': print 'warning, 40h job with fwid ',fwid
    if misc.ask(uQuestion):
        inc = listOfIncompleteJobStrs()

        for fwid in unconvergeds:
            p = fwid2params(fwid)
            p['sigma']  += 0.05
            p['mixing'] *= 0.5
            job         = jobs.Job(p)
            lpad.archive_wf(fwid)
            if job.new(inc): job.submit()

    misc.launch()

def listOfIncompleteJobStrs(): 
    """Use this for checking whether a proposed job is a duplicate with an unfinished job"""
    manageSharedDatabase.updateDB('strjob',deps=True)
    return map(fwid2strjob,incompleteFWIDS())

def printErr():
    for e in map(getErr,fizzledFWIDS()): print e


