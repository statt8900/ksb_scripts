#!/nfs/slac/g/suncatfs/sw/py2.7.13/bin/pythonwrapper
        # ! /usr/bin/env python276 <- joe's version had this
import sys
import os
import time
from subprocess import check_output
from glob import glob
sys.path.insert(0,'/afs/slac.stanford.edu/u/if/jgauth32/PythonModules/')
from prettytable import PrettyTable


def glob_sort(list):
    nums=[]
    for item in list:
        string = item[-7:-4]
        nums.append(int(float(string)))
    nums.sort()
    return nums[-1]

def analyze_LSFout():
    cmd = ['bjobs','-l']
    stdout = check_output('bjobs -l',shell=True)
	#stdout = process.stdout.read()

    monthdict = {'Jan':1,'Feb':2,'Mar':3,'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug':8,'Sep':9,'Oct':10,'Nov':11,'Dec':12}
    renames = {'Job':'Number','Job Priority':'Priority','Submitted from host':'Host','CWD':'Path','Output File':'StdOut','Error File':'StdErr','Execution CWD':'ExecutionPath'}

    jobdict = {}
    stdout = stdout.replace('\t','')
    stdout = stdout.replace('\n'+' '*21,'')
    jobs = stdout.split('------')

    jobs = [j for j in jobs if j.replace('-','').strip()]

    job_order = []
    for j in jobs:
        "NEW JOB"
        vardict = {}
        lines = [l for l in j.split('\n') if l.strip()]
        line1 = lines.pop(0)
        defs = line1.split(',')

        line2 = lines.pop(0)
        submittime = ':'.join(line2.split(':')[0:3])
        year = time.strftime('%Y',time.localtime())
        submittime_dt = time.strptime(submittime + " " + year,'%a %b %d %H:%M:%S %Y')
        line2 = line2.replace(submittime+': ','')
        defs += line2.split(',')
        
        line3 = lines.pop(0)
        if 'Started' not in line3:
            lines.pop(0)
            line3 = lines.pop(0)

        if 'Started' in line3:
            starttime = ':'.join(line3.split(':')[0:3])
            starttime_dt = time.strptime(starttime + " " + year,'%a %b %d %H:%M:%S %Y')
            day,month,date,tm = starttime.split(None,3)
            month = monthdict[month]
            hr,min,sec = tm.split(':')
            now = time.localtime()
            then = list(time.localtime())
            then[1] = int(month)
            then[2] = int(date)
            then[3] = int(hr)
            then[4] = int(min)
            then = tuple(then)
            diff = time.mktime(now) - time.mktime(then)
            runtime = round(diff/3600.,2)
            vardict['Runtime'] = runtime
            vardict['Qtime'] = round((time.mktime(starttime_dt) - time.mktime(submittime_dt))/3600.,1)
            line3 = line3.replace(starttime+': ','')
            defs += line3.split(',')
        else:
            vardict['Runtime'] = '-'
            vardict['Qtime'] = round((time.mktime(time.localtime()) - time.mktime(submittime_dt))/3600.,1)
        
        for d in defs:
            if 'Processors Requested' in d:
                num,rubbish = d.split(' ',1)
                vardict['Cores'] = num
                vardict['NodeIDs'] = '-'
            elif 'Started on' in d:
                rubbish,val = d.split('<',1)
                val = val.replace('< >',',')
                val = val.replace('>','')
                if '<' in val:
                    val = val.split('<')
                    node_ids = []
                    for i in val:
                        node_ids.append(i.rstrip()[-4:])
                    val = ','.join(node_ids)
                else:
                    val = val[-4:]
                vardict['NodeIDs'] = val

            else:
                key,val = d.rsplit(' ',1)
                key = key.strip()
                val = val.replace('<','')
                val = val.replace('>','')
                if key in renames:
                    key = renames[key]
                vardict[key] = val
        
        #print vardict
        try:
            logfiles = glob(vardict['Path']+'/*qn*.log')
            lognum=glob_sort(logfiles)
            logfile=vardict['Path']+'/qn%05i.log' % lognum
            if logfile:
                log = open(logfile).readlines()
                for i,line in enumerate(log):
                    if ' 0' in line[5:10]:
                        start = i
                    if 'BFGS' in line:
                        end = i
                """
                lines = log.readlines()
                n_iter = len(lines)
                log.close()
                last = lines[-1]
                fmax = round(float(last.split()[-1]),2)
                """
                vardict['fmax'] = str(float(log[i][-7:-1]))
                #print vardict['fmax']
                vardict['N_iter'] = str(float(len(log)-start-1))
                #if 'Runtime' in vardict:
                if vardict['Runtime'] != '-':
                    vardict['AvgIter'] = str(round(float(vardict['Runtime'])/float(len(log)-start-1),2))
                else:
                    vardict['AvgIter'] = '-'
            else:
                vardict['fmax'] = '-'
                vardict['AvgIter'] = '-'
                vardict['N_iter'] = '-'
        except:
            vardict['fmax'] = '-'
            vardict['AvgIter'] = '-'
            vardict['N_iter'] = '-'
    
        jobdict[vardict['Number']] = vardict
        job_order.append(vardict['Number'])

    return jobdict,job_order


filter = {}
showmax = 75
if len(sys.argv) > 1:
    showmax = int(sys.argv[1])
    if len(sys.argv) > 2:
        filter = {'User':sys.argv[2]}

job_dict,job_order = globals()['analyze_LSFout']()

##Check prev_jobs.txt and compare to current jobs to determine which jobs are new and which have expired.
current_ids = job_order
current_jobs = check_output('bjobs',shell=True).split('\n')[0:]
prev_jobs_path = '/u/if/ksb/usr/bin/prev_jobs.txt'
exp_jobs = []; new_jobs = []
f = open(prev_jobs_path,'r')
prev_jobs = f.read().splitlines()[3:-1]
prev_ids = [prev_job[2:8].strip() for prev_job in prev_jobs]
f.close()
prev_time = time.ctime(os.path.getmtime(prev_jobs_path))
"""
for current_job,current_id in zip(current_jobs[1:],current_ids):
    #print current_job
    if current_id not in prev_ids:
        new_jobs.append(current_job)
"""
for prev_job,prev_id in zip(prev_jobs,prev_ids):
    if prev_id not in current_ids:
        exp_jobs.append([prev_job for prev_job in prev_jobs if prev_id in prev_job][0])

print "Last checked on %s" %prev_time
if exp_jobs != []:
    print 'Jobs completed or expired since last check:'
    for exp_job in exp_jobs: print exp_job
    print '\n'
"""
if new_jobs != []:
    print 'Jobs submitted since last check:\n'
    for new_job in new_jobs: print new_job
"""
##Print Table
"""
if len(job_order) > showmax:
    full = len(job_order)
    job_order = job_order[:showmax]
    print 'Showing '+str(showmax)+' of '+str(full)+' jobs.'
"""
show_fields = ['Number','NodeIDs','Queue','Qtime','Job Name','Runtime','Path']

"""
labels = [t for t in show_fields]
#labels[1] += '\t' #align runtime properly
print '\t'.join(labels)
"""
table = PrettyTable(show_fields)
table.align['NodeIDs'] = 'l'
table.align['Path'] = 'l'

for job in job_order:
    info = []
    for field in show_fields:
        try:
            info.append(job_dict[job][field])
        except KeyError:
            info.append('')
    info = [str(i) for i in info]
    show = True
    for key in filter:
        if job_dict[job][key] != filter[key]:
            show = False
    if show == True:
        #print '\t'.join(info)
        table.add_row(info)

print table
with open(prev_jobs_path,'w') as file:
    file.write(table.get_string())
