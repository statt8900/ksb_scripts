#!/home/vossj/suncat/bin/python
from subprocess import check_output
import argparse
import os.path
import time
import sys
from ase.io import read
from ase import Atoms
from glob import glob
sys.path.insert(0,'/home/ksb/PythonModules/')
from prettytable import PrettyTable


parser = argparse.ArgumentParser()
#parser.add_argument('-f','--fresh',action='store_true')
#parser.add_argument('-s','--short',action='store_true',help='Use this option for abbreviated output that will fit most small screen sizes')
parser.add_argument('-p','--path',default=None,help='Return path only, used by fj in .bashrc')
parser.add_argument('-n','--norecord',action='store_true',help='Use if you dont want job check to be recorded')
args = parser.parse_args()
new_jobs = []
exp_jobs = []
prev_jobs_path = '/home/ksb/usr/bin/prev_jobs.txt'
user = 'ksb'

if args.path != None:
    job_info = check_output("scontrol show jobid %s" %args.path,shell=True).split('\n')
    path = job_info[19][11:]
    if path[-4:]=='.err': #job hasn't started yet, different scontrol dialog
        path = '/'+path[:-9]
    print path
    sys.exit()


##Check prev_jobs.txt and compare to current jobs to determine which jobs are new and which have expired.

current_jobs = check_output(['squeue','-u','%s'%user,'-o','%.18i %.9P %.24j %.8u %.2t %.10M %.6D %.20V %.20S %Z']).split('\n')
current_jobs = filter(None,current_jobs)
current_ids = [str(int(current_job[:18])) for current_job in current_jobs[1:]]
try:
    f = open(prev_jobs_path,'r')
    prev_jobs = f.read().splitlines()[3:-1]
    #for prev_job in prev_jobs:
	#print str(int(prev_job.split('|')[1]))
    prev_ids = [str(int(prev_job.split('|')[1])) for prev_job in prev_jobs]
    f.close()
    prev_time = time.ctime(os.path.getmtime(prev_jobs_path))
    for current_job,current_id in zip(current_jobs[1:],current_ids):
        #print current_job
        if current_id not in prev_ids:
            new_jobs.append(current_job)

    for prev_job,prev_id in zip(prev_jobs,prev_ids):
        if prev_id not in current_ids:
            exp_jobs.append([prev_job for prev_job in prev_jobs if prev_id in prev_job][0])

    print "Last checked on %s" %prev_time

    if exp_jobs != []:
        print 'Jobs completed or expired since last check:\n'
        for exp_job in exp_jobs: print exp_job
    if new_jobs != []:
        print 'Jobs submitted since last check:\n'
        for new_job in new_jobs: print new_job


except:
    print "No prev_jobs.txt found. Writing prev_jobs.txt now"

##Adding extra arguments (fmax, etc)
jobs = current_jobs[1:]
ids = [str(int(job[:18])) for job in jobs]
table = PrettyTable(['Job ID','Partition','Name','ST','Run Time','#','Qtime','Path'])
table.align['Path'] = 'l'
table.align['Name'] = 'l'

for (job,id) in zip(jobs,ids):
    partition = job.split()[1].strip()
    if len(job.split()) == 10:
        name = job.split()[2].strip()
    else:
        name = ''

    num_nodes = job.split()[-4].strip()

    ST = job.split()[-6].strip()

    path = job.split()[-1].strip()
    if path[0:13] == '/home/%s'%user:
        path_disp = '~' + path[13:]
    else:
        path_disp = path

    runtime = job.split()[-5].strip()
    submittime = time.strptime(job.split()[-3].strip().replace('T',' '),'%Y-%m-%d %H:%M:%S')
    if runtime == "0:00":
        qtime = round((time.mktime(time.localtime()) - time.mktime(submittime))/3600.,2)
    else:
        starttime = time.strptime(job.split()[-2].strip().replace('T',' '),'%Y-%m-%d %H:%M:%S')
        qtime = round((time.mktime(starttime) - time.mktime(submittime))/3600.,2)

    days = 0
    hours = 0
    minutes = 0
    seconds = float(runtime.split(':')[-1])
    if len(runtime.split(':')) > 1:
        minutes = float(runtime.split(':')[-2])
    if len(runtime.split(':')) > 2:
        if '-' in runtime:
            hours = float(runtime.split(':')[-3].split('-')[-1])
            days = float(runtime.split('-')[0])
        else:
            hours = float(runtime.split(':')[-3])
    runtime_h = round(days*24 + hours + minutes/60. + seconds/3600.,2)

    try:
        qnlog = open(path + '/qn.log','r').readlines()
        start = 0
        for i,line in enumerate(qnlog):
            if 'kpts' in line:
                start = i #most recent starting point
            if 'BFGS' in line:
                end = i #most recent BFGS step
        iter_num = float(len(qnlog)-start-1)
        fmax = round(float(qnlog[i][57:63]),2)
        iter_avg = round(runtime_h/iter_num,2)
        #print '%s\t%3.2f\t%4.3f\t%s' %(iter_avg,fmax,path)
    except:
        if os.path.isfile(path+'/OUTCAR'):
            iter_avg = '-'
            try:
                atoms = read(path+'/OUTCAR')
                unconstrained = [atom.index for atom in atoms if atom.index not in list(atoms.constraints[0].index)]
                ftemp = []
                for ind in unconstrained:
                    f = atoms.get_forces()[ind]
                    ftemp.append((f[0]**2+f[1]**2+f[2]**2)**0.5)
                ftemp.sort()
                fmax = round(ftemp[-1],2)
            except:
                fmax = '-'
        else:
            fmax = '-'
            iter_avg = '-'
    table.add_row([id,partition,name,ST,runtime_h,num_nodes,qtime,path_disp])

print table.get_string(sortby="Path")
os.system('mv /home/%s/usr/bin/prev_jobs.txt /home/%s/usr/bin/prev_jobs.txt.bak'%(user,user))
with open('/home/%s/usr/bin/prev_jobs.txt'%user,'w') as file:
    file.write(table.get_string(sortby="Path"))
