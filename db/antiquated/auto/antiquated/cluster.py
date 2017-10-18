import os
from objectClass import Object
from printParse  import multiLine,printTime
from standardScripts import getCluster

class Cluster(Object):
	def __init__(self,name,header,runFile,root,result,python,timeout,unconverged,fail,errorFile):
		self.name        = name     # String
		self.header      = header   # Time     -> String
		self.runFile     = runFile  # FilePath -> String
		self.root        = root     # FilePath
		self.result      = result   # FilePath
		self.python      = python   # String
		self.timeout     = timeout  # String -> Bool
		self.unconverged = unconverged # STring -> Bool
		self.fail        = fail
		self.errorFile   = errorFile

	def runFileHeader(self,job): 
		return (multiLine(	['#!/bin/bash'
						,sherlock.header(job.timeLim)+'\n'
						,"#Slurm parameters"
						,"NTASKS=`echo $SLURM_TASKS_PER_NODE|tr '(' ' '|awk '{print $1}'`"
						,"NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`"
						,'NCPU=`echo " $NTASKS * $NNODES " | bc`'
						,'echo "$NCPU"'
						,"#load gpaw-specific paths"
						,"source /scratch/users/ksb/gpaw/paths.bash"
						,'echo "$1"'
						,"#run parallel gpaw"
						,"mpirun -n $NCPU gpaw-python sherlock_opt.py\n"])
				 if job.dftCode == 'gpaw' else self.runFile)

sherlock = Cluster ('sherlock',lambda t: multiLine(	['#SBATCH -p iric'
													,'#SBATCH -x gpu-14-1,sh-20-35'
													,'#SBATCH --job-name=myjob'
													,'#SBATCH --output=myjob.out'
													,'#SBATCH --error=myjob.err'
													,'#SBATCH --time={0}:00'.format(printTime(t))
													,'#SBATCH --qos=iric'
													,'#SBATCH --nodes=1'
													,'#SBATCH --mem-per-cpu=4000'
													,'#SBATCH --mail-type=END,FAIL'
													,'#SBATCH  --mail-user=ksb@stanford.edu'
													,'#SBATCH --ntasks-per-node=16\n'])
						, "sbatch sherlock_opt.py\n"
						,"/scratch/users/ksb/auto"
						,"/scratch/users/ksb/autoResult"
						,'#!/scratch/PI/suncat/sw/bin/python'
						,lambda x:'TIMEOUT' in x
						,lambda x:'Did not converge' in x
						,lambda x: 'IndexError' in x
						,'myjob.err')

suncat2 = Cluster ('suncat2',lambda t: multiLine(	['#!/usr/bin/env python'
													,'#LSF -q suncat2'
													,'#LSF -n 12'
													,'#LSF -W {0}'.format(printTime(t))
													,'#LSF -o opt.log'
													,'#LSF -e err.log'
													,'#LSF -sp 90'
													,'#LSF -N\n'])
						,"esp-ver-bsub 20 python suncat2_opt.py\n"
						,"/a/suncatfs1/u1/ksb/auto"
						,"/a/suncatfs1/u1/ksb/autoResult"
						,"#!/usr/bin/env python"
						,lambda x:'TIMEOUT' in x
						,lambda x:'Did not converge' in x
						,lambda x: 'IndexError' in x
						,'err.log')




clusterDict = {'sherlock':sherlock,'suncat2':suncat2}

def getCluster(n):
	if 'sh' in os.environ['HOSTNAME'].lower(): return sherlock
	elif n == 1 : return suncat1 #important to distinguish suncat2 and 3?
	elif n == 2 : return suncat2
	elif n == 3 : return suncat3
	else: raise ValueError, "getCluster needs 1,2, or 3 to choose suncat partition"
