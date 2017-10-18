import printParse 	as pp

""" 
The Cluster object and its only instances are declared

NEED TO ADD:
	- SHERLOCK2
	- CORI
	- EDISON
"""

class Cluster(object):
	def __init__(self,name,fworker,qfunc,launchdir):
		self.name 		= name
		self.fworker 	= fworker
		self.qfunc  	= qfunc
		self.launchdir  = launchdir
	def __repr__(self): return self.name

###########################################################################
###########################################################################

def sherlockQ(timeInHours):
	return {'_fw_name': 		'CommonAdapter'
			,'_fw_q_type': 		'SLURM'
			,'rocket_launch': 	'cd /scratch/users/ksb/fireworks;rlaunch singleshot' #'/scratch/PI/suncat/sw/lib/python2.7/site-packages/FireWorks-1.4.8-py2.7.egg/fireworks/scripts/rlaunch_run.py -w $SCRATCH/fireworks/my_fworker.yaml -l $SCRATCH/fireworks/my_launchpad.yaml singleshot'
			,'nodes': 			1
			,'ntasks_per_node': 16
			,'walltime': 		pp.printTime(timeInHours)+':00'
			,'queue': 			'owners,iric,normal'
			,'pre_rocket': 		'source /scratch/PI/suncat/sw/env.bash;source ~/scripts/rc/RCsher.sh' # 
			,'logdir': 			'/scratch/users/ksb/fireworks/logs/'}

def sherlockQGPAW(timeInHours):
	return {'_fw_name': 		'CommonAdapter'
			,'_fw_q_type': 		'SLURM'
			,'rocket_launch': 	'cd /scratch/users/ksb/fireworks;rlaunch singleshot'
			,'nodes': 			1
			,'ntasks_per_node': 16
			,'walltime': 		pp.printTime(timeInHours)+':00'
			,'queue': 			'owners,iric,normal'
			,'pre_rocket': 		';'.join(["source /scratch/PI/suncat/sw/env.bash"
										,"source ~/scripts/rc/RCsher.sh"
										,"export OMP_NUM_THREADS=1"
										,"export PYTHONPATH=/scratch/users/ksb/gpaw/ggafirst/install/lib/python2.7/site-packages:/scratch/users/ksb/gpaw/gpaw_sg15/lib/python2.7/site-packages:$PYTHONPATH"
										,"export PATH=/scratch/users/ksb/gpaw/ggafirst/install/bin:$PATH"
										,"export GPAW_SETUP_PATH=/scratch/users/ksb/gpaw/gpaw_sg15/norm_conserving_setups"])
			,'logdir': 			'/scratch/users/ksb/fireworks/logs/'}

###########################################################################
def sherlock2Q(timeInHours):
	return {'_fw_name': 		'CommonAdapter'
			,'_fw_q_type': 		'SLURM'
			,'rocket_launch': 	'python $SCRATCH/fireworks/fireworks_env/lib/python2.7/site-packages/fireworks/scripts/rlaunch_run.py -w $SCRATCH/fireworks/my_fworker.yaml -l $SCRATCH/fireworks/my_launchpad.yaml  singleshot'
			,'nodes': 			1
			,'ntasks_per_node': 16
			,'walltime': 		pp.printTime(timeInHours)+':00'
			,'queue': 			'owners,iric,normal'
			,'pre_rocket': 		'export PATH=$PATH:$SCRATCH/fireworks/fireworks_env/bin/' # 
			,'logdir': 			'/scratch/users/ksb/fireworks/logs/'}
def sherlock2QGPAW(timeInHours):
	return {'_fw_name': 		'CommonAdapter'
			,'_fw_q_type': 		'SLURM'
			,'rocket_launch': 	'python $SCRATCH/fireworks/fireworks_env/lib/python2.7/site-packages/fireworks/scripts/rlaunch_run.py -w $SCRATCH/fireworks/my_fworker.yaml -l $SCRATCH/fireworks/my_launchpad.yaml  singleshot'
			,'nodes': 			1
			,'ntasks_per_node': 16
			,'walltime': 		pp.printTime(timeInHours)+':00'
			,'queue': 			'owners,iric,normal'
			,'pre_rocket': 		'export PATH=$PATH:$SCRATCH/fireworks/fireworks_env/bin/;source /scratch/PI/suncat/sw/env.bash;export OMP_NUM_THREADS=1;export PYTHONPATH=/scratch/users/ksb/gpaw/gpaw_sg15/lib/python2.7/site-packages:$PYTHONPATH;export PATH=/scratch/users/ksb/gpaw/gpaw_sg15/bin:$PATH;export GPAW_SETUP_PATH=/scratch/users/ksb/gpaw/gpaw_sg15/norm_conserving_setups' # 
			,'logdir': 			'/scratch/users/ksb/fireworks/logs/'}
###########################################################################
###########################################################################
def suncatQ(timeInHours):
	return {'_fw_name': 		'CommonAdapter'
			,'_fw_q_type': 		'LoadSharingFacility'
			,'rocket_launch':  	'python /nfs/slac/g/suncatfs/ksb/fireworks/fireworks_virtualenv/lib/python2.7/site-packages/fireworks/scripts/rlaunch_run.py -w /nfs/slac/g/suncatfs/${USER}/fireworks/my_fworker.yaml -l /nfs/slac/g/suncatfs/${USER}/fireworks/my_launchpad.yaml singleshot'
			,'nodes': 			1
			,'ntasks_per_node':	8
			,'walltime': 		pp.printTime(timeInHours)
			,'queue': 			'suncat'
			,'pre_rocket':  	'unset LS_COLORS;source /nfs/slac/g/suncatfs/sw/espv20/setupenv;setenv PYTHONPATH /nfs/slac/g/suncatfs/${USER}/fireworks/fireworks_virtualenv/lib/python2.7/site-packages:/nfs/slac/g/suncatfs/fireworks_scripts/standard_tasks/:${PYTHONPATH};setenv PATH /afs/slac.stanford.edu/package/lsf/9.1.2/linux2.6-glibc2.3-x86_64/bin:${PATH}'
			,'logdir': 			'/nfs/slac/g/suncatfs/ksb/fireworks/logs/'}
###########################################################################
def suncat2Q(timeInHours):
	return {'_fw_name': 		'CommonAdapter'
			,'_fw_q_type': 		'LoadSharingFacility'
			,'rocket_launch':  	'python /nfs/slac/g/suncatfs/ksb/fireworks/fireworks_virtualenv/lib/python2.7/site-packages/fireworks/scripts/rlaunch_run.py -w /nfs/slac/g/suncatfs/${USER}/fireworks/my_fworker.yaml -l /nfs/slac/g/suncatfs/${USER}/fireworks/my_launchpad.yaml singleshot'
			,'nodes': 			1
			,'ntasks_per_node':	12
			,'walltime': 		pp.printTime(timeInHours)
			,'queue': 			'suncat2'
			,'pre_rocket':  	'unset LS_COLORS;source /nfs/slac/g/suncatfs/sw/espv20/setupenv;setenv PYTHONPATH /nfs/slac/g/suncatfs/${USER}/fireworks/fireworks_virtualenv/lib/python2.7/site-packages:/nfs/slac/g/suncatfs/fireworks_scripts/standard_tasks/:${PYTHONPATH};setenv PATH /afs/slac.stanford.edu/package/lsf/9.1.2/linux2.6-glibc2.3-x86_64/bin:${PATH}'
			,'logdir': 			'/nfs/slac/g/suncatfs/ksb/fireworks/logs/'}
###########################################################################

###########################################################################
###########################################################################
# 							Name 		fworker 	qadapter 		launchdir
sherlock 		= Cluster('sherlock',	'sherlock',	sherlockQ,		'/scratch/users/ksb/fireworks/jobs/')
sherlockgpaw 	= Cluster('sherlock',	'sherlock',	sherlockQGPAW,	'/scratch/users/ksb/fireworks/jobs/')

sherlock2 		= Cluster('sherlock',	'sherlock',	sherlock2Q,		'/scratch/users/ksb/fireworks/jobs/')
sherlock2gpaw 	= Cluster('sherlock',	'sherlock',	sherlock2QGPAW,	'/scratch/users/ksb/fireworks/jobs/')

suncat 			= Cluster('suncat',		'suncat',	suncatQ,		'/nfs/slac/g/suncatfs/ksb/fireworks/jobs/')
suncat2 		= Cluster('suncat2',	'suncat',	suncat2Q,		'/nfs/slac/g/suncatfs/ksb/fireworks/jobs/')
