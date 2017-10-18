from fabric.api import run, env

env.hosts = ['ksb@suncatls1.slac.stanford.edu']

def suncat(): run('/u/if/ksb/scripts/fireworks/suncat_launcher.sh')

#WARNING, this is not executed from /home/ksb/scripts... --- must be copied into /scratch/users/ksb/fireworks to take effect