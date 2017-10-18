#!/bin/bash
source /scratch/PI/suncat/sw/env.bash
source ~/scripts/rc/RCsher.sh
cd /scratch/users/ksb/fireworks

export HOSTNAME=sherlock #These are needed in order to run qlaunch, which expects to be able to import standardScripts for some reason
export SCRATCH=/scratch/users/ksb

#qlaunch -r rapidfire
#fab suncat

/scratch/users/ksb/fireworks/fireworks_env/bin/lpad admin maintain
lpad detect_lostruns --time 20000 --fizzle # falsely accuses??? maybe make daemon less frequent and time longer
lpad detect_unreserved --time 20000 --rerun


echo "DONE $(date)"