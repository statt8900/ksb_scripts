#!/bin/bash
source /scratch/PI/suncat/sw/env.bash
source ~/scripts/rc/RCsher.sh
cd /scratch/users/ksb/fireworks

export HOSTNAME=sherlock #These are needed in order to run qlaunch, which expects to be able to import standardScripts for some reason
export SCRATCH=/scratch/users/ksb

qlaunch -r rapidfire
fab suncat
