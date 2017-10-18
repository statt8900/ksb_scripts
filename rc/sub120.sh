#!/bin/bash 
bsub -W20:09 -qsuncat -n8 -N -oomyjob.out -emyjob.err "$(pwd)/$1"
