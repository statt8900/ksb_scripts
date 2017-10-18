#!/bin/bash 
bsub -W20:09 -qsuncat2 -n12 -N -oomyjob.out -emyjob.err "$(pwd)/$1"
