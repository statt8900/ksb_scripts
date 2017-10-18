#!/bin/bash 
bsub -W10:09 -qsuncat2 -n12 -N -oomyjob.out -emyjob.err "$(pwd)/$1"