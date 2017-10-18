#!/bin/bash 
bsub -W10:09 -qsuncat -n8 -N -oomyjob.out -emyjob.err "$(pwd)/$1"