#!/bin/bash 
bsub -W00:09 -qsuncat-test -n4 -N -oomyjob.out -emyjob.err "$(pwd)/$1"

