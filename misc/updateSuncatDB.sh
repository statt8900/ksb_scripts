#!/bin/bash
rsync -rtqvu -e "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null" /scratch/users/ksb/db ksb@suncatls2.slac.stanford.edu:/nfs/slac/g/suncatfs/ksb/
echo "DONE $(date)"