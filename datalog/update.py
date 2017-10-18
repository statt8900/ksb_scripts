#!/scratch/PI/suncat/sw/bin/python
#SBATCH -p iric
#SBATCH --output=opt.log
#SBATCH --error=err.log
#SBATCH --time=01:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ksb@stanford.edu
#SBATCH --ntasks-per-node=16
#SBATCH -x sh-20-35,gpu-14-1

import manageSharedDatabase
.resetAll()
print 'updating '
manageSharedDatabase.updateDB()