#!/bin/bash
#SBATCH -p iric,owners
#SBATCH --job-name=name
#SBATCH --output=opt.log
#SBATCH --error=err.log
#SBATCH --time=01:00:00
#SBATCH --qos=normal
#SBATCH --nodes=7
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=FAIL
#SBATCH  --mail-user=jgauth32@stanford.edu
#SBATCH --ntasks-per-node=16
#SBATCH -x sh-20-35,gpu-14-1

export OMP_NUM_THREADS=1
LD_LIBRARY_PATH=/share/sw/free/openmpi/1.10.2/intel/2016/lib:/share/sw/non-free/intel/2016/u1/compilers_and_libraries_2016.1.150/linux/mkl/lib/intel64:/share/sw/non-free/intel/2016/u1/compilers_and_libraries_2016.1.150/linux/compiler/lib/intel64_lin:/home/vossj/suncat/lib:$LD_LIBRARY_PATH
PATH=/home/vossj/suncat/vbin:/share/sw/free/openmpi/1.10.2/intel/2016/bin:$PATH
NTASKS=`echo $SLURM_TASKS_PER_NODE|tr '(' ' '|awk '{print $1}'`
NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`
NCPU=`echo " $NTASKS * $NNODES " | bc`
echo $VASP_PP_PATH
#mpirun -n $NCPU /home/vossj/suncat/vbin/vasp
mpirun -n $NCPU /scratch/groups/suncat/vasp_std_ringe_073117

