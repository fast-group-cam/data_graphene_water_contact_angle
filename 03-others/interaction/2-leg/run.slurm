#!/bin/bash
#SBATCH --job-name=single-point
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --account=XXXXX
#SBATCH --partition=standard
#SBATCH --qos=short

# Load the relevant CP2K module
# Ensure OMP_NUM_THREADS is consistent with cpus-per-task above

module load cp2k/cp2k-9.1.0
module load PrgEnv-gnu cpe/21.09 mkl/19.5-281
module load eigen/3.4.0
module load gsl/2.7
module load cray-python
module load genmaskcpu
module unload cray-python

export OMP_NUM_THREADS=1
WDIR=$(pwd)
FILE=cp2k
EXEC=/work/y07/shared/apps/core/cp2k/cp2k-9.1.0/exe/ARCHER2/cp2k.popt


###############################

# Part added with bash through loops
srun ${EXEC} ${FILE}_sub.inp >> ${FILE}_sub.out || exit 1
