#!/bin/bash
#SBATCH --job-name=slab-mace
#SBATCH --nodes=7
#SBATCH --ntasks=896
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=XXXXX
#SBATCH --partition=standard
#SBATCH --qos=standard

# Ensure OMP_NUM_THREADS is consistent with cpus-per-task above
module load PrgEnv-gnu
module load cray-python
export SLURM_CPU_FREQ_REQ=2250000
LMP_EXE=../../../../00-common/software/lammps/lmp

cd run

echo "Run started at $(date '+%Y-%m-%d %H:%M:%S')"
srun --overlap --oversubscribe --ntasks=896 --unbuffered $LMP_EXE -in ../../inp_prod.lmp
rm log.lammps
egrep '^(keywords:|data:$|---$|\.\.\.$|  - \[)' nvt_prod_log_0.lammps > nvt_prod_0.yaml
echo "Run completed at $(date '+%Y-%m-%d %H:%M:%S')"
