#!/bin/bash
#SBATCH --job-name=s-0.25
#SBATCH --nodes=8
#SBATCH --ntasks=1024
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
LMP_EXE=../../../00-common/software/lammps/lmp

RUN_NUM=11
cd run_prod

echo "Run started at $(date '+%Y-%m-%d %H:%M:%S')"
RESTART_NUM=$(( RUN_NUM * 80000 ))
sed -e "s/RESTART_NUM/${RESTART_NUM}/1;s/RUN_NUM/${RUN_NUM}/1" ../../inp_prod_continue_template.lmp > continue_${RUN_NUM}.lmp
srun --overlap --oversubscribe --ntasks=1024 --unbuffered $LMP_EXE -in continue_${RUN_NUM}.lmp
rm continue_${RUN_NUM}.lmp
echo "Run completed at $(date '+%Y-%m-%d %H:%M:%S')"
