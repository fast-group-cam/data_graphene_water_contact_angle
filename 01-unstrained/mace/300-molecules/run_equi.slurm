#!/bin/bash
#SBATCH --job-name=lmp-mace
#SBATCH --nodes=4
#SBATCH --ntasks=512
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
source ../../../00-common/software/venv/bin/activate
LMP_EXE=../../../../00-common/software/lammps/lmp

rm -rf run_equi
mkdir run_equi
cd run_equi

echo "Run started at $(date '+%Y-%m-%d %H:%M:%S')"
srun --overlap --oversubscribe --ntasks=512 --unbuffered $LMP_EXE -in ../../inp_equi.lmp
echo "Run completed at $(date '+%Y-%m-%d %H:%M:%S')"

cd ..
cp run_equi/nvt_equi_final.lammps-data equi_final.lammps-data

