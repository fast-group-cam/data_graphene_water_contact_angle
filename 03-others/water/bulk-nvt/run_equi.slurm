#!/bin/bash
#SBATCH --job-name=water-nvt
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=02:30:00
#SBATCH --account=XXXXX
#SBATCH --partition=standard
#SBATCH --qos=standard

# Ensure OMP_NUM_THREADS is consistent with cpus-per-task above
module load PrgEnv-gnu
module load cray-python
export SLURM_CPU_FREQ_REQ=2250000
LMP_EXE=../../../../00-common/software/lammps/lmp

rm -rf run
mkdir run
cd run

echo "Run started at $(date '+%Y-%m-%d %H:%M:%S')"
srun --overlap --oversubscribe --ntasks=128 --unbuffered $LMP_EXE -in ../inp_equi.lmp
rm log.lammps
cd ..
mv run/nvt_equi_final.lammps-data equi_final.lammps-data
echo "Run completed at $(date '+%Y-%m-%d %H:%M:%S')"
