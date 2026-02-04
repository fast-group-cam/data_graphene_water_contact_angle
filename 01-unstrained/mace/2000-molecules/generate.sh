source ../../../00-common/software/venv/bin/activate
python ../../../00-common/software/scripts/delete-waters.py \
    ../2500-molecules/equi_final.lammps-data \
    -o start.lammps-data \
    -N 2000
deactivate
