source ../../../00-common/software/venv/bin/activate
python ../../../00-common/software/scripts/delete-waters.py \
    ../3000-molecules/equi_final.lammps-data \
    -o start.lammps-data \
    -N 2500
deactivate
