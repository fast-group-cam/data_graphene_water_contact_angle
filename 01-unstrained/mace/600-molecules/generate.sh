source ../../../00-common/software/venv/bin/activate
python ../../../00-common/software/scripts/delete-waters.py \
    ../800-molecules/equi_final.lammps-data \
    -o start.lammps-data \
    -N 600
deactivate
