source ../../../00-common/software/venv/bin/activate
python ../../../00-common/software/scripts/flatten-graphene.py \
    ../2000-molecules/equi_final.lammps-data \
    -o start.lammps-data
deactivate
