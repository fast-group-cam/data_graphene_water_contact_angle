source ../../00-common/software/venv/bin/activate
python ../../00-common/software/scripts/stretch-graphene.py \
    ../../01-unstrained/mace/fixed-1000-molecules/equi_final.lammps-data \
    -s="-0.02" \
    -o start.lammps-data
deactivate
