source ../../../00-common/software/venv/bin/activate
python ../../../00-common/software/scripts/bond-waters.py \
    ../../mace/fixed-2000-molecules/equi_final.lammps-data \
    --fix \
    -o start.lammps-data
echo "Added bonds to water molecules"
deactivate
