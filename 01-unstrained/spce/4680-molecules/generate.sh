source ../../../00-common/software/venv/bin/activate
python ../../../00-common/software/scripts/flatten-graphene.py \
    ../../mace/4680-molecules/equi_final.lammps-data \
    -o tmp_0.lammps-data
echo "Flattened graphene sheet"
python ../../../00-common/software/scripts/bond-waters.py \
    tmp_0.lammps-data \
    --fix \
    -o start.lammps-data
echo "Added bonds to water molecules"
rm tmp_0.lammps-data
deactivate
