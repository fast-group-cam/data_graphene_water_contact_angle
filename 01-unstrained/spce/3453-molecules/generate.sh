source ../../../00-common/software/venv/bin/activate
python ../../../00-common/software/scripts/flatten-graphene.py \
    ../../mace/3500-molecules/equi_final.lammps-data \
    -o tmp_0.lammps-data
echo "Flattened graphene sheet"
python ../../../00-common/software/scripts/delete-waters.py \
    tmp_0.lammps-data \
    -o tmp_1.lammps-data \
    -N 3453
echo "Deleted water molecules"
python ../../../00-common/software/scripts/bond-waters.py \
    tmp_1.lammps-data \
    --fix \
    -o start.lammps-data
echo "Added bonds to water molecules"
rm tmp_0.lammps-data
rm tmp_1.lammps-data
deactivate
