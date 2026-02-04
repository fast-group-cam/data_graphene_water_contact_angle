echo ""
ulimit -s unlimited

source ../../../00-common/software/venv/bin/activate

for root_folder in "s+0.00" "s+1.00" "s+2.00"; do
    for dist in `seq 2.50 0.25 4.50`; do
        cd "${root_folder}/d+${dist}"
        cp ../../cp2k.inp .
        ase convert geometry.xyz init.pdb

        dims=$(awk 'match($0,/Lattice="([^"]+)"/,m){ print m[1]; exit }' geometry.xyz)
        IFS=' ' read -r -a A <<< "$dims"
        DIM1=${A[0]}
        DIM2=${A[4]}
        DIM3=${A[8]}
        sed -e "s/XXXDIM1XXX/${DIM1}/" -e "s/XXXDIM2XXX/${DIM2}/" -e "s/XXXDIM3XXX/${DIM3}/" cp2k.inp > cp2k_sub.inp

        cd ../..
    done
done

deactivate

