Build LAMMPS (patched with symmetrix) in this directory.

Compile options depend on the computation architecture, but an example would be:

>    cd ..
>    git clone -b release https://github.com/lammps/lammps
>    git clone --recursive https://github.com/wcwitt/symmetrix
>    cd symmetrix/pair_symmetrix
>    ./install.sh ../../lammps
>    cd ../..
>    cd lammps
>    cmake \
>        -B build \
>        -D CMAKE_BUILD_TYPE=Release \
>        -D CMAKE_C_COMPILER=cc \
>        -D CMAKE_CXX_COMPILER=CC \
>        -D CMAKE_Fortran_COMPILER=ftn \
>        -D CMAKE_CXX_STANDARD=20 \
>        -D CMAKE_CXX_STANDARD_REQUIRED=ON \
>        -D CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -march=native -ffast-math" \
>        -D BUILD_SHARED_LIBS=ON \
>        -D PKG_KOKKOS=ON \
>        -D Kokkos_ENABLE_SERIAL=ON \
>        -D Kokkos_ENABLE_OPENMP=ON \
>        -D BUILD_OMP=ON \
>        -D Kokkos_ARCH_NATIVE=ON \
>        -D Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=ON \
>        -D SYMMETRIX_KOKKOS=ON \
>        cmake
>    cmake --build build -j 128
>    cd ../..

See https://github.com/wcwitt/symmetrix/blob/main/pair_symmetrix/README.md for more details.
