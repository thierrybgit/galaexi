#!/bin/bash

current=$PWD
module load LUMI PrgEnv-cray partition/G buildtools cray-hdf5-parallel rocm
module li

cd galaexi
rm -rf build;mkdir build
cd build
cmake --fresh -DLIBS_USE_MPI=OFF -DLIBS_USE_OPENMP=OFF -DLIBS_USE_ACCEL=AMD -DFLEXI_NODETYPE=GAUSS-LOBATTO   \
              -DFLEXI_SPLIT_DG=ON -DFLEXI_UNITTESTS=OFF -DFLEXI_PARABOLIC=OFF -DLIBS_BUILD_HDF5=OFF          \
              -DCMAKE_HIP_COMPILER="/opt/rocm-6.0.3/bin/amdclang++" -DCMAKE_HIP_ARCHITECTURES=gfx90a         \
              -DCMAKE_HIP_FLAGS="-Rpass-analysis=kernel-resource-usage" -DFLEXI_TESTCASE="taylorgreenvortex" \
              -DCMAKE_Fortran_COMPILER=ftn  ..
make clean
make -j
