#pragma once
/*================================================================================================================================
  Copyright (c) 2022-2024 Prof. Andrea Beck
  Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
  This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
  For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
 
  FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 
  FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
 
  You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
================================================================================================================================*/

#include <unordered_map>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#if USE_MPI
#include "mpi.h"
#endif
#include "flexi.h"

// Include the runtime header of the vendor programming model
#if (USE_ACCEL == ACCEL_CUDA)
    #include <cuda_runtime.h>
#elif (USE_ACCEL == ACCEL_HIP)
    #include <hip/hip_runtime.h>
#ifdef USE_NVTX
    #include <nvtx3/nvToolsExt.h>
#endif
#endif /* USE_ACCEL == ACCEL_CUDA */

// Global variables for device memory management
//----------------------------------------------------------------------------------------------------------------------------------
// This hash map stores pointers to all variables in device memory
// It is keyed by and integer value stored in a variable with the name of the device variable with the pattern:
// host variable       --> var
// device variable key --> d_var
// Instantiated in device_mem.cpp
// Only has a scope within libflexif90 compile unit (can't be used in the entry point method from flexi.f90)
extern std::unordered_map<int, void*> DeviceVars;

// Simple struct that stores return values of the a call to cudaGetMemInfo
// May expand later to accomodate more device information
// Mirrored on host side by Fortran type declared in device.f90
struct DeviceMemoryState
{
    size_t FreeMemory;
    size_t TotalMemory;
};

// Declarations for streams
//----------------------------------------------------------------------------------------------------------------------------------
#if (USE_ACCEL == ACCEL_CUDA)
#define STREAM_T cudaStream_t
#elif (USE_ACCEL == ACCEL_HIP)
#define STREAM_T hipStream_t
#else
#define STREAM_T int
#endif
// streams[0] is the default stream
extern STREAM_T streams[7];


// Define macros to make device handling easier.
//----------------------------------------------------------------------------------------------------------------------------------
#if (USE_ACCEL == ACCEL_CUDA)
#define DEVICE_ERR_CHECK( call )                                                                                          \
{                                                                                                                  \
  cudaError_t err = call;                                                                                          \
  if ( cudaSuccess != err)                                                                                         \
  {                                                                                                                \
    fprintf(stderr, "CUDA error %d for %s in %d of %s : %s.\n", err, #call , __LINE__ , __FILE__ ,cudaGetErrorString(err)); \
    exit(1);                                                                                                         \
  }                                                                                                                  \
}

#define INVOKE_KERNEL( func, blocks, threads, shared, stream, ... ) func<<<blocks,threads,shared,stream>>>(__VA_ARGS__)
#elif (USE_ACCEL == ACCEL_HIP)
#define DEVICE_ERR_CHECK( call )                                                                                          \
{                                                                                                                  \
  hipError_t err = call;                                                                                          \
  if ( hipSuccess != err)                                                                                         \
  {                                                                                                                \
    fprintf(stderr, "HIP error %d for %s in %d of %s : %s.\n", err, #call , __LINE__ , __FILE__ ,hipGetErrorString(err)); \
    exit(1);                                                                                                         \
  }                                                                                                                  \
}

#define INVOKE_KERNEL( func, blocks, threads, shared, stream, ... ) func<<<dim3(blocks),dim3(threads),shared,stream>>>(__VA_ARGS__)
#endif
