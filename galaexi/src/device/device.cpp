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
#include "device.h"
#include <algorithm>

/***********************************************************************************************************************
 * @file device.cpp
 * @brief C++ implementations for GPU management method interfaces in device.f90
***********************************************************************************************************************/

STREAM_T streams[7];

extern "C"
{
   void InitStreams(bool useStreams);
   void FinalizeGPU(); 
}

/***********************************************************************************************************************
 * @brief Initialize device streams
***********************************************************************************************************************/
void InitStreams(bool useStreams)
{

    // Create the default stream. We do this even if there are no streams so we don't have to init streams fully
    // for unit testing.
    streams[0] = (STREAM_T)0;

    if (useStreams)
    {
        // Get available stream priority range of current device
        int min_priority, max_priority;
#if (USE_ACCEL == ACCEL_CUDA)
        DEVICE_ERR_CHECK( cudaDeviceGetStreamPriorityRange(&min_priority, &max_priority) );
        int mid_priority = std::min(max_priority+1, min_priority);

        DEVICE_ERR_CHECK( cudaStreamCreateWithPriority(&streams[1], cudaStreamDefault, min_priority) );
        DEVICE_ERR_CHECK( cudaStreamCreateWithPriority(&streams[2], cudaStreamDefault, mid_priority) );
        DEVICE_ERR_CHECK( cudaStreamCreateWithPriority(&streams[3], cudaStreamDefault, max_priority) );
        DEVICE_ERR_CHECK( cudaStreamCreateWithPriority(&streams[4], cudaStreamDefault, mid_priority) );
        DEVICE_ERR_CHECK( cudaStreamCreateWithPriority(&streams[5], cudaStreamDefault, mid_priority) );
        DEVICE_ERR_CHECK( cudaStreamCreateWithPriority(&streams[6], cudaStreamDefault, mid_priority) );
#elif (USE_ACCEL == ACCEL_HIP)
        DEVICE_ERR_CHECK( hipDeviceGetStreamPriorityRange(&min_priority, &max_priority) );
        int mid_priority = std::min(max_priority+1, min_priority);

        DEVICE_ERR_CHECK( hipStreamCreateWithPriority(&streams[1], hipStreamDefault, min_priority) );
        DEVICE_ERR_CHECK( hipStreamCreateWithPriority(&streams[2], hipStreamDefault, mid_priority) );
        DEVICE_ERR_CHECK( hipStreamCreateWithPriority(&streams[3], hipStreamDefault, max_priority) );
        DEVICE_ERR_CHECK( hipStreamCreateWithPriority(&streams[4], hipStreamDefault, mid_priority) );
        DEVICE_ERR_CHECK( hipStreamCreateWithPriority(&streams[5], hipStreamDefault, mid_priority) );
        DEVICE_ERR_CHECK( hipStreamCreateWithPriority(&streams[6], hipStreamDefault, mid_priority) );
#endif
    }
    else
    {
        streams[1] = streams[0];
        streams[2] = streams[0];
        streams[3] = streams[0];
        streams[4] = streams[0];
        streams[5] = streams[0];
        streams[6] = streams[0];
    }
}

/***********************************************************************************************************************
 * @brief Finalize device and destroy any initialize device streams
***********************************************************************************************************************/
void FinalizeGPU()
{
#if (USE_ACCEL == ACCEL_CUDA)
    DEVICE_ERR_CHECK( cudaStreamDestroy(streams[1]) );
    DEVICE_ERR_CHECK( cudaStreamDestroy(streams[2]) );
    DEVICE_ERR_CHECK( cudaStreamDestroy(streams[3]) );
    DEVICE_ERR_CHECK( cudaStreamDestroy(streams[4]) );
    DEVICE_ERR_CHECK( cudaStreamDestroy(streams[5]) );
    DEVICE_ERR_CHECK( cudaStreamDestroy(streams[6]) );
#elif (USE_ACCEL == ACCEL_HIP)
    DEVICE_ERR_CHECK( hipStreamDestroy(streams[1]) );
    DEVICE_ERR_CHECK( hipStreamDestroy(streams[2]) );
    DEVICE_ERR_CHECK( hipStreamDestroy(streams[3]) );
    DEVICE_ERR_CHECK( hipStreamDestroy(streams[4]) );
    DEVICE_ERR_CHECK( hipStreamDestroy(streams[5]) );
    DEVICE_ERR_CHECK( hipStreamDestroy(streams[6]) );
#endif
}
