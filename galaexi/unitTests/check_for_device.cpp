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
================================================================================================================================*/\

#include <iostream>
#include "../src/device/device.h"

/**
 * @brief Test program to test for the existence of functioning devices
 * 
 * The result of this program is used to by CMake to know whether
 * there are available accelerators to run accelerator related tests.
 * 
 * @returns 1 if working devices found, 0 otherwise.
 */
int main()
{
    int deviceCount, device;
    int gpuDeviceCount = 0;

#if (USE_ACCEL == ACCEL_CUDA)
    cudaDeviceProp properties;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    if (err != cudaSuccess)
    {
        std::cout << "Call to cudaGetDeviceCount FAILED with error: " << cudaGetErrorString(err) << std::endl;
        return 0;
    }

    // Machines with no GPUs can still report one emulation device
    for (device = 0; device < deviceCount; ++device)
    {
        DEVICE_ERR_CHECK( cudaGetDeviceProperties(&properties, device) );
        if (properties.major != 9999) // 9999 means emulation only
            ++gpuDeviceCount;
    }
#elif (USE_ACCEL == ACCEL_HIP)
    hipDeviceProp_t properties;
    hipError_t err = hipGetDeviceCount(&deviceCount);
    if (err != hipSuccess)
    {
        std::cout << "Call to hipGetDeviceCount FAILED with error: " << hipGetErrorString(err) << std::endl;
        return 0;
    }

    // Machines with no GPUs can still report one emulation device
    for (device = 0; device < deviceCount; ++device)
    {
        DEVICE_ERR_CHECK( hipGetDeviceProperties(&properties, device) );
        if (properties.major != 9999) // 9999 means emulation only
            ++gpuDeviceCount;
    }
#endif

    std::cout << gpuDeviceCount << " Devices found." << std::endl;
    // Return code to CMake
    if (gpuDeviceCount > 0)
        return 1;
    else
        return 0;
}
