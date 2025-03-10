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

#include <stddef.h>
#include "device.h"

std::unordered_map<int, void*> DeviceVars;

// Declare memory management functions as external so they can bind to Fortran interfaces
extern "C"
{
    void AllocateDeviceMemory_Device(int dVarKey, size_t typeSize_bytes, int arraySize);
    //void CopyDoubleToDevice(int dVarKey, double* hostMemory, int arraySize);
    //void CopyIntToDevice(int dVarKey, int* hostMemory, int arraySize);
    void CopyToDevice(int dVarKey, void* hostMemory, size_t typeSize_bytes, int arraySize);
    //void CopyDoubleFromDevice(double* hostMemory,  int dVarKey, int arraySize);
    //void CopyIntFromDevice(int* hostMemory,  int dVarKey, int arraySize);
    void CopyFromDevice(void* hostMemory,  int dVarKey, size_t typeSize_bytes, int arraySize);
    void FreeDeviceMemory(int dVarKey);
    void CheckDeviceMemoryState(size_t* FreeMemory, size_t* TotalMemory);
}

/***********************************************************************************************************************
 * @brief Allocate a block of memory on the device and store the pointer in the device vars map using the requested key.
 * @param dVarKey Desired device variable key
 * @param typeSize_bytes Size in bytes of the data type to be allocated
 * @param arraySize Number of elements of type to allocate (1 for scalars)
***********************************************************************************************************************/
void AllocateDeviceMemory_Device(int dVarKey, size_t typeSize_bytes, int arraySize)
{
    void* d_arr;

#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
    // NVIDIA GH chips allocate host managed memory on "first touch", meaning whichever chip (CPU or GPU) uses it first gets it.
    // Because we have almost certainly initialized the memory on the CPU, it has touched it and we need to explicitly
    // allocate and copy it over to the GPU now. So only use pure hybrid handling for AMD chips.
    d_arr = NULL;
#else /* USE_HYBRID */
#if (USE_ACCEL == ACCEL_CUDA)
    DEVICE_ERR_CHECK( cudaMalloc(&d_arr, typeSize_bytes*arraySize) );
#elif (USE_ACCEL == ACCEL_HIP)
    DEVICE_ERR_CHECK( hipMalloc(&d_arr, typeSize_bytes*arraySize) );
#endif
#endif /* USE_HYBRID */

    DeviceVars[dVarKey] = d_arr;
}

/***********************************************************************************************************************
 * @brief Copy double host data to the device
 * @param dVarKey Device variable key
 * @param hostMemory Pointer to host copy of the data
 * @param arraySize Number of elements allocated for the block being copied (1 for scalars)
***********************************************************************************************************************/
void CopyDoubleToDevice(int dVarKey, double* hostMemory, int arraySize)
{
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
    // See note in AllocatedDeviceMemory_Device on NVIDIA hybrid memory copies
    DeviceVars[dVarKey] = (void*)hostMemory;
#else /* USE_HYBRID */
#if (USE_ACCEL == ACCEL_CUDA)
    DEVICE_ERR_CHECK( cudaMemcpy(DeviceVars[dVarKey], hostMemory, arraySize*sizeof(double), cudaMemcpyHostToDevice) );
#elif (USE_ACCEL == ACCEL_HIP)
    DEVICE_ERR_CHECK( hipMemcpy(DeviceVars[dVarKey], hostMemory, arraySize*sizeof(double), hipMemcpyHostToDevice) );
#endif
#endif /* USE_HYBRID */
}

/***********************************************************************************************************************
 * @brief Copy integer host data to the device
 * @param dVarKey Device variable key
 * @param hostMemory Pointer to host copy of the data
 * @param arraySize Number of elements allocated for the block being copied (1 for scalars)
***********************************************************************************************************************/
void CopyIntToDevice(int dVarKey, int* hostMemory, int arraySize)
{
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
    // See note in AllocatedDeviceMemory_Device on NVIDIA hybrid memory copies
    DeviceVars[dVarKey] = (void*)hostMemory;
#else /* USE_HYBRID */
#if (USE_ACCEL == ACCEL_CUDA)
    DEVICE_ERR_CHECK( cudaMemcpy(DeviceVars[dVarKey], hostMemory, arraySize*sizeof(int), cudaMemcpyHostToDevice) );
#elif (USE_ACCEL == ACCEL_HIP)
    DEVICE_ERR_CHECK( hipMemcpy(DeviceVars[dVarKey], hostMemory, arraySize*sizeof(int), hipMemcpyHostToDevice) );
#endif
#endif /* USE_HYBRID */
}

void CopyToDevice(int dVarKey, void* hostMemory, size_t typeSize_bytes, int arraySize)
{
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
    // See note in AllocatedDeviceMemory_Device on NVIDIA hybrid memory copies
    DeviceVars[dVarKey] = (void*)hostMemory;
#else /* USE_HYBRID */
#if (USE_ACCEL == ACCEL_CUDA)
    DEVICE_ERR_CHECK( cudaMemcpy(DeviceVars[dVarKey], hostMemory, typeSize_bytes*arraySize, cudaMemcpyHostToDevice) );
#elif (USE_ACCEL == ACCEL_HIP)
    DEVICE_ERR_CHECK( hipMemcpy(DeviceVars[dVarKey], hostMemory, typeSize_bytes*arraySize, hipMemcpyHostToDevice) );
#endif
#endif /* USE_HYBRID */
}

/***********************************************************************************************************************
 * @brief Copy double data from the device to the host
 * @param hostMemory Pointer to host copy of the data
 * @param dVarKey Device variable key
 * @param arraySize Number of elements allocated for the block being copied (1 for scalars)
***********************************************************************************************************************/
void CopyDoubleFromDevice(double* hostMemory, int dVarKey, int arraySize)
{
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
    // For NVIDIA hybrid architectures, memory has to be explicitly copied back.
    // For AMD hybrid architectures, we have to call synch API the for host to see device changes to memory
    DEVICE_ERR_CHECK( hipDeviceSynchronize() )
#else /* USE_HYBRID */
#if (USE_ACCEL == ACCEL_CUDA)
    DEVICE_ERR_CHECK( cudaMemcpy(hostMemory, DeviceVars[dVarKey], arraySize*sizeof(double), cudaMemcpyDeviceToHost) );
#elif (USE_ACCEL == ACCEL_HIP)
    DEVICE_ERR_CHECK( hipMemcpy(hostMemory, DeviceVars[dVarKey], arraySize*sizeof(double), hipMemcpyDeviceToHost) );
#endif
#endif /* USE_HYBRID */
}

/***********************************************************************************************************************
 * @brief Copy integer data from the device to the host
 * @param hostMemory Pointer to host copy of the data
 * @param dVarKey Device variable key
 * @param arraySize Number of elements allocated for the block being copied (1 for scalars)
***********************************************************************************************************************/
void CopyIntFromDevice(int* hostMemory, int dVarKey, int arraySize)
{
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
    // For NVIDIA hybrid architectures, memory has to be explicitly copied back.
    // For AMD hybrid architectures, we have to call synch API the for host to see device changes to memory
    DEVICE_ERR_CHECK( hipDeviceSynchronize() )
#else /* USE_HYBRID */
#if (USE_ACCEL == ACCEL_CUDA)
    DEVICE_ERR_CHECK( cudaMemcpy(hostMemory, DeviceVars[dVarKey], arraySize*sizeof(int), cudaMemcpyDeviceToHost) );
#elif (USE_ACCEL == ACCEL_HIP)
    DEVICE_ERR_CHECK( hipMemcpy(hostMemory, DeviceVars[dVarKey], arraySize*sizeof(int), hipMemcpyDeviceToHost) );
#endif
#endif /* USE_HYBRID */
}

void CopyFromDevice(void* hostMemory, int dVarKey, size_t typeSize_bytes, int arraySize)
{
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
    // For NVIDIA hybrid architectures, memory has to be explicitly copied back.
    // For AMD hybrid architectures, we have to call synch API the for host to see device changes to memory
    DEVICE_ERR_CHECK( hipDeviceSynchronize() )
#else /* USE_HYBRID */
#if (USE_ACCEL == ACCEL_CUDA)
    DEVICE_ERR_CHECK( cudaMemcpy(hostMemory, DeviceVars[dVarKey], typeSize_bytes*arraySize, cudaMemcpyDeviceToHost) );
#elif (USE_ACCEL == ACCEL_HIP)
    DEVICE_ERR_CHECK( hipMemcpy(hostMemory, DeviceVars[dVarKey], typeSize_bytes*arraySize, hipMemcpyDeviceToHost) );
#endif
#endif /* USE_HYBRID */
}

/***********************************************************************************************************************
 * @brief Deallocate memory on the device
 * @param dVarKey Device variable key
***********************************************************************************************************************/
void FreeDeviceMemory(int dVarKey)
{

#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
    // Do nothing. We want the Fortran code (Host) to control allocation of memory in this case.
    // In fact, we explicitly do NOT want to free memory here, because we would lose the Host owned allocation too.
#else
#if (USE_ACCEL == ACCEL_CUDA)
    DEVICE_ERR_CHECK( cudaFree(DeviceVars[dVarKey]) );
#elif (USE_ACCEL == ACCEL_HIP)
    DEVICE_ERR_CHECK( hipFree(DeviceVars[dVarKey]) );
#endif
#endif /* USE_HYBRID */

    // Need to remove the key from the map to make sure that it isn't illegally referenced
    // We remove the key from the dictionary, regardless if we are using a hybrid chip.
    // We can still access the allocated memory from the Fortran side via the host pointer.
    DeviceVars.erase(dVarKey);
}

/***********************************************************************************************************************
 * @brief Check the current memory quota of the device
***********************************************************************************************************************/
void CheckDeviceMemoryState(size_t* FreeMemory, size_t* TotalMemory)
{
#if (USE_ACCEL == ACCEL_CUDA)
    DEVICE_ERR_CHECK( cudaMemGetInfo(FreeMemory, TotalMemory) );
#elif (USE_ACCEL == ACCEL_HIP)
    DEVICE_ERR_CHECK( hipMemGetInfo(FreeMemory, TotalMemory) );
#endif

}
