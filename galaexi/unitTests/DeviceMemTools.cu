#include "device.h"

///================================================================================================================================
/// Unit test 'DeviceMemTools'
/// Test the routines from module: 'DeviceMem'.
/// This file holds only the device side code. Host side code in DeviceMemTools.f90
/// \author Spencer Starr - 06.2024
//================================================================================================================================

__device__ int devErrCode;

extern "C"
{
    int CheckDeviceVarValues(int realArrKey, int intArrKey);
    void ModifyValsOnDevice(int realArrKey, int intArrKey);
}

//----------------------------------------------------------------------------------------------------------------------
// KERNEL METHODS
//----------------------------------------------------------------------------------------------------------------------

/***********************************************************************************************************************
 * @brief Kernel for checking values were copied to the device correctly
 * @param realArr Pointer to device copy of the double float test array
 * @param intArr Pointer to the device copy of the int test array
***********************************************************************************************************************/
__global__ void CheckDeviceVarValues_Kernel(double* realArr, int* intArr)
{
    devErrCode = 0;

    if (realArr[4562] != 11.0)
    {
        devErrCode += 1;
    }

    if (intArr[999345] != 51)
    {
        devErrCode += 2;
    }
    
}

/***********************************************************************************************************************
 * @brief Modify the two test arrays on the device such that we can check the copy back to host is working.
 * @param realArr Pointer to device copy of the double float test array
 * @param intArr Pointer to the device copy of the int test array 
***********************************************************************************************************************/
__global__ void ModifyValsOnDevice_Kernel(double* realArr, int* intArr)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    realArr[i] = 25.0;
    intArr[i] = 45;
}


//----------------------------------------------------------------------------------------------------------------------
// INTERFACES TO FORTRAN METHODS
//----------------------------------------------------------------------------------------------------------------------

/***********************************************************************************************************************
 * @brief C interface to method for checking values were copied to the device correctly.
 * @param realArrKey Key for the double float test array
 * @param intArrKey Key for the int test array
 * @return Error code from kernel
***********************************************************************************************************************/
int CheckDeviceVarValues(int realArrKey, int intArrKey)
{

    INVOKE_KERNEL( CheckDeviceVarValues_Kernel, 1, 1, 0, 0, (double*)DeviceVars[realArrKey], (int*)DeviceVars[intArrKey] );

    // Get the error code from the device
    int hostErrCode = 999;
#if (USE_ACCEL == ACCEL_CUDA)
    DEVICE_ERR_CHECK( cudaMemcpyFromSymbol(&hostErrCode, devErrCode, sizeof(int)) );
#elif (USE_ACCEL == ACCEL_HIP)
    DEVICE_ERR_CHECK( hipMemcpyFromSymbol(&hostErrCode, HIP_SYMBOL(devErrCode), sizeof(int), 0) );
#endif

    return hostErrCode;
}

/***********************************************************************************************************************
 * @brief C interface to method for modifying values on the device in order to test copy back to host.
 * @param realArrKey Key for the double float test array
 * @param intArrKey Key for the int test array
***********************************************************************************************************************/
void ModifyValsOnDevice(int realArrKey, int intArrKey)
{
    INVOKE_KERNEL( ModifyValsOnDevice_Kernel, 1, 100, 0, 0, (double*)DeviceVars[realArrKey], (int*)DeviceVars[intArrKey] );
}
