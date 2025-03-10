/*================================================================================================================================
  Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
  Copyright (c) 2022-2024 Prof. Andrea Beck
  This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
  For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
 
  FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 
  FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
 
  You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
================================================================================================================================*/

// Declare device management functions as external so they can bind to Fortran interfaces
extern "C"
{
    void SetDevice(int deviceID);
    void SynchronizeDevice();
}


/***********************************************************************************************************************
 * @brief Wrapper for the set device methods in CUDA and HIP. For now just used to create device contexts in unit tests.
 * @param deviceID ID of the device being chosen
***********************************************************************************************************************/
void SetDevice(int deviceID)
{
#if (USE_ACCEL == ACCEL_CUDA)
        DEVICE_ERR_CHECK( cudaSetDevice(deviceID) )
#elif (USE_ACCEL == ACCEL_HIP)
        DEVICE_ERR_CHECK( hipSetDevice(deviceID) )
#endif
}

/***********************************************************************************************************************
 * @brief Wrapper for device synchronize methods
***********************************************************************************************************************/
void SynchronizeDevice()
{
#if (USE_ACCEL == ACCEL_CUDA)
        DEVICE_ERR_CHECK( cudaPeekAtLastError() )
        DEVICE_ERR_CHECK( cudaDeviceSynchronize() )
#elif (USE_ACCEL == ACCEL_HIP)
        DEVICE_ERR_CHECK( hipPeekAtLastError() )
        DEVICE_ERR_CHECK( hipDeviceSynchronize() )
#endif
}
