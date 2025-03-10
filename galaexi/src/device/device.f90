!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!==================================================================================================================================
!> \brief Fortran interfaces to C++ device variables in device.h needed on the Fortran side
!>
!> The members of this module are those C++ device-side variables that require a Fortran copy. Typically this is for API calls,
!> passing data to/from the device.
!==================================================================================================================================
MODULE MOD_Device
USE ISO_C_BINDING
IMPLICIT NONE
PRIVATE
    
!----------------------------------------------------------------------------------------------------------------------------------
! PUBLIC VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: DeviceMemoryState
PUBLIC :: CurrentDeviceMemState
PUBLIC :: STREAM_DEFAULT, STREAM_LOW, STREAM_MID, STREAM_HIGH

PUBLIC :: DefineParametersDevice
#if (USE_ACCEL != ACCEL_OFF)
PUBLIC :: InitDevice
PUBLIC :: FinalizeDevice
PUBLIC :: GenerateDeviceKey
#endif

! Fortran interfaces for memory management methods in device_mem.cpp
!----------------------------------------------------------------------------------------------------------------------
! Mirrors the C struct of the same name defined in device.h
! Only used as the return value for CheckDeviceMemoryState
TYPE :: DeviceMemoryState
    INTEGER(C_SIZE_T) :: FreeMemory
    INTEGER(C_SIZE_T) :: TotalMemory
END TYPE

! Instance of the DeviceMemoryState type that can be used globally to check state
! and then keep it stored, so it can be checked without having to repoll if we that
! state hasn't changed.
TYPE(DeviceMemoryState) :: CurrentDeviceMemState

! Values for readable passing streams. These values are indexes for the streams array found in device.h/device.cpp
INTEGER, PARAMETER :: STREAM_DEFAULT = 0
INTEGER, PARAMETER :: STREAM_LOW = 1
INTEGER, PARAMETER :: STREAM_MID = 2
INTEGER, PARAMETER :: STREAM_HIGH = 3

! Counter for number of device variables, used to generate sequential list of integers for device pointer keys
INTEGER :: keyCount = 0

!----------------------------------------------------------------------------------------------------------------------------------
! INTERFACES TO GENERAL DEVICE METHODS IN device.cpp
!----------------------------------------------------------------------------------------------------------------------------------
#if (USE_ACCEL != ACCEL_OFF)
INTERFACE InitStreams
    SUBROUTINE InitStreams(useStreams) bind(C, name="InitStreams")
        USE ISO_C_BINDING
        IMPLICIT NONE
        LOGICAL, VALUE :: useStreams
    END SUBROUTINE InitStreams
END INTERFACE

INTERFACE FinalizeDevice
    SUBROUTINE FinalizeDevice() bind(C, name="FinalizeDevice")
        USE ISO_C_BINDING
        IMPLICIT NONE
    END SUBROUTINE FinalizeDevice
END INTERFACE
#endif

CONTAINS

!==================================================================================================================================
!> Define device input parameters
!==================================================================================================================================
SUBROUTINE DefineParametersDevice()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

    CALL prms%SetSection("Device")
    CALL prms%CreateLogicalOption('useStreams', "Set true to use concurrent streaming on device.", '.TRUE.')

END SUBROUTINE DefineParametersDevice

#if (USE_ACCEL != ACCEL_OFF)
!==================================================================================================================================
!> Initialize the device
!==================================================================================================================================
SUBROUTINE InitDevice()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: GETLOGICAL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL               :: useStreams
!==================================================================================================================================

    SWRITE(UNIT_stdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' INIT DEVICE...'

    useStreams = GETLOGICAL("useStreams")
    CALL InitStreams(useStreams)

    SWRITE(UNIT_stdOut,'(A)')' INIT DEVICE DONE!'
    SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitDevice

!==================================================================================================================================
!> Generate a device key value and then increment the counter
!==================================================================================================================================
INTEGER(C_INT) FUNCTION GenerateDeviceKey()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

    GenerateDeviceKey = keyCount
    keyCount = keyCount + 1

END FUNCTION GenerateDeviceKey
#endif

END MODULE MOD_Device