!=================================================================================================================================
! Copyright (c) 2022-2024 Prof. Andrea Beck
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
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
!> \brief Fortran interfaces to C++ device memory methods in device_mem.cpp
!==================================================================================================================================
MODULE MOD_DeviceMem
USE MOD_Device
USE ISO_C_BINDING
IMPLICIT NONE

PRIVATE

!----------------------------------------------------------------------------------------------------------------------------------
! PUBLIC INTERFACES
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: AllocateDeviceMemory
PUBLIC :: CopyToDevice
PUBLIC :: CopyFromDevice
PUBLIC :: FreeDeviceMemory
PUBLIC :: CheckDeviceMemoryState
PUBLIC :: SIZE_C_INT, SIZE_C_DOUBLE

! Values for readable passing of variable sizes to AllocateDeviceMemory
! These are declared in as integer variables becuase the current implementation
! of Fortran enums doesn't allow any other type that INTEGER(4), and we need 
! INTEGER(C_SIZE_T)
INTEGER(C_SIZE_T), PARAMETER :: SIZE_C_INT = 4_C_SIZE_T
INTEGER(C_SIZE_T), PARAMETER :: SIZE_C_DOUBLE = 8_C_SIZE_T

INTERFACE
    SUBROUTINE AllocateDeviceMemory_Device(dVarKey, typeSize_bytes, arraySize) &
                bind(C, name="AllocateDeviceMemory_Device")
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTEGER(C_INT),VALUE     :: dVarKey
        INTEGER(C_SIZE_T), VALUE :: typeSize_bytes
        INTEGER(C_INT), VALUE    :: arraySize
    END SUBROUTINE AllocateDeviceMemory_Device
END INTERFACE

INTERFACE CopyToDevice
    SUBROUTINE CopyToDeviceGeneric(dVarKey, hostMemory, typeSize_bytes, arraySize) &
                bind(C, name="CopyToDevice")
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTEGER(C_INT),VALUE   :: dVarKey
        TYPE(C_PTR),VALUE   :: hostMemory
        INTEGER(C_SIZE_T), VALUE :: typeSize_bytes
        INTEGER(C_INT), VALUE  :: arraySize
    END SUBROUTINE CopyToDeviceGeneric
END INTERFACE CopyToDevice

INTERFACE CopyFromDevice
    SUBROUTINE CopyFromDeviceGeneric(hostMemory, dVarKey, typeSize_bytes, arraySize) &
        bind(C, name="CopyFromDevice")
    USE ISO_C_BINDING
    IMPLICIT NONE
        TYPE(C_PTR),VALUE     :: hostMemory
        INTEGER(C_INT),VALUE  :: dVarKey
        INTEGER(C_SIZE_T), VALUE :: typeSize_bytes
        INTEGER(C_INT),VALUE  :: arraySize
    END SUBROUTINE CopyFromDeviceGeneric
END INTERFACE CopyFromDevice

INTERFACE
    SUBROUTINE FreeDeviceMemory(dVarKey) &
                bind(C, name="FreeDeviceMemory")
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTEGER(C_INT),VALUE :: dVarKey
    END SUBROUTINE FreeDeviceMemory
END INTERFACE

INTERFACE
    SUBROUTINE CheckDeviceMemoryState(FreeMemory, TotalMemory) &
                bind(C, name="CheckDeviceMemoryState")
    USE ISO_C_BINDING, ONLY: C_SIZE_T
    IMPLICIT NONE
        INTEGER(C_SIZE_T) :: FreeMemory
        INTEGER(C_SIZE_T) :: TotalMemory
    END SUBROUTINE CheckDeviceMemoryState
END INTERFACE

CONTAINS

    !==================================================================================================================================
    !> Middle man method for device memory allocation that generate key value and calls C backend
    !==================================================================================================================================
    SUBROUTINE AllocateDeviceMemory(dVarKey, typeSize_bytes, arraySize)
    USE ISO_C_BINDING
    IMPLICIT NONE
    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT / OUTPUT VARIABLES
    INTEGER(C_INT), INTENT(INOUT) :: dVarKey
    INTEGER(C_SIZE_T),INTENT(IN)  :: typeSize_bytes
    INTEGER(C_INT), INTENT(IN)    :: arraySize
    !==================================================================================================================================

        ! Generate a unique integer key for the variable being allocated
        dVarKey = GenerateDeviceKey()

        ! Call device API wrapper to allocate device memory and store device pointer
        CALL AllocateDeviceMemory_Device(dVarKey, typeSize_bytes, arraySize)
    
    END SUBROUTINE AllocateDeviceMemory

END MODULE MOD_DeviceMem
