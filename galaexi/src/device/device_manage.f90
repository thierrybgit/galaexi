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
!> \brief Fortran interfaces to C++ methods used to issues commands to the device or change its operational state
!==================================================================================================================================
MODULE MOD_DeviceManage
USE MOD_Device
USE ISO_C_BINDING
IMPLICIT NONE

PRIVATE

!----------------------------------------------------------------------------------------------------------------------------------
! PUBLIC INTERFACES
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: SetDevice
PUBLIC :: SynchronizeDevice


INTERFACE
    SUBROUTINE SetDevice(deviceID) bind(C, name="SetDevice")
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTEGER(C_INT), VALUE    :: deviceID
    END SUBROUTINE SetDevice
END INTERFACE

INTERFACE
    SUBROUTINE SynchronizeDevice() bind(C, name="SynchronizeDevice")
    IMPLICIT NONE
    END SUBROUTINE SynchronizeDevice
END INTERFACE

END MODULE MOD_DeviceManage