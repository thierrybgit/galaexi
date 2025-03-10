#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Unit test 'TestApplyJacobian'
!> Test the routines from module: 'MOD_ApplyJacobian'.
!> @author Spencer Starr - 08.2024
!==================================================================================================================================

MODULE TestApplyJacobian
USE MOD_ApplyJacobian
USE MOD_Device
USE MOD_DG_Vars, ONLY: U, d_U
USE MOD_PreProc
USE MOD_Mesh_Vars, ONLY: nElems, sJ, d_sJ
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_DeviceMem
#endif
IMPLICIT NONE

CONTAINS

    !==================================================================================================================================
    !> @brief Initialize the domain for the test
    !==================================================================================================================================
    SUBROUTINE InitTest

        nElems = 16

        ! Allocate and initialize host side arrays
        ALLOCATE( U(PP_nVar,0:1,0:1,0:1,nElems) )
        ALLOCATE( sJ(0:1,0:1,0:1,nElems,0:FV_SIZE) )

        U(:,:,:,:,:) = 10.d0
        sJ(:,:,:,:,:) = 2.d0

        ! If this is a device test, allocate and copy to device
#if (USE_ACCEL != ACCEL_OFF)
        CALL AllocateDeviceMemory(d_U, SIZE_C_DOUBLE, SIZE(U))
        CALL CopyToDevice(d_U, U, SIZE(U))

        CALL AllocateDeviceMemory(d_sJ, SIZE_C_DOUBLE, SIZE(sJ))
        CALL CopyToDevice(d_sJ, sJ, SIZE(sJ))
#endif

    END SUBROUTINE InitTest

    !==================================================================================================================================
    !> @brief Test ApplyJacobian method
    !==================================================================================================================================
    SUBROUTINE TestMethod

    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER :: iElem,i,j,k,n
    !=========================================================================================

        ! Test WITHOUT toPhysical
        CALL ApplyJacobian(PP_nVar, 1, U, d_U, .FALSE., streamID=STREAM_DEFAULT)

#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyFromDevice(U, d_U, size(U))
#endif

        DO iElem=1,nElems
            DO k=0,1; DO j=0,1; DO i=0,1
                DO n=1,PP_nVar
                    IF (.NOT. ALMOSTEQUALABSORREL(U(n,i,j,k,iElem),5.d0,200*PP_RealTolerance)) THEN
                        PRINT '(A,i1,A,3(i2),A,i2,A,f12.6,A)', "Value for variable ",n," in U for point",i,j,k, " in element ", &
                                                            iElem, " is: ", U(n,i,j,k,iElem), " should be: 5.0."
                        ERROR STOP
                    END IF
                END DO
            END DO; END DO; END DO
        END DO

        ! Test WITH toPhysical
        CALL ApplyJacobian(PP_nVar, 1, U, d_U, .TRUE., streamID=STREAM_DEFAULT)

#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyFromDevice(U, d_U, size(U))
#endif

        DO iElem=1,nElems
            DO k=0,1; DO j=0,1; DO i=0,1
                DO n=1,PP_nVar
                    IF (.NOT. ALMOSTEQUALABSORREL(U(n,i,j,k,iElem),10.d0,200*PP_RealTolerance)) THEN
                        PRINT '(A,i1,A,3(i2),A,i2,A,f12.6,A)', "Value for variable ",n," in U for point",i,j,k, " in element ", &
                                                            iElem," is: ",U(n,i,j,k,iElem)," should be: 10.0."
                        ERROR STOP
                    END IF
                END DO
            END DO; END DO; END DO
        END DO

    END SUBROUTINE TestMethod

END MODULE TestApplyJacobian


!==================================================================================================================================
!> @brief Test driver program
!==================================================================================================================================
PROGRAM ApplyJacobian_TestDriver
USE TestApplyJacobian
IMPLICIT NONE

    CALL InitTest
    CALL TestMethod

END PROGRAM ApplyJacobian_TestDriver
