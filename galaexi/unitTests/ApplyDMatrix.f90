#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Unit test 'ApplyDMatrix'
!> Test the routines from module: 'MOD_Flux'.
!> @author Spencer Starr - 08.2024
!==================================================================================================================================

MODULE TestApplyDMatrix
USE ISO_C_BINDING
USE MOD_DG_Vars, ONLY: Ut, F, G, H, D_Hat_T
USE MOD_Mesh_Vars, ONLY: nElems
USE MOD_ApplyDMatrix
USE MOD_PreProc
USE MOD_Device
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_DeviceMem
USE MOD_DG_Vars, ONLY: d_Ut, d_F, d_G, d_H, d_D_Hat_T
#endif
IMPLICIT NONE

CONTAINS

    !==================================================================================================================================
    !> @brief Initialize the domain for the test
    !==================================================================================================================================
    SUBROUTINE InitTest

        nElems = 16

        ! Allocate and init host side arrays - use 16 5th order elements
        ALLOCATE( Ut(PP_nVar,0:1,0:1,0:1,1:nElems) )
        ALLOCATE( F(PP_nVar,0:1,0:1,0:1,1:nElems) )
        ALLOCATE( G(PP_nVar,0:1,0:1,0:1,1:nElems) )
        ALLOCATE( H(PP_nVar,0:1,0:1,0:1,1:nElems) )
        ALLOCATE( D_Hat_T(0:1,0:1) )

        F(DENS,:,:,:,:) = 1.0
        F(MOM1,:,:,:,:) = 2.0
        F(MOM2,:,:,:,:) = 3.0
        F(MOM3,:,:,:,:) = 4.0
        F(ENER,:,:,:,:) = 5.0

        G(DENS,:,:,:,:) = 1.0
        G(MOM1,:,:,:,:) = 2.0
        G(MOM2,:,:,:,:) = 3.0
        G(MOM3,:,:,:,:) = 4.0
        G(ENER,:,:,:,:) = 5.0

        H(DENS,:,:,:,:) = 1.0
        H(MOM1,:,:,:,:) = 2.0
        H(MOM2,:,:,:,:) = 3.0
        H(MOM3,:,:,:,:) = 4.0
        H(ENER,:,:,:,:) = 5.0

        D_Hat_T = 1.0
        Ut = 0.0

        
        ! If device test, allocate and copy on device
#if (USE_ACCEL != ACCEL_OFF)
        CALL AllocateDeviceMemory(d_Ut, SIZE_C_DOUBLE, SIZE(Ut))
        CALL CopyToDevice(d_Ut, Ut, SIZE(Ut))

        CALL AllocateDeviceMemory(d_F, SIZE_C_DOUBLE, SIZE(F))
        CALL CopyToDevice(d_F, F, SIZE(F))

        CALL AllocateDeviceMemory(d_G, SIZE_C_DOUBLE, SIZE(G))
        CALL CopyToDevice(d_G, G, SIZE(G))

        CALL AllocateDeviceMemory(d_H, SIZE_C_DOUBLE, SIZE(H))
        CALL CopyToDevice(d_H, H, SIZE(H))

        CALL AllocateDeviceMemory(d_D_Hat_T, SIZE_C_DOUBLE, SIZE(D_Hat_T))
        CALL CopyToDevice(d_D_Hat_T, D_Hat_T, SIZE(D_Hat_T))
#endif

    END SUBROUTINE InitTest


    !==================================================================================================================================
    !> @brief Test ApplyDMatrix methods
    !==================================================================================================================================
    SUBROUTINE TestMethod
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER :: iElem,i,j,k,n
    REAL    :: ans(PP_nVar)
    !=========================================================================================

        ! Invoke the method
        CALL ApplyDMatrix(1,nElems, Ut, F, G, H, D_Hat_T, streamID=STREAM_DEFAULT, doOverwrite=.TRUE.)

#if (USE_ACCEL != ACCEL_OFF)
        ! Copy the solution from the device
        CALL CopyFromDevice(Ut, d_Ut, SIZE(Ut))
#endif

        ans = [6.0, 12.0, 18.0, 24.0, 30.0]

        ! Check the solution based on which method was selected
        DO iElem = 1,16
            DO k = 0,1; DO j = 0,1; DO i = 0,1
                DO n = 1,PP_nVar
                    IF (.NOT. ALMOSTEQUALABSORREL(Ut(n,i,j,k,iElem),ans(n),200*PP_RealTolerance)) THEN
                        PRINT *, "Value for Ut from ApplyDMatrix for variable: ", n, " in DOF: ",i,j,k, &
                                 " of element: ", iElem, " is: ", Ut(n,i,j,k,iElem), " should be: ", ans(n)
                        ERROR STOP
                    END IF
                END DO
            END DO; END DO; END DO
        END DO

    END SUBROUTINE TestMethod


END MODULE TestApplyDMatrix



!==================================================================================================================================
!> @brief Test driver program
!==================================================================================================================================
PROGRAM ApplyDMatrix_TestDriver
USE TestApplyDMatrix
IMPLICIT NONE

    CALL InitTest
    CALL TestMethod

END PROGRAM ApplyDMatrix_TestDriver
