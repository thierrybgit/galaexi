#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Unit test 'TestVAX'
!> Test the routines from module: 'MOD_Vector'.
!> @author Spencer Starr - 08.2024
!==================================================================================================================================

MODULE TestVAX
USE ISO_C_BINDING
USE MOD_Vector
USE MOD_PreProc
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_DeviceMem
#endif
IMPLICIT NONE

INTEGER :: nTotal
REAL,ALLOCATABLE :: VecIn(:)
REAL,ALLOCATABLE :: VecOut(:)
REAL,ALLOCATABLE :: VecOut2(:)
REAL :: Const1
REAL :: Const2

INTEGER(C_INT) :: d_VecIn 
INTEGER(C_INT) :: d_VecOut 
INTEGER(C_INT) :: d_VecOut2 
INTEGER(C_INT) :: d_Const1 
INTEGER(C_INT) :: d_Const2 

CONTAINS

    !==================================================================================================================================
    !> @brief Initialize the arrays for the test
    !==================================================================================================================================
    SUBROUTINE InitTest

        nTotal = 1000000
        ALLOCATE( VecIn(nTotal) )
        ALLOCATE( VecOut(nTotal) )
        ALLOCATE( VecOut2(nTotal) )

        Const1 = 2.d0
        Const2 = 3.d0

#if (USE_ACCEL != ACCEL_OFF)
        CALL AllocateDeviceMemory(d_VecIn, SIZE_C_DOUBLE, SIZE(VecIn))
        CALL AllocateDeviceMemory(d_VecOut, SIZE_C_DOUBLE, SIZE(VecOut))
        CALL AllocateDeviceMemory(d_VecOut2, SIZE_C_DOUBLE, SIZE(VecOut2))
#endif

    END SUBROUTINE InitTest

    !==================================================================================================================================
    !> @brief Test VAX methods
    !==================================================================================================================================
    SUBROUTINE TestMethods
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER :: i
    !=========================================================================================

        ! Test VAXPB_2STEP
        VecIn = 5.d0
        VecOut = 4.d0
        VecOut2 = 7.d0
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_VecIn, VecIn, SIZE(VecIn))
        CALL CopyToDevice(d_VecOut, VecOut, SIZE(VecOut))
        CALL CopyToDevice(d_VecOut2, VecOut2, SIZE(VecOut2))
#endif
        CALL VAXPB_2STEP(nTotal,VecOut,VecOut2,VecIn,d_VecOut,d_VecOut2,d_VecIn,Const1,Const2)
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyFromDevice(VecOut, d_VecOut, SIZE(VecOut))
        CALL CopyFromDevice(VecOut2, d_VecOut2, SIZE(VecOut2))
#endif
        DO i = 1,nTotal
            IF (.NOT. ALMOSTEQUALABSORREL(VecOut(i),13.d0,200*PP_RealTolerance)) THEN
                PRINT '(A,i6,A,F12.6,A)', "Test of VAXPB_2STEP FAILED for element:",i,". Value is:",VecOut(i),", should be: 13.0."
                ERROR STOP
            END IF

            IF (.NOT. ALMOSTEQUALABSORREL(VecOut2(i),46.d0,200*PP_RealTolerance)) THEN
                PRINT '(A,i6,A,F12.6,A)', "Test of VAXPB_2STEP FAILED for element:",i,". Value is:",VecOut2(i),", should be: 46.0."
                ERROR STOP
            END IF
        END DO
        
        ! Test VAXPB_IN
        VecOut = 4.d0
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_VecOut, VecOut, SIZE(VecOut))
#endif
        CALL VAXPB_IN(nTotal,VecOut,VecIn,d_VecOut,d_VecIn,Const1)
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyFromDevice(VecOut, d_VecOut, SIZE(VecOut))
#endif
        DO i = 1,nTotal
            IF (.NOT. ALMOSTEQUALABSORREL(VecOut(i),14.d0,200*PP_RealTolerance)) THEN
                PRINT '(A,i6,A,F12.6,A)', "Test of VAXPB_IN FAILED for element:",i,". Value is:",VecOut(i),", should be: 14.0."
                ERROR STOP
            END IF
        END DO

        ! Test VAXPB_OUT
        VecOut = 4.d0
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_VecOut, VecOut, SIZE(VecOut))
#endif
        CALL VAXPB_OUT(nTotal,VecOut,VecIn,d_VecOut,d_VecIn,Const1)
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyFromDevice(VecOut, d_VecOut, SIZE(VecOut))
#endif
        DO i = 1,nTotal
            IF (.NOT. ALMOSTEQUALABSORREL(VecOut(i),13.d0,200*PP_RealTolerance)) THEN
                PRINT '(A,i6,A,F12.6,A)', "Test of VAXPB_OUT FAILED for element:",i,". Value is:",VecOut(i),", should be: 13.0."
                ERROR STOP
            END IF
        END DO
        
        ! VAXPB_IN_OUT
        VecOut = 4.d0
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_VecOut, VecOut, SIZE(VecOut))
#endif
        CALL VAXPB_IN_OUT(nTotal,VecOut,VecIn,d_VecOut,d_VecIn,Const1,Const2)
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyFromDevice(VecOut, d_VecOut, SIZE(VecOut))
#endif
        DO i = 1,nTotal
            IF (.NOT. ALMOSTEQUALABSORREL(VecOut(i),23.d0,200*PP_RealTolerance)) THEN
                PRINT '(A,i6,A,F12.6,A)', "Test of VAXPB_IN_OUT FAILED for element:",i,". Value is:",VecOut(i),", should be: 23.0."
                ERROR STOP
            END IF
        END DO

        ! VAXPB_ADD
        VecOut = 4.d0
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_VecOut, VecOut, SIZE(VecOut))
#endif
        CALL VAXPB_ADD(nTotal,VecOut,VecIn,d_VecOut,d_VecIn)
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyFromDevice(VecOut, d_VecOut, SIZE(VecOut))
#endif
        DO i = 1,nTotal
            IF (.NOT. ALMOSTEQUALABSORREL(VecOut(i),9.d0,200*PP_RealTolerance)) THEN
                PRINT '(A,i6,A,F12.6,A)', "Test of VAXPB_ADD FAILED for element:",i,". Value is:",VecOut(i),", should be: 9.0."
                ERROR STOP
            END IF
        END DO

        ! VAXPB_CONST
        VecOut = 4.d0
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_VecOut, VecOut, SIZE(VecOut))
#endif
        CALL VAXPB_CONST(nTotal,VecOut,d_VecOut,Const1)
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyFromDevice(VecOut, d_VecOut, SIZE(VecOut))
#endif
        DO i = 1,nTotal
            IF (.NOT. ALMOSTEQUALABSORREL(VecOut(i),8.d0,200*PP_RealTolerance)) THEN
                PRINT '(A,i6,A,F12.6,A)', "Test of VAXPB_CONST FAILED for element:",i,". Value is:",VecOut(i),", should be: 8.0."
                ERROR STOP
            END IF
        END DO

    END SUBROUTINE TestMethods

END MODULE TestVAX

!==================================================================================================================================
!> @brief Test driver program
!==================================================================================================================================
PROGRAM VAX_TestDriver
USE TestVAX
IMPLICIT NONE

    CALL InitTest
    CALL TestMethods

END PROGRAM VAX_TestDriver
