#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Unit test 'TestRiemann'
!> Test the routines from module: 'MOD_Riemann'.
!> This file holds the host side code. Device side code in Riemann.cu
!> @author Spencer Starr - 08.2024
!==================================================================================================================================

MODULE TestRiemann
USE ISO_C_BINDING
USE MOD_Riemann
USE MOD_PreProc
USE MOD_EOS_Vars
#ifdef SPLIT_DG
USE MOD_SplitFlux, ONLY: InitSplitDG
#endif
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_Device
USE MOD_DeviceMem
USE MOD_DeviceManage
USE MOD_DG_Vars, ONLY: d_U, d_UPrim, d_U_master, d_UPrim_master
USE MOD_Mesh_Vars, ONLY: d_Metrics_fTilde, d_Metrics_gTilde, d_Metrics_hTilde, d_Ja_Face
#endif
IMPLICIT NONE

! Interfaces to C++ methods needed for tests on device side
#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
    SUBROUTINE TestMethods_Device(riemannID, kappa, d_U_L, d_U_R, d_UPrim_L, d_UPrim_R, d_nv, d_t1, d_t2, d_Flux) &
                                                            BIND(C, NAME="TestMethods_Device")
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(C_INT), VALUE :: riemannID
        REAL, VALUE :: kappa
        INTEGER(C_INT),VALUE :: d_U_L
        INTEGER(C_INT),VALUE :: d_U_R
        INTEGER(C_INT),VALUE :: d_UPrim_L
        INTEGER(C_INT),VALUE :: d_UPrim_R
        INTEGER(C_INT),VALUE :: d_nv, d_t1, d_t2
        INTEGER(C_INT),VALUE :: d_Flux
    END SUBROUTINE TestMethods_Device
END INTERFACE
#endif

! MODULE SCOPE VARIABLES
INTEGER :: nSolvers
INTEGER, ALLOCATABLE :: solver_codes(:)


REAL :: ans(PP_nVar)

REAL :: Flux(PP_nVar)
REAL :: U_L(PP_nVar)
REAL :: U_R(PP_nVar)
REAL :: UPrim_L(PP_nVarPrim)
REAL :: UPrim_R(PP_nVarPrim)
REAL :: nv(3), t1(3), t2(3)

#if (USE_ACCEL != ACCEL_OFF)
INTEGER :: dev_err
INTEGER(C_INT) :: d_Flux 
INTEGER(C_INT) :: d_U_L 
INTEGER(C_INT) :: d_U_R 
INTEGER(C_INT) :: d_UPrim_L 
INTEGER(C_INT) :: d_UPRim_R 
INTEGER(C_INT) :: d_nv 
INTEGER(C_INT) :: d_t1 
INTEGER(C_INT) :: d_t2 
#endif

CONTAINS

    !==================================================================================================================================
    !> @brief Initialize the domain for the test
    !==================================================================================================================================
    SUBROUTINE InitTest

        ! Allocate and init host side arrays
        U_L = [0.1d0, 2.d0, 3.d0, 4.d0, 5.d0]
        U_R = [0.1d0, 2.d0, 3.d0, 4.d0, 5.d0]
        UPrim_L = [0.1d0, 7.d0, 8.d0, 9.d0, 10.d0, 11.d0]
        UPrim_R = [0.1d0, 7.d0, 8.d0, 9.d0, 10.d0, 11.d0]
        nv = [1.d0, 0.d0, 0.d0]
        t1 = [0.d0, 1.d0, 0.d0]
        t2 = [0.d0, 0.d0, 1.d0]

        ! Initialize test solution
        Flux = 0.d0

        ! Set up array of solver codes
#ifndef SPLIT_DG
        nSolvers = 8
        ALLOCATE( solver_codes(nSolvers) )
        solver_codes = [1,2,3,32,33,4,5,6]
#else
        nSolvers = 6
        ALLOCATE( solver_codes(nSolvers) )
        solver_codes = [1,3,32,33,7,0]
#endif

        ! Use the standard split flux methods throughout
#ifdef SPLIT_DG
        CALL InitSplitDG(PRM_SPLITDG_SD)
#endif

        ! Initialize EOS Vars
        Kappa = 1.4d0
        KappaM1  = Kappa-1.d0
        
        ! If device test, allocate and copy on device
#if (USE_ACCEL != ACCEL_OFF)
        CALL AllocateDeviceMemory(d_U_L, SIZE_C_DOUBLE, SIZE(U_L))
        CALL CopyToDevice(d_U_L, U_L, SIZE(U_L))

        CALL AllocateDeviceMemory(d_U_R, SIZE_C_DOUBLE, SIZE(U_R))
        CALL CopyToDevice(d_U_R, U_R, SIZE(U_R))

        CALL AllocateDeviceMemory(d_UPrim_L, SIZE_C_DOUBLE, SIZE(UPrim_L))
        CALL CopyToDevice(d_UPrim_L, UPrim_L, SIZE(UPrim_L))

        CALL AllocateDeviceMemory(d_UPrim_R, SIZE_C_DOUBLE, SIZE(UPrim_R))
        CALL CopyToDevice(d_UPrim_R, UPrim_R, SIZE(UPrim_R))

        CALL AllocateDeviceMemory(d_nv, SIZE_C_DOUBLE, SIZE(nv))
        CALL CopyToDevice(d_nv, nv, SIZE(nv))

        CALL AllocateDeviceMemory(d_t1, SIZE_C_DOUBLE, SIZE(t1))
        CALL CopyToDevice(d_t1, t1, SIZE(t1))

        CALL AllocateDeviceMemory(d_t2, SIZE_C_DOUBLE, SIZE(t2))
        CALL CopyToDevice(d_t2, t2, SIZE(t2))

        CALL AllocateDeviceMemory(d_Flux, SIZE_C_DOUBLE, SIZE(Flux))
        CALL CopyToDevice(d_Flux, Flux, SIZE(Flux))
#endif

    END SUBROUTINE InitTest


    !==================================================================================================================================
    !> @brief Test Riemann solver methods
    !==================================================================================================================================
    SUBROUTINE TestMethods(riemannID)
    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT / OUTPUT VARIABLES
    INTEGER, INTENT(IN) :: riemannID ! ID of the formulation to test

    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER :: i
    !=========================================================================================

        ! Set up the device pointers
        CALL InitRiemann(riemannID, riemannID)

        ! Initialize the Split DG func pointers on the device side as well
#if (USE_ACCEL == ACCEL_OFF)

        ! Invoke the method
        CALL Riemann(Flux, U_L, U_R, UPrim_L, UPrim_R, nv, t1, t2, .FALSE.)

#else
        CALL TestMethods_Device(riemannID, Kappa, d_U_L, d_U_R, d_UPrim_L, d_UPrim_R, d_nv, d_t1, d_t2, d_Flux)
        CALL SynchronizeDevice

        ! Copy the solution from the device
        CALL CopyFromDevice(Flux, d_Flux, SIZE(Flux))
#endif

        ! Check the solution based on which method was selected
        ans = [0.7d0, 14.9d0, 5.6d0, 6.3d0, 105.d0]
        DO i = 1,PP_nVar
            SELECT CASE(riemannID)
            CASE(PRM_RIEMANN_LF)
                IF (.NOT. ALMOSTEQUALABSORREL(Flux(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for Lax-Friedrichs is: ", Flux(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_RIEMANN_ROE)
                IF (.NOT. ALMOSTEQUALABSORREL(Flux(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for Roe is: ", Flux(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_RIEMANN_ROEL2)
                IF (.NOT. ALMOSTEQUALABSORREL(Flux(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for RoeL2 is: ", Flux(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_RIEMANN_ROEENTROPYFIX)
                IF (.NOT. ALMOSTEQUALABSORREL(Flux(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for Roe Entropy Fix is: ", Flux(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
#ifndef SPLIT_DG
            CASE(PRM_RIEMANN_HLL)
                IF (.NOT. ALMOSTEQUALABSORREL(Flux(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for HLL is: ", Flux(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_RIEMANN_HLLC)
                IF (.NOT. ALMOSTEQUALABSORREL(Flux(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for HLLC is: ", Flux(i), " should be: ", ans(i)
                    ERROR STOP
                END IF

            CASE(PRM_RIEMANN_HLLE)
                IF (.NOT. ALMOSTEQUALABSORREL(Flux(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for HLLE is: ", Flux(i), " should be: ", ans(i)
                    ERROR STOP
                END IF

            CASE(PRM_RIEMANN_HLLEM)
                IF (.NOT. ALMOSTEQUALABSORREL(Flux(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for HLLEM is: ", Flux(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
#else
            CASE(PRM_RIEMANN_Average)
                IF (.NOT. ALMOSTEQUALABSORREL(Flux(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for Average is: ", Flux(i), " should be: ", ans(i)
                    ERROR STOP
                END IF

            CASE(PRM_RIEMANN_CH)
                IF (.NOT. ALMOSTEQUALABSORREL(Flux(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Flux value for Chandrashekar formulation is: ", Flux(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
#endif /* SPLIT_DG */
            END SELECT
        END DO

    END SUBROUTINE TestMethods    

END MODULE TestRiemann



!==================================================================================================================================
!> @brief Test driver program
!==================================================================================================================================
PROGRAM SplitFlux_TestDriver
USE TestRiemann
IMPLICIT NONE

INTEGER :: i

    CALL InitTest
    DO i = 1,nSolvers
        CALL TestMethods(solver_codes(i))
    END DO

END PROGRAM SplitFlux_TestDriver
