#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Unit test 'EvalFlux'
!> Test the routines from module: 'MOD_Flux'.
!> This file holds the host side code. Device side code in EvalFlux.cu
!> @author Spencer Starr - 08.2024
!==================================================================================================================================

MODULE TestEvalFlux
USE ISO_C_BINDING
USE MOD_DG_Vars, ONLY: U, UPrim, F, G, H
USE MOD_Mesh_Vars, ONLY: Metrics_fTilde, Metrics_gTilde, Metrics_hTilde
USE MOD_Mesh_Vars, ONLY: nElems
USE MOD_Flux
USE MOD_PreProc
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: gradUx,gradUy,gradUz
#endif
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_Device
USE MOD_DeviceMem
USE MOD_DG_Vars, ONLY: d_U, d_UPrim, d_F, d_G, d_H
USE MOD_Mesh_Vars, ONLY: d_Metrics_fTilde, d_Metrics_gTilde, d_Metrics_hTilde
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: d_gradUx,d_gradUy,d_gradUz
#endif
#endif
IMPLICIT NONE

! Interfaces to CUDA C++ methods needed for tests on device side (EvalFlux.cu)
#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
    SUBROUTINE TestEvalEulerFlux1D_Device(d_U_Long, d_Flux) BIND(C, NAME="TestEvalEulerFlux1D_Device")
        USE ISO_C_BINDING, ONLY: C_INT
        IMPLICIT NONE
        INTEGER(C_INT),VALUE :: d_U_Long
        INTEGER(C_INT),VALUE :: d_Flux
    END SUBROUTINE
END INTERFACE
#endif

! MODULE SCOPE VARIABLES
REAL :: Flux(PP_nVar)
REAL :: U_Long(PP_2Var)

#if (USE_ACCEL != ACCEL_OFF)
INTEGER :: dev_err
INTEGER(C_INT) :: d_Flux 
INTEGER(C_INT) :: d_U_Long 
#endif

CONTAINS

    !==================================================================================================================================
    !> @brief Initialize the domain for the test
    !==================================================================================================================================
    SUBROUTINE InitTest

        nElems = 16

        ! Allocate and init host side arrays - use 16 5th order elements
        ALLOCATE( U(PP_nVar,0:1,0:1,0:1,1:nElems) )
        ALLOCATE( UPrim(PP_nVarPrim,0:1,0:1,0:1,1:nElems) )
        ALLOCATE( F(PP_nVar,0:1,0:1,0:1,1:nElems) )
        ALLOCATE( G(PP_nVar,0:1,0:1,0:1,1:nElems) )
        ALLOCATE( H(PP_nVar,0:1,0:1,0:1,1:nElems) )
        ALLOCATE( Metrics_fTilde(1:3,0:1,0:1,0:1,nElems,0:FV_SIZE) )
        ALLOCATE( Metrics_gTilde(1:3,0:1,0:1,0:1,nElems,0:FV_SIZE) )
        ALLOCATE( Metrics_hTilde(1:3,0:1,0:1,0:1,nElems,0:FV_SIZE) )
#if PARABOLIC
        ALLOCATE( gradUx(PP_nVarLifting,0:1,0:1,0:1,nElems) )
        ALLOCATE( gradUy(PP_nVarLifting,0:1,0:1,0:1,nElems) )
        ALLOCATE( gradUz(PP_nVarLifting,0:1,0:1,0:1,nElems) )
        gradUx=0.0
        gradUy=0.0
        gradUz=0.0
#endif

        U(DENS,:,:,:,:) = 2.0
        U(MOM1,:,:,:,:) = 3.0
        U(MOM2,:,:,:,:) = 4.0
        U(MOM3,:,:,:,:) = 5.0
        U(ENER,:,:,:,:) = 6.0

        UPrim(DENS,:,:,:,:) = 7.0
        UPrim(VEL1,:,:,:,:) = 8.0
        UPrim(VEL2,:,:,:,:) = 9.0
        UPrim(VEL3,:,:,:,:) = 10.0
        UPrim(PRES,:,:,:,:) = 11.0
        UPrim(TEMP,:,:,:,:) = 12.0

        Metrics_fTilde(1,:,:,:,:,:) = 13.0
        Metrics_fTilde(2,:,:,:,:,:) = 14.0
        Metrics_fTilde(3,:,:,:,:,:) = 15.0

        Metrics_gTilde(1,:,:,:,:,:) = 13.0
        Metrics_gTilde(2,:,:,:,:,:) = 14.0
        Metrics_gTilde(3,:,:,:,:,:) = 15.0

        Metrics_hTilde(1,:,:,:,:,:) = 13.0
        Metrics_hTilde(2,:,:,:,:,:) = 14.0
        Metrics_hTilde(3,:,:,:,:,:) = 15.0

        U_Long = [2.d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0, 11.d0, 12.d0]

        ! Initialize test solutions
        F = 0.0
        G = 0.0
        H = 0.0
        Flux = 0.0
        
        ! If device test, allocate and copy on device
#if (USE_ACCEL != ACCEL_OFF)
        CALL AllocateDeviceMemory(d_U, SIZE_C_DOUBLE, SIZE(U))
        CALL CopyToDevice(d_U, U, SIZE(U))

        CALL AllocateDeviceMemory(d_UPrim, SIZE_C_DOUBLE, SIZE(UPrim))
        CALL CopyToDevice(d_UPrim, UPrim, SIZE(UPrim))

        CALL AllocateDeviceMemory(d_F, SIZE_C_DOUBLE, SIZE(F))
        CALL CopyToDevice(d_F, F, SIZE(F))

        CALL AllocateDeviceMemory(d_G, SIZE_C_DOUBLE, SIZE(G))
        CALL CopyToDevice(d_G, G, SIZE(G))

        CALL AllocateDeviceMemory(d_H, SIZE_C_DOUBLE, SIZE(H))
        CALL CopyToDevice(d_H, H, SIZE(H))

        CALL AllocateDeviceMemory(d_Metrics_fTilde, SIZE_C_DOUBLE, SIZE(Metrics_fTilde))
        CALL CopyToDevice(d_Metrics_fTilde, Metrics_fTilde, SIZE(Metrics_fTilde))

        CALL AllocateDeviceMemory(d_Metrics_gTilde, SIZE_C_DOUBLE, SIZE(Metrics_gTilde))
        CALL CopyToDevice(d_Metrics_gTilde, Metrics_gTilde, SIZE(Metrics_gTilde))

        CALL AllocateDeviceMemory(d_Metrics_hTilde, SIZE_C_DOUBLE, SIZE(Metrics_hTilde))
        CALL CopyToDevice(d_Metrics_hTilde, Metrics_hTilde, SIZE(Metrics_hTilde))

        CALL AllocateDeviceMemory(d_Flux, SIZE_C_DOUBLE, SIZE(Flux))
        CALL CopyToDevice(d_Flux, Flux, SIZE(Flux))

        CALL AllocateDeviceMemory(d_U_Long, SIZE_C_DOUBLE, SIZE(U_Long))
        CALL CopyToDevice(d_U_Long, U_Long, SIZE(U_Long))

#if PARABOLIC
        CALL AllocateDeviceMemory(d_gradUx, SIZE_C_DOUBLE, SIZE(gradUx))
        CALL CopyToDevice(d_gradUx, gradUx, SIZE(gradUx))

        CALL AllocateDeviceMemory(d_gradUy, SIZE_C_DOUBLE, SIZE(gradUy))
        CALL CopyToDevice(d_gradUy, gradUy, SIZE(gradUy))

        CALL AllocateDeviceMemory(d_gradUz, SIZE_C_DOUBLE, SIZE(gradUz))
        CALL CopyToDevice(d_gradUz, gradUz, SIZE(gradUz))
#endif
#endif

    END SUBROUTINE InitTest


    !==================================================================================================================================
    !> @brief Test SplitFlux volume methods
    !==================================================================================================================================
    SUBROUTINE TestEvalTransformedFlux
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER :: iElem,i,j,k,n
    REAL :: ans(PP_nVar)
    !=========================================================================================

        ! Invoke the method
        CALL EvalTransformedFlux3D(1, 16 &
                                    ,U(    :,:,:,:,1:16) &
                                    ,UPrim(:,:,:,:,1:16) &
#if PARABOLIC
                                    ,gradUx(:,:,:,:,1:16) &
                                    ,gradUy(:,:,:,:,1:16) &
                                    ,gradUz(:,:,:,:,1:16) &
#endif
                                    ,F( :,:,:,:,1:16)  &
                                    ,G( :,:,:,:,1:16)  &
                                    ,H( :,:,:,:,1:16)  &
                                    ,Metrics_fTilde(:,:,:,:,1:16,0) &
                                    ,Metrics_gTilde(:,:,:,:,1:16,0) &
                                    ,Metrics_hTilde(:,:,:,:,1:16,0) &
                                    ,0)

#if (USE_ACCEL != ACCEL_OFF)
        ! Copy the solution from the device
        CALL CopyFromDevice(F, d_F, SIZE(F))
        CALL CopyFromDevice(G, d_G, SIZE(G))
        CALL CopyFromDevice(H, d_H, SIZE(H))
#endif

        ! Check the solution based on which method was selected
#if PARABOLIC
        ans = [0.d0, 0.d0, 0.d0, 0.d0, 0.d0]
#else
        ans = [170.d0, 1503.d0, 1684.d0, 1865.d0, 1445.d0]
#endif
        DO iElem = 1,16
            DO k = 0,1; DO j = 0,1; DO i = 0,1
                DO n = 1,PP_nVar
                    IF (.NOT. ALMOSTEQUALABSORREL(F(n,i,j,k,iElem),ans(n),200*PP_RealTolerance)) THEN
                        PRINT *, "Value for F from EvalTransformedFlux3D for variable: ", n, " in DOF: ",i,j,k, &
                                 " of element: ", iElem, " is: ", F(n,i,j,k,iElem), " should be: ", ans(n)
                        ERROR STOP
                    END IF

                    IF (.NOT. ALMOSTEQUALABSORREL(G(n,i,j,k,iElem),ans(n),200*PP_RealTolerance)) THEN
                        PRINT *, "Value for G from EvalTransformedFlux3D for variable: ", n, " in DOF: ",i,j,k, &
                                 " of element: ", iElem, " is: ", G(n,i,j,k,iElem), " should be: ", ans(n)
                        ERROR STOP
                    END IF

                    IF (.NOT. ALMOSTEQUALABSORREL(H(n,i,j,k,iElem),ans(n),200*PP_RealTolerance)) THEN
                        PRINT *, "Value for H from EvalTransformedFlux3D for variable: ", n, " in DOF: ",i,j,k, &
                                 " of element: ", iElem, " is: ", H(n,i,j,k,iElem), " should be: ", ans(n)
                        ERROR STOP
                    END IF
                END DO
            END DO; END DO; END DO
        END DO

    END SUBROUTINE TestEvalTransformedFlux


    !==================================================================================================================================
    !> @brief Test host and device backends for EvalEulerFlux1D_fast method
    !==================================================================================================================================
    SUBROUTINE TestEvalEulerFlux
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER :: i
    REAL :: ans(PP_nVar)
    !=========================================================================================

        ! Initialize the Split DG func pointers on the device side as well
#if (USE_ACCEL == ACCEL_OFF)

        ! Invoke the method
        CALL EvalEulerFlux1D_fast(U_Long, Flux)

#else
        CALL TestEvalEulerFlux1D_Device(d_U_Long, d_Flux)

        ! Copy the solution from the device
        CALL CopyFromDevice(Flux, d_Flux, SIZE(Flux))
#endif

        ! Check the solution based on which method was selected
        ans = [3.d0, 35.d0, 27.d0, 30.d0, 136.d0]
        DO i = 1,PP_nVar
            IF (.NOT. ALMOSTEQUALABSORREL(Flux(i),ans(i),200*PP_RealTolerance)) THEN
                PRINT *, "EvalEulerFlux1D value for variable: ", i, " is: ", Flux(i), " should be: ", ans(i)
                ERROR STOP
            END IF
        END DO


    END SUBROUTINE TestEvalEulerFlux
    

END MODULE TestEvalFlux



!==================================================================================================================================
!> @brief Test driver program
!==================================================================================================================================
PROGRAM EvalFlux_TestDriver
USE TestEvalFlux
IMPLICIT NONE

    CALL InitTest
    CALL TestEvalTransformedFlux
    CALL TestEvalEulerFlux

END PROGRAM EvalFlux_TestDriver
