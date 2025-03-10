#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Unit test 'SplitFlux'
!> Test the routines from module: 'MOD_SplitFlux'.
!> This file holds the host side code. Device side code in 
!> @author Spencer Starr - 07.2024
!==================================================================================================================================

MODULE TestSplitFlux
USE ISO_C_BINDING
USE MOD_DG_Vars, ONLY: U, UPrim, U_master, UPrim_master
USE MOD_Mesh_Vars, ONLY: Metrics_fTilde, Metrics_gTilde, Metrics_hTilde, Ja_Face
USE MOD_Mesh_Vars, ONLY: nElems, nSides
USE MOD_SplitFlux
USE MOD_PreProc
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
    SUBROUTINE TestVolumeMethods_Device(splitDG_ID, d_U, d_UPrim, d_U_master, d_UPrim_master, d_Metrics, d_Ja_Face, d_Flux) &
                                                            BIND(C, NAME="TestVolumeMethods_Device")
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(C_INT), VALUE :: splitDG_ID
        INTEGER(C_INT),VALUE :: d_U
        INTEGER(C_INT),VALUE :: d_UPrim
        INTEGER(C_INT),VALUE :: d_U_master
        INTEGER(C_INT),VALUE :: d_UPrim_master
        INTEGER(C_INT),VALUE :: d_Metrics
        INTEGER(C_INT),VALUE :: d_Ja_Face
        INTEGER(C_INT),VALUE :: d_Flux
    END SUBROUTINE
END INTERFACE

INTERFACE
    SUBROUTINE TestSurfaceMethods_Device(splitDG_ID, d_U_LL, d_U_RR, d_Flux) BIND(C, NAME="TestSurfaceMethods_Device")
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(C_INT), VALUE :: splitDG_ID
        INTEGER(C_INT),VALUE :: d_U_LL
        INTEGER(C_INT),VALUE :: d_U_RR
        INTEGER(C_INT),VALUE :: d_Flux
    END SUBROUTINE
END INTERFACE
#endif

! MODULE SCOPE VARIABLES
REAL :: Flux_Vol(PP_nVar, 0:1)
REAL :: ans(PP_nVar)

REAL :: Flux_Surf(PP_nVar)
REAL :: U_LL(PP_2Var)
REAL :: U_RR(PP_2Var)

#if (USE_ACCEL != ACCEL_OFF)
INTEGER :: dev_err
INTEGER(C_INT) :: d_Flux_Vol 
INTEGER(C_INT) :: d_Flux_Surf 
INTEGER(C_INT) :: d_U_LL 
INTEGER(C_INT) :: d_U_RR 
#endif

CONTAINS

    !==================================================================================================================================
    !> @brief Initialize the domain for the test
    !==================================================================================================================================
    SUBROUTINE InitTest

        nElems = 1
        nSides = 1

        ! Allocate and init host side arrays - use 16 5th order elements
        ALLOCATE( U(PP_nVar,0:1,0:1,0:1,1:nElems) )
        ALLOCATE( UPrim(PP_nVarPrim,0:1,0:1,0:1,1:nElems) )
        ALLOCATE( U_master(PP_nVar,0:1,0:1,1:nSides) )
        ALLOCATE( UPrim_master(PP_nVarPrim,0:1,0:1,1:nSides) )
        ALLOCATE( Metrics_fTilde(1:3,0:1,0:1,0:1,nElems,0:FV_SIZE) )
        ALLOCATE( Ja_Face(1:3,1:3,0:1,0:1,1:nSides) )

        U(DENS,:,:,:,:) = 2.0
        U(MOM1,:,:,:,:) = 3.0
        U(MOM2,:,:,:,:) = 4.0
        U(MOM3,:,:,:,:) = 5.0
        U(ENER,:,:,:,:) = 6.0

        U_master(DENS,:,:,:) = 2.0
        U_master(MOM1,:,:,:) = 3.0
        U_master(MOM2,:,:,:) = 4.0
        U_master(MOM3,:,:,:) = 5.0
        U_master(ENER,:,:,:) = 6.0

        UPrim(DENS,:,:,:,:) = 7.0
        UPrim(VEL1,:,:,:,:) = 8.0
        UPrim(VEL2,:,:,:,:) = 9.0
        UPrim(VEL3,:,:,:,:) = 10.0
        UPrim(PRES,:,:,:,:) = 11.0
        UPrim(TEMP,:,:,:,:) = 12.0

        UPrim_master(DENS,:,:,:) = 7.0
        UPrim_master(VEL1,:,:,:) = 8.0
        UPrim_master(VEL2,:,:,:) = 9.0
        UPrim_master(VEL3,:,:,:) = 10.0
        UPrim_master(PRES,:,:,:) = 11.0
        UPrim_master(TEMP,:,:,:) = 12.0

        Metrics_fTilde(1,:,:,:,:,:) = 13.0
        Metrics_fTilde(2,:,:,:,:,:) = 14.0
        Metrics_fTilde(3,:,:,:,:,:) = 15.0

        Ja_Face(:,1,:,:,:) = 16.0
        Ja_Face(:,2,:,:,:) = 17.0
        Ja_Face(:,3,:,:,:) = 18.0

        U_LL = [2.d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0, 11.d0, 12.d0]
        U_RR = [2.d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0, 11.d0, 12.d0]

        ! Initialize test solution
        Flux_Vol = 0.0
        Flux_Surf = 0.0
        
        ! If device test, allocate and copy on device
#if (USE_ACCEL != ACCEL_OFF)
        CALL AllocateDeviceMemory(d_U, SIZE_C_DOUBLE, SIZE(U))
        CALL CopyToDevice(d_U, U, SIZE(U))

        CALL AllocateDeviceMemory(d_UPrim, SIZE_C_DOUBLE, SIZE(UPrim))
        CALL CopyToDevice(d_UPrim, UPrim, SIZE(UPrim))

        CALL AllocateDeviceMemory(d_U_master, SIZE_C_DOUBLE, SIZE(U_master))
        CALL CopyToDevice(d_U_master, U_master, SIZE(U_master))

        CALL AllocateDeviceMemory(d_UPrim_master, SIZE_C_DOUBLE, SIZE(UPrim_master))
        CALL CopyToDevice(d_UPrim_master, UPrim_master, SIZE(UPrim_master))

        CALL AllocateDeviceMemory(d_Metrics_fTilde, SIZE_C_DOUBLE, SIZE(Metrics_fTilde))
        CALL CopyToDevice(d_Metrics_fTilde, Metrics_fTilde, SIZE(Metrics_fTilde))

        CALL AllocateDeviceMemory(d_Ja_Face, SIZE_C_DOUBLE, SIZE(Ja_Face))
        CALL CopyToDevice(d_Ja_Face, Ja_Face, SIZE(Ja_Face))

        CALL AllocateDeviceMemory(d_Flux_Vol, SIZE_C_DOUBLE, SIZE(Flux_Vol))
        CALL CopyToDevice(d_Flux_Vol, Flux_Vol, SIZE(Flux_Vol))

        CALL AllocateDeviceMemory(d_Flux_Surf, SIZE_C_DOUBLE, SIZE(Flux_Surf))
        CALL CopyToDevice(d_Flux_Surf, Flux_Surf, SIZE(Flux_Surf))

        CALL AllocateDeviceMemory(d_U_LL, SIZE_C_DOUBLE, SIZE(U_LL))
        CALL CopyToDevice(d_U_LL, U_LL, SIZE(U_LL))

        CALL AllocateDeviceMemory(d_U_RR, SIZE_C_DOUBLE, SIZE(U_RR))
        CALL CopyToDevice(d_U_RR, U_RR, SIZE(U_RR))
#endif

    END SUBROUTINE InitTest


    !==================================================================================================================================
    !> @brief Test SplitFlux volume methods
    !==================================================================================================================================
    SUBROUTINE TestVolumeMethods(splitDG_ID)
    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT / OUTPUT VARIABLES
    INTEGER, INTENT(IN) :: splitDG_ID ! ID of the formulation to test

    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER :: i
    !=========================================================================================

        ! Set up the device pointers
        CALL InitSplitDG(splitDG_ID)

        ! Initialize the Split DG func pointers on the device side as well
#if (USE_ACCEL == ACCEL_OFF)

        ! Invoke the method
        CALL SplitDGVolume_pointer(U(:,0,0,0,1) ,UPrim(:,0,0,0,1), U_master(:,0,0,1),UPrim_master(:,0,0,1),   &
                                   Metrics_fTilde(:,0,0,0,1,0), Ja_face(1,:,0,0,1), Flux_Vol(:,0))

#else
        CALL TestVolumeMethods_Device(splitDG_ID, d_U, d_UPrim, d_U_master, d_UPrim_master, d_Metrics_fTilde, d_Ja_Face, d_Flux_Vol)

        ! Copy the solution from the device
        CALL CopyFromDevice(Flux_Vol, d_Flux_Vol, SIZE(Flux_Vol))
#endif

        ! Check the solution based on which method was selected
        DO i = 1,PP_nVar
            SELECT CASE(splitDG_ID)
            CASE(PRM_SPLITDG_SD)
                ans = [376.d0, 2842.d0, 3560.d0, 4123.d0, 14297.d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Vol(i,0),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for standard flux splitting formulation for element: ", i, " is: ", Flux_Vol(i,0), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_SPLITDG_MO)
                ans = [376.d0, 3327.d0, 3725.d0, 4123.d0, -145688.d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Vol(i,0),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for Morinishi formulation for element: ", i, " is: ", Flux_Vol(i,0), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_SPLITDG_DU)
                ans = [1682.d0, 2842.d0, 3705.d0, 4568.d0, 14297.d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Vol(i,0),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for Ducros formulation for element: ", i, " is: ", Flux_Vol(i,0), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_SPLITDG_KG)
                ans = [1682.d0, 13775.d0, 15479.d0, 17183.d0, 14297.d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Vol(i,0),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for Kennedy and Gruber formulation for element: ", i, " is: ", Flux_Vol(i,0), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_SPLITDG_PI)
                ans = [1682.d0, 13775.d0, 15479.d0, 17183.d0, 14297.d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Vol(i,0),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for Pirozzoli formulation for element: ", i, " is: ", Flux_Vol(i,0), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_SPLITDG_CH)
                ans = [1682.d0, 13775.d0, 15479.d0, 17183.d0, 215296.d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Vol(i,0),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Volume flux value for Chandrashekar formulation for element: ", i, " is: ", Flux_Vol(i,0), " should be: ", ans(i)
                    ERROR STOP
                END IF
            END SELECT
        END DO

    END SUBROUTINE TestVolumeMethods


    !==================================================================================================================================
    !> @brief Test SplitFlux volume methods
    !==================================================================================================================================
    SUBROUTINE TestSurfaceMethods(splitDG_ID)
    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT / OUTPUT VARIABLES
    INTEGER, INTENT(IN) :: splitDG_ID ! ID of the formulation to test

    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER :: i
    !=========================================================================================
        
        ! Set up the device pointers
        CALL InitSplitDG(splitDG_ID)

        ! Initialize the Split DG func pointers on the device side as well
#if (USE_ACCEL == ACCEL_OFF)

        ! Invoke the method
        CALL SplitDGSurface_pointer(U_LL, U_RR, Flux_Surf)

#else
        CALL TestSurfaceMethods_Device(splitDG_ID, d_U_LL, d_U_RR, d_Flux_Surf)

        ! Copy the solution from the device
        CALL CopyFromDevice(Flux_Surf, d_Flux_Surf, SIZE(Flux_Surf))
#endif

        ! Check the solution based on which method was selected
        DO i = 1,PP_nVar
            SELECT CASE(splitDG_ID)
            CASE(PRM_SPLITDG_SD)
                ans = [3.d0, 35.d0, 27.d0, 30.d0, 136.d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Surf(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Surface flux value for standard flux splitting formulation for element: ", i, " is: ", Flux_Surf(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_SPLITDG_MO)
                ans = [3.d0, 35.d0, 27.d0, 30.d0, -1456.5d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Surf(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Surface flux value for Morinishi formulation for element: ", i, " is: ", Flux_Surf(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_SPLITDG_DU)
                ans = [16.d0, 35.d0, 32.d0, 40.d0, 136.d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Surf(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Surface flux value for Ducros formulation for element: ", i, " is: ", Flux_Surf(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_SPLITDG_KG)
                ans = [16.d0, 139.d0, 144.d0, 160.d0, 136.d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Surf(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Surface flux value for Kennedy and Gruber formulation for element: ", i, " is: ", Flux_Surf(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_SPLITDG_PI)
                ans = [16.d0, 139.d0, 144.d0, 160.d0, 136.d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Surf(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Surface flux value for Pirozzoli formulation for element: ", i, " is: ", Flux_Surf(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
            CASE(PRM_SPLITDG_CH)
                ans = [16.d0, 139.d0, 144.d0, 160.d0, 2048.d0]
                IF (.NOT. ALMOSTEQUALABSORREL(Flux_Surf(i),ans(i),200*PP_RealTolerance)) THEN
                    PRINT *, "Surface flux value for Chandrashekar formulation for element: ", i, " is: ", Flux_Surf(i), " should be: ", ans(i)
                    ERROR STOP
                END IF
            END SELECT
        END DO


    END SUBROUTINE TestSurfaceMethods
    

END MODULE TestSplitFlux



!==================================================================================================================================
!> @brief Test driver program
!==================================================================================================================================
PROGRAM SplitFlux_TestDriver
USE TestSplitFlux
IMPLICIT NONE

INTEGER :: i

    CALL InitTest
    DO i = PRM_SPLITDG_SD, PRM_SPLITDG_CH
        CALL TestVolumeMethods(i)
        CALL TestSurfaceMethods(i)
    END DO

END PROGRAM SplitFlux_TestDriver