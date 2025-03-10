#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Unit test 'TestEOS'
!> Test the routines from module: 'EOS'.
!> @author Spencer Starr - 08.2024
!==================================================================================================================================

MODULE TestEOS
USE ISO_C_BINDING
USE MOD_EOS
USE MOD_PreProc
USE MOD_DG_Vars, ONLY: U, UPrim, d_U, d_UPrim
#if PP_EntropyVars == 1
USE MOD_DG_Vars, ONLY: V, d_V
#endif
USE MOD_Mesh_Vars, ONLY: nElems
USE MOD_EOS_Vars,ONLY: Kappa, sKappaM1, KappaM1, R
USE MOD_Device
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_DeviceMem
USE MOD_DeviceManage
#endif /* USE_ACCEL != ACCEL_OFF */ 
IMPLICIT NONE

CONTAINS

    !==================================================================================================================================
    !> @brief Initialize data for the test
    !==================================================================================================================================
    SUBROUTINE InitTest
    !=========================================================================================

        ! Allocate host arrays as if using an 4x4x4 element grid with 5th order elements
        nElems = 16
        ALLOCATE( U(1:PP_nVar,0:5,0:5,0:5,1:nElems) )
        ALLOCATE( UPrim(1:PP_nVarPrim,0:5,0:5,0:5,1:nElems) )
#if PP_EntropyVars == 1
        ALLOCATE( V(1:PP_nVar,0:5,0:5,0:5,1:nElems) )
#endif

        ! Set them with simple initial data with which a response from the kernel can be measured
        U(DENS,:,:,:,:) = 2.d0
        U(MOM1,:,:,:,:) = 3.d0
        U(MOM2,:,:,:,:) = 4.d0
        U(MOM3,:,:,:,:) = 5.d0
        U(ENER,:,:,:,:) = 15.d0
        UPrim(:,:,:,:,:) = 100.d0
        
        ! Initialize everything else
        Kappa = 1.d0
        KappaM1 = 0.5d0
        sKappaM1 = 1.d0/KappaM1
        R = 10.d0

        ! If testing accelerated methods, get flow vars onto the device
#if (USE_ACCEL != ACCEL_OFF)
        CALL AllocateDeviceMemory(d_U, SIZE_C_DOUBLE, size(U))
        CALL CopyToDevice(d_U, U, SIZE(U))

        CALL AllocateDeviceMemory(d_UPrim, SIZE_C_DOUBLE, size(UPrim))
        CALL CopyToDevice(d_Uprim, UPrim, SIZE(UPrim))
#if PP_EntropyVars == 1
        CALL AllocateDeviceMemory(d_V, SIZE_C_DOUBLE, size(V))
        CALL CopyToDevice(d_V, V, SIZE(V))
#endif /* PP_EntropyVars */
#endif /* USE_ACCEL */

    END SUBROUTINE InitTest

    !==================================================================================================================================
    !> @brief Test ConsToPrim_Volume method
    !==================================================================================================================================
    SUBROUTINE TestConsToPrim_Volume
    ! LOCAL VARIABLES
    INTEGER :: i,j,k, iElem
    CHARACTER(LEN=200) :: err_msg
    !=========================================================================================

        ! Reset the solution array
        UPrim = 100.d0
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_UPrim, UPrim, SIZE(UPrim))
#endif

        CALL ConsToPrim(5, UPrim, U, d_UPrim, d_U, STREAM_DEFAULT)

        ! If running on the device, copy data back over
#if (USE_ACCEL != ACCEL_OFF)
        CALL SynchronizeDevice()
        CALL CopyFromDevice(UPrim, d_UPrim, SIZE(UPrim))
#endif

        ! Check values are correct
        DO iElem=1,nElems
            DO k=0,5; DO j=0,5; DO i=0,5
                IF (UPrim(DENS,i,j,k,iElem) /= 2.0) THEN
                    write(err_msg, '("Calc of density in ConsToPrim_Volume is: ", f12.6, " for ",4(i3),". Should be: 2.0.")') UPrim(DENS,i,j,k,iElem),i,j,k,iElem
                    ERROR STOP err_msg
                END IF

                IF (UPrim(VEL1,i,j,k,iElem) /= 1.5) THEN
                    write(err_msg, '("Calc of x-velocity in ConsToPrim_Volume is: ", f12.6, " for ",4(i3),". Should be: 1.5.")') UPrim(VEL1,i,j,k,iElem),i,j,k,iElem
                    ERROR STOP err_msg
                END IF
        
                IF (UPrim(VEL2,i,j,k,iElem) /= 2.0) THEN
                    write(err_msg, '("Calc of y-velocity in ConsToPrim_Volume is: ", f12.6, " for ",4(i3),". Should be: 2.0.")') UPrim(VEL2,i,j,k,iElem),i,j,k,iElem
                    ERROR STOP err_msg
                END IF
        
                IF (UPrim(VEL3,i,j,k,iElem) /= 2.5) THEN
                    write(err_msg, '("Calc of z-velocity in ConsToPrim_Volume is: ", f12.6, " for ",4(i3),". Should be: 2.5.")') UPrim(VEL3,i,j,k,iElem),i,j,k,iElem
                    ERROR STOP err_msg
                END IF
        
                IF (UPrim(PRES,i,j,k,iElem) /= 1.25) THEN
                    write(err_msg, '("Calc of pressure in ConsToPrim_Volume is: ", f12.6, " for ",4(i3),". Should be: 1.25.")') UPrim(PRES,i,j,k,iElem),i,j,k,iElem
                    ERROR STOP err_msg
                END IF
        
                IF (UPrim(TEMP,i,j,k,iElem) /= 0.0625) THEN
                    write(err_msg, '("Calc of temperature in ConsToPrim_Volume is: ", f12.6, " for ",4(i3),". Should be: 0.0625.")') UPrim(ENER,i,j,k,iElem),i,j,k,iElem
                    ERROR STOP err_msg
                END IF
            END DO; END DO; END DO! i,j,k=0,Nloc
        END DO ! iElem        

    END SUBROUTINE TestConsToPrim_Volume

    !==================================================================================================================================
    !> @brief Test ConsToPrim_Elem method
    !==================================================================================================================================
    SUBROUTINE TestConsToPrim_Elem
    ! LOCAL VARIABLES
    INTEGER :: i,j,k
    INTEGER :: elemIdx = 12
    CHARACTER(LEN=200) :: err_msg
    !=========================================================================================

        ! Reset the solution array
        UPrim = 100.d0
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_UPrim, UPrim, SIZE(UPrim))
#endif

        CALL ConsToPrim(5, UPrim(:,:,:,:,elemIdx), U(:,:,:,:,elemIdx), d_UPrim, d_U, elemIdx)

        ! If running on the device, copy data back over
#if (USE_ACCEL != ACCEL_OFF)
        CALL SynchronizeDevice()
        CALL CopyFromDevice(UPrim, d_UPrim, SIZE(UPrim))
#endif

        ! Check values are correct
        DO k=0,5; DO j=0,5; DO i=0,5
            IF (UPrim(DENS,i,j,k,elemIdx) /= 2.0) THEN
                write(err_msg, '("Calc of density in ConsToPrim_Elem is: ", f12.6, " for ",3(i3),". Should be: 2.0.")') UPrim(DENS,i,j,k,elemIdx),i,j,k
                ERROR STOP err_msg
            END IF

            IF (UPrim(VEL1,i,j,k,elemIdx) /= 1.5) THEN
                write(err_msg, '("Calc of x-velocity in ConsToPrim_Elem is: ", f12.6, " for ",3(i3),". Should be: 1.5.")') UPrim(VEL1,i,j,k,elemIdx),i,j,k
                ERROR STOP err_msg
            END IF
    
            IF (UPrim(VEL2,i,j,k,elemIdx) /= 2.0) THEN
                write(err_msg, '("Calc of y-velocity in ConsToPrim_Elem is: ", f12.6, " for ",3(i3),". Should be: 2.0.")') UPrim(VEL2,i,j,k,elemIdx),i,j,k
                ERROR STOP err_msg
            END IF
    
            IF (UPrim(VEL3,i,j,k,elemIdx) /= 2.5) THEN
                write(err_msg, '("Calc of z-velocity in ConsToPrim_Elem is: ", f12.6, " for ",3(i3),". Should be: 2.5.")') UPrim(VEL3,i,j,k,elemIdx),i,j,k
                ERROR STOP err_msg
            END IF
    
            IF (UPrim(PRES,i,j,k,elemIdx) /= 1.25) THEN
                write(err_msg, '("Calc of pressure in ConsToPrim_Elem is: ", f12.6, " for ",3(i3),". Should be: 1.25.")') UPrim(PRES,i,j,k,elemIdx),i,j,k
                ERROR STOP err_msg
            END IF
    
            IF (UPrim(TEMP,i,j,k,elemIdx) /= 0.0625) THEN
                write(err_msg, '("Calc of temperature in ConsToPrim_Elem is: ", f12.6, " for ",3(i3),". Should be: 0.0625.")') UPrim(ENER,i,j,k,elemIdx),i,j,k
                ERROR STOP err_msg
            END IF
        END DO; END DO; END DO! i,j,k=0,Nloc

    END SUBROUTINE TestConsToPrim_Elem

    !==================================================================================================================================
    !> @brief Test ConsToPrim_Side method
    !==================================================================================================================================
    SUBROUTINE TestConsToPrim_Side
    ! LOCAL VARIABLES
    INTEGER :: i,j
    INTEGER :: sideIdx = 3
    INTEGER :: elemIdx = 12
    CHARACTER(LEN=200) :: err_msg
    !=========================================================================================

        ! Reset the solution array
        UPrim = 100.d0
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_UPrim, UPrim, SIZE(UPrim))
#endif

        ! Call single side method
        CALL ConsToPrim(5, UPrim(:,:,:,sideIdx,elemIdx), U(:,:,:,sideIdx,elemIdx), d_UPrim, d_U, elemIdx, sideIdx)

        ! If running on the device, copy data back over
#if (USE_ACCEL != ACCEL_OFF)
        CALL SynchronizeDevice()        
        CALL CopyFromDevice(UPrim, d_UPrim, SIZE(UPrim))
#endif

        ! Check values are correct
        DO j=0,5; DO i=0,5
            IF (UPrim(DENS,i,j,sideIdx,elemIdx) /= 2.0) THEN
                write(err_msg, '("Calc of density in ConsToPrim_Side is: ", f12.6, " for ",2(i3),". Should be: 2.0.")') UPrim(DENS,i,j,sideIdx,elemIdx),i,j
                ERROR STOP err_msg
            END IF

            IF (UPrim(VEL1,i,j,sideIdx,elemIdx) /= 1.5) THEN
                write(err_msg, '("Calc of x-velocity in ConsToPrim_Side is: ", f12.6, " for ",2(i3),". Should be: 1.5.")') UPrim(VEL1,i,j,sideIdx,elemIdx),i,j
                ERROR STOP err_msg
            END IF
    
            IF (UPrim(VEL2,i,j,sideIdx,elemIdx) /= 2.0) THEN
                write(err_msg, '("Calc of y-velocity in ConsToPrim_Side is: ", f12.6, " for ",2(i3),". Should be: 2.0.")') UPrim(VEL2,i,j,sideIdx,elemIdx),i,j
                ERROR STOP err_msg
            END IF
    
            IF (UPrim(VEL3,i,j,sideIdx,elemIdx) /= 2.5) THEN
                write(err_msg, '("Calc of z-velocity in ConsToPrim_Side is: ", f12.6, " for ",2(i3),". Should be: 2.5.")') UPrim(VEL3,i,j,sideIdx,elemIdx),i,j
                ERROR STOP err_msg
            END IF
    
            IF (UPrim(PRES,i,j,sideIdx,elemIdx) /= 1.25) THEN
                write(err_msg, '("Calc of pressure in ConsToPrim_Side is: ", f12.6, " for ",2(i3),". Should be: 1.25.")') UPrim(PRES,i,j,sideIdx,elemIdx),i,j
                ERROR STOP err_msg
            END IF
    
            IF (UPrim(TEMP,i,j,sideIdx,elemIdx) /= 0.0625) THEN
                write(err_msg, '("Calc of temperature in ConsToPrim_Side is: ", f12.6, " for ",2(i3),". Should be: 0.0625.")') UPrim(ENER,i,j,sideIdx,elemIdx),i,j
                ERROR STOP err_msg
            END IF
        END DO; END DO ! i,j=0,Nloc

    END SUBROUTINE TestConsToPrim_Side

    !==================================================================================================================================
    !> @brief Test ConsToPrim_Side method
    !==================================================================================================================================
    SUBROUTINE TestConsToPrim_Sides
    ! LOCAL VARIABLES
    INTEGER :: i,j, iSide
    INTEGER :: sideIdx = 1
    INTEGER :: elemIdx = 12
    INTEGER :: nSides = 3
    CHARACTER(LEN=200) :: err_msg
    !=========================================================================================

        ! Reset the solution array
        UPrim = 100.d0
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_UPrim, UPrim, SIZE(UPrim))
#endif

        ! Call sides method
        CALL ConsToPrim(5, UPrim(:,:,:,sideIdx:sideIdx+nSides-1,elemIdx), U(:,:,:,sideIdx:sideIdx+nSides-1,elemIdx), d_UPrim, d_U, elemIdx, sideIdx, nSides, STREAM_DEFAULT)

        ! If running on the device, copy data back over
#if (USE_ACCEL != ACCEL_OFF)
        CALL SynchronizeDevice()
        CALL CopyFromDevice(UPrim, d_UPrim, SIZE(UPrim))
#endif

        ! Check values are correct
        DO iSide = sideIdx,sideIdx+nSides-1
            DO j=0,5; DO i=0,5
                IF (UPrim(DENS,i,j,iSide,elemIdx) /= 2.0) THEN
                    write(err_msg, '("Calc of density in ConsToPrim_Sides is: ", f12.6, " for ",3(i3),". Should be: 2.0.")') UPrim(DENS,i,j,iSide,elemIdx),i,j,iSide
                    ERROR STOP err_msg
                END IF

                IF (UPrim(VEL1,i,j,iSide,elemIdx) /= 1.5) THEN
                    write(err_msg, '("Calc of x-velocity in ConsToPrim_Sides is: ", f12.6, " for ",3(i3),". Should be: 1.5.")') UPrim(VEL1,i,j,iSide,elemIdx),i,j,iSide
                    ERROR STOP err_msg
                END IF
        
                IF (UPrim(VEL2,i,j,iSide,elemIdx) /= 2.0) THEN
                    write(err_msg, '("Calc of y-velocity in ConsToPrim_Sides is: ", f12.6, " for ",3(i3),". Should be: 2.0.")') UPrim(VEL2,i,j,iSide,elemIdx),i,j,iSide
                    ERROR STOP err_msg
                END IF
        
                IF (UPrim(VEL3,i,j,iSide,elemIdx) /= 2.5) THEN
                    write(err_msg, '("Calc of z-velocity in ConsToPrim_Sides is: ", f12.6, " for ",3(i3),". Should be: 2.5.")') UPrim(VEL3,i,j,iSide,elemIdx),i,j,iSide
                    ERROR STOP err_msg
                END IF
        
                IF (UPrim(PRES,i,j,iSide,elemIdx) /= 1.25) THEN
                    write(err_msg, '("Calc of pressure in ConsToPrim_Sides is: ", f12.6, " for ",3(i3),". Should be: 1.25.")') UPrim(PRES,i,j,iSide,elemIdx),i,j,iSide
                    ERROR STOP err_msg
                END IF
        
                IF (UPrim(TEMP,i,j,iSide,elemIdx) /= 0.0625) THEN
                    write(err_msg, '("Calc of temperature in ConsToPrim_Sides is: ", f12.6, " for ",3(i3),". Should be: 0.0625.")') UPrim(ENER,i,j,iSide,elemIdx),i,j,iSide
                    ERROR STOP err_msg
                END IF
            END DO; END DO ! i,j=0,Nloc
        END DO

    END SUBROUTINE TestConsToPrim_Sides

    !==================================================================================================================================
    !> @brief Test ConsToPrim_Point method
    !==================================================================================================================================
    SUBROUTINE TestConsToPrim_Point
    ! LOCAL VARIABLES
    INTEGER :: i = 3, j = 2, k = 4
    INTEGER :: elemIdx = 12
    CHARACTER(LEN=200) :: err_msg
    !=========================================================================================

        ! Reset the solution array
        UPrim = 100.d0
#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_UPrim, UPrim, SIZE(UPrim))
#endif

        ! Call point method
        CALL ConsToPrim(5, UPrim(:,i,j,k,elemIdx), U(:,i,j,k,elemIdx), d_UPrim, d_U, elemIdx, i, j, k)

        ! If running on the device, copy data back over
#if (USE_ACCEL != ACCEL_OFF)
        CALL SynchronizeDevice()
        CALL CopyFromDevice(UPrim, d_UPrim, SIZE(UPrim))
#endif

        ! Check values are correct
        IF (UPrim(DENS,i,j,k,elemIdx) /= 2.0) THEN
            write(err_msg, '("Calc of density in ConsToPrim_Point is: ", f12.6,". Should be: 2.0.")') UPrim(DENS,i,j,k,elemIdx)
            ERROR STOP err_msg
        END IF

        IF (UPrim(VEL1,i,j,k,elemIdx) /= 1.5) THEN
            write(err_msg, '("Calc of x-velocity in ConsToPrim_Point is: ", f12.6,". Should be: 1.5.")') UPrim(VEL1,i,j,k,elemIdx)
            ERROR STOP err_msg
        END IF

        IF (UPrim(VEL2,i,j,k,elemIdx) /= 2.0) THEN
            write(err_msg, '("Calc of y-velocity in ConsToPrim_Point is: ", f12.6,". Should be: 2.0.")') UPrim(VEL2,i,j,k,elemIdx)
            ERROR STOP err_msg
        END IF

        IF (UPrim(VEL3,i,j,k,elemIdx) /= 2.5) THEN
            write(err_msg, '("Calc of z-velocity in ConsToPrim_Point is: ", f12.6,". Should be: 2.5.")') UPrim(VEL3,i,j,k,elemIdx)
            ERROR STOP err_msg
        END IF

        IF (UPrim(PRES,i,j,k,elemIdx) /= 1.25) THEN
            write(err_msg, '("Calc of pressure in ConsToPrim_Point is: ", f12.6,". Should be: 1.25.")') UPrim(PRES,i,j,k,elemIdx)
            ERROR STOP err_msg
        END IF

        IF (UPrim(TEMP,i,j,k,elemIdx) /= 0.0625) THEN
            write(err_msg, '("Calc of temperature in ConsToPrim_Point is: ", f12.6,". Should be: 0.0625.")') UPrim(ENER,i,j,k,elemIdx)
            ERROR STOP err_msg
        END IF

    END SUBROUTINE TestConsToPrim_Point


#if PP_EntropyVars == 1
    !==================================================================================================================================
    !> @brief Test ConsToEntropy_Volume method
    !==================================================================================================================================
    SUBROUTINE TestConsToEntropy_Volume
    ! LOCAL VARIABLES
    INTEGER :: i,j,k,n,iElem
    REAL :: ans(PP_nVar)
    CHARACTER(LEN=200) :: err_msg
    !=========================================================================================

        ! Init the entropy array for this test
        V(:,:,:,:,:) = 100.d0

#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_V, V, SIZE(V))
#endif

        ! Call method
        Call ConsToEntropy(5,V,U,d_V,d_U)

        ! If running on the device, copy data back over
#if (USE_ACCEL != ACCEL_OFF)
        CALL SynchronizeDevice()
        CALL CopyFromDevice(V, d_V, SIZE(V))
#endif

        ! Check the solution
        ans = (/ -7.0599927415085, 2.4d0, 3.2d0, 4.d0, -1.6d0 /)
        DO iElem=1,nElems
            DO k=0,5; DO j=0,5; DO i=0,5
                DO n = 1,PP_nVar
                    IF (.NOT. ALMOSTEQUALABSORREL(V(n,i,j,k,iElem),ans(n),200*PP_RealTolerance)) THEN
                        write(err_msg, '("Calc of variable ",i1," in ConsToEntropy_Volume is: ", f12.6, " for ",4(i3),". Should be:", f12.6)') n,V(n,i,j,k,iElem),i,j,k,iElem,ans(n)
                        ERROR STOP err_msg
                    END IF
                END DO
            END DO; END DO; END DO! i,j,k=0,Nloc
        END DO ! iElem 


    END SUBROUTINE TestConsToEntropy_Volume


    !==================================================================================================================================
    !> @brief Test EntropyToCons_Side method
    !==================================================================================================================================
    SUBROUTINE TestEntropyToCons_Side
    ! LOCAL VARIABLES
    INTEGER :: i,j,n
    REAL :: ans(PP_nVar)
    CHARACTER(LEN=200) :: err_msg
    !=========================================================================================

        ! Reset U and init V
        V(DENS,:,:,:,:) = 2.d0
        V(MOM1,:,:,:,:) = 3.d0
        V(MOM2,:,:,:,:) = 4.d0
        V(MOM3,:,:,:,:) = 5.d0
        V(ENER,:,:,:,:) = 15.d0
        U(:,:,:,:,:) = 100.d0

#if (USE_ACCEL != ACCEL_OFF)
        CALL CopyToDevice(d_V, V, SIZE(V))
        CALL CopyToDevice(d_U, U, SIZE(U))
#endif

        ! Call method
        ! Pass SideID of 4 instead of 3 to SideID because device backend assumes SideID is passed as 1 indexed
        ! but is 0 indexed here
        CALL EntropyToCons(5,4,1,V(:,:,:,3,1),U(:,:,:,3,1),d_V,d_U)

        ! If running on the device, copy data back over
#if (USE_ACCEL != ACCEL_OFF)
        CALL SynchronizeDevice()
        CALL CopyFromDevice(U, d_U, SIZE(U))
#endif

        ! Check solution
        ans = (/ -6.2958534279187d-3, 1.2591706855837d-3, 1.6788942474449d-3, 2.0986178093062d-3, 1.3990785395374d-4 /)
        DO j=0,5; DO i=0,5
            DO n = 1,PP_nVar
                IF (.NOT. ALMOSTEQUALABSORREL(U(n,i,j,3,1),ans(n),200*PP_RealTolerance)) THEN
                    write(err_msg, '("Calc of var ",i1," in EntropyToCons_Side is: ", f12.6, " for ",2(i3),". Should be: ", f12.6)') n,U(n,i,j,3,1),i,j,ans(n)
                    ERROR STOP err_msg
                END IF
            END DO
        END DO; END DO ! i,j=0,Nloc

    END SUBROUTINE TestEntropyToCons_Side
#endif /* PP_EntropyVars == 1 */

END MODULE TestEOS


!==================================================================================================================================
!> @brief Test driver program
!==================================================================================================================================
PROGRAM EOS_TestDriver
USE TestEOS
IMPLICIT NONE

    CALL InitTest
    CALL TestConsToPrim_Volume
    CALL TestConsToPrim_Elem
    CALL TestConsToPrim_Side
    CALL TestConsToPrim_Sides
    CALL TestConsToPrim_Point
#if PP_EntropyVars == 1
    CALL TestConsToEntropy_Volume
    CALL TestEntropyToCons_Side
#endif

END PROGRAM EOS_TestDriver
