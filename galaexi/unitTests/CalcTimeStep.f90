#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Unit test 'CalcTimeStep'
!> Test the routines from module: 'MOD_CalcTimeStep'.
!> @author Spencer Starr - 08.2024
!==================================================================================================================================

MODULE TestCalcTimeStep
USE ISO_C_BINDING
USE MOD_PreProc
USE MOD_Unittest_Vars
USE MOD_Unittest,           ONLY: ReadInReferenceElementData
USE MOD_CalcTimeStep
USE MOD_Device
USE MOD_DG_Vars,            ONLY: U, d_U
USE MOD_Mesh_Vars,          ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,Elem_xGP,nElems
USE MOD_TimeDisc_Vars,      ONLY: CFLScale,dtElem,d_CFLScale
USE MOD_EOS_Vars,           ONLY: R
#if PARABOLIC
USE MOD_TimeDisc_Vars,      ONLY: DFLScale, d_DFLScale
USE MOD_EOS_Vars,           ONLY: KappasPr
#endif /*PARABOLIC*/
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems
#if FV_ENABLED == 2
USE MOD_FV_Vars,            ONLY: FV_alpha,FV_alpha_min
#endif /* FV_ENABLED == 2 */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,      ONLY: muSGS
#endif
! Device stuff
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_DeviceMem
USE MOD_DG_Vars,            ONLY: d_U
USE MOD_Mesh_Vars,          ONLY: d_Metrics_fTilde,d_Metrics_gTilde,d_Metrics_hTilde
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: d_FV_Elems
#if FV_ENABLED == 2
USE MOD_FV_Vars,            ONLY: d_FV_alpha
#endif /* FV_ENABLED == 2 */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,      ONLY: d_muSGS
#endif
#endif /* USE_ACCEL */
IMPLICIT NONE


CONTAINS

    !==================================================================================================================================
    !> @brief Initialize the domain for the test
    !==================================================================================================================================
    SUBROUTINE InitTest
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    REAL                           :: Uvol(0:NRef,0:NRef,0:NRef,nElemsRef)
    INTEGER                        :: i
    CHARACTER(LEN=*),PARAMETER     :: BinaryUString='ProlongToFaceUvol.bin'
    !==================================================================================================================================

        ! Allocate and init host side arrays - use 16 5th order elements
        nElems = nElemsRef
        ALLOCATE( U(PP_nVar,0:NRef,0:NRef,0:NRef,nElemsRef) )
        ALLOCATE( Metrics_fTilde(3,0:NRef,0:NRef,0:NRef,nElemsRef,0:FV_SIZE) )
        ALLOCATE( Metrics_gTilde(3,0:NRef,0:NRef,0:NRef,nElemsRef,0:FV_SIZE) )
        ALLOCATE( Metrics_hTilde(3,0:NRef,0:NRef,0:NRef,nElemsRef,0:FV_SIZE) )
        ALLOCATE( Elem_xGP(3,0:NRef,0:NRef,0:NRef,nElemsRef) )
        ALLOCATE( dtElem(nElemsRef) )
#if EDDYVISCOSITY
        ALLOCATE( muSGS(1,0:NRef,0:NRef,0:NRef,nElemsRef) )
        muSGS = 1.d0
#endif

#if FV_ENABLED
        ALLOCATE(FV_Elems(nElemsRef))
        FV_Elems = 0
#endif

        ! Intialize the host side arrays
        ! Read in the random Uvol. Just use the ProlongToFace solution that already exists for simplicity

        OPEN(UNIT=11, STATUS='old',FILE=TRIM(BinaryUString),FORM='unformatted')  ! open an existing file
        READ(11) Uvol
        CLOSE(11) ! close the file
        ! Build Uvol on PP_nVar as expected by ProlongToFace
        DO i=1,PP_nVar
            U(i,:,:,:,:) = Uvol(:,:,:,:)
        END DO

        ! Read in data from single curved element 
        CALL ReadInReferenceElementData()

        Metrics_fTilde = 1.d0
        Metrics_gTilde = 1.d0
        Metrics_hTilde = 1.d0

        DO i = 0,FV_SIZE
            CFLScale(i) = 0.5d0
#if PARABOLIC
            DFLScale(i) = 0.5d0
#endif /*PARABOLIC*/
        END DO

        R = 1.d0
#if PARABOLIC
        KappasPr = 1.d0
#endif /*PARABOLIC*/

        ! If device test, allocate and copy on device
#if (USE_ACCEL != ACCEL_OFF)
        CALL AllocateDeviceMemory(d_U, SIZE_C_DOUBLE, SIZE(U))
        CALL CopyToDevice(d_U, U, SIZE(U))

        CALL AllocateDeviceMemory(d_Metrics_fTilde, SIZE_C_DOUBLE, SIZE(Metrics_fTilde))
        CALL CopyToDevice(d_Metrics_fTilde, Metrics_fTilde, SIZE(Metrics_fTilde))
        CALL AllocateDeviceMemory(d_Metrics_gTilde, SIZE_C_DOUBLE, SIZE(Metrics_gTilde))
        CALL CopyToDevice(d_Metrics_gTilde, Metrics_gTilde, SIZE(Metrics_fTilde))
        CALL AllocateDeviceMemory(d_Metrics_hTilde, SIZE_C_DOUBLE, SIZE(Metrics_hTilde))
        CALL CopyToDevice(d_Metrics_hTilde, Metrics_hTilde, SIZE(Metrics_hTilde))

        CALL AllocateDeviceMemory(d_CFLScale, SIZE_C_DOUBLE, SIZE(CFLScale))
        CALL CopyToDevice(d_CFLScale, CFLScale, SIZE(CFLScale))

#if PARABOLIC
        CALL AllocateDeviceMemory(d_DFLScale, SIZE_C_DOUBLE, SIZE(DFLScale))
        CALL CopyToDevice(d_DFLScale, DFLScale, SIZE(DFLScale))
#endif

#if FV_ENABLED
        CALL AllocateDeviceMemory(d_FV_Elems, SIZE_C_INT, SIZE(FV_Elems))
        CALL CopyToDevice(d_FV_Elems, FV_Elems, SIZE(FV_Elems))

        CALL AllocateDeviceMemory(d_FV_alpha, SIZE_C_DOUBLE, SIZE(FV_alpha))
        CALL CopyToDevice(d_FV_alpha, FV_alpha, SIZE(FV_alpha))
#endif
#if EDDYVISCOSITY
        CALL AllocateDeviceMemory(d_muSGS, SIZE_C_DOUBLE, SIZE(muSGS))
        CALL CopyToDevice(d_muSGS, muSGS, SIZE(muSGS))
#endif 
#endif /* USE_ACCEL */
 
    END SUBROUTINE InitTest


    !==================================================================================================================================
    !> @brief Test calculation of the time step
    !==================================================================================================================================
    SUBROUTINE TestMethod
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    REAL :: ans
    REAL :: timeStep
    INTEGER :: err_code
    !==================================================================================================================================
 
        ! Set the expected solution based on physics
        ans = 1.3744260421477193d-2
        ! Call init method and then calc timestep
        CALL InitCalcTimestep
        CALL CalcTimeStep(timeStep,err_code,PP_N)

        ! Check solution
        IF (.NOT. ALMOSTEQUALABSORREL(timeStep,ans,200*PP_RealTolerance)) THEN
            PRINT *, "Value for time step from CalcTimeStep for variable is: ", timeStep, " should be: ", ans
            ERROR STOP
        END IF

    END SUBROUTINE TestMethod

END MODULE TestCalcTimeStep



!==================================================================================================================================
!> @brief Test driver program
!==================================================================================================================================
PROGRAM CalcTimeStep_TestDriver
USE TestCalcTimeStep
USE MOD_Globals
#if USE_MPI
USE MPI
#endif
IMPLICIT NONE

INTEGER :: ierr

#if USE_MPI
	CALL MPI_Init(ierr)
	MPI_COMM_FLEXI = MPI_COMM_WORLD
#endif

    CALL InitTest
    CALL TestMethod

#if USE_MPI
	CALL MPI_Finalize(ierr)
#endif

END PROGRAM CalcTimeStep_TestDriver
