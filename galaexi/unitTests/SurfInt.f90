#include "flexi.h"

!==================================================================================================================================
!> Unit test 'SurfIntUnitTest'
!> Test the routine: 'SurfInt', from module: 'SurfInt'.
!> Compare against precomputed and stored values.
!==================================================================================================================================
PROGRAM SurfIntUnitTest
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Unittest_Vars
USE MOD_Unittest,           ONLY: ReadInReferenceElementData
USE MOD_SurfIntCons,        ONLY: SurfIntCons
! Modules needed to read in reference element
USE MOD_DG_Vars,            ONLY: L_HatPlus,L_HatMinus
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_w,FV_w_inv, FV_Elems, FV_Elems_master,FV_Elems_slave
#endif
#if ((PP_NodeType==1) && defined(SPLIT_DG))
USE MOD_DG_Vars,            ONLY: U,UPrim,U_master,UPrim_master
USE MOD_Mesh_Vars,          ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems,Ja_Face,Ja_slave
USE MOD_SplitFlux,          ONLY: InitSplitDG
#endif /*((PP_NodeType==1) && defined(SPLIT_DG))*/

#if (USE_ACCEL != ACCEL_OFF)
USE MOD_Device
USE MOD_DeviceMem
USE MOD_DG_Vars      ,ONLY: d_Ut,d_Flux_master, d_Flux_slave
#if FV_ENABLED
USE MOD_FV_Vars      , ONLY: d_FV_Elems_master, d_FV_Elems_slave
#endif
#if ((PP_NodeType==1) && defined(SPLIT_DG))
USE MOD_Mesh_Vars    ,ONLY: d_Metrics_fTilde,d_Metrics_gTilde,d_Metrics_hTilde,d_Ja_face,d_Ja_slave
USE MOD_DG_Vars      ,ONLY: d_U,d_UPrim,d_U_master,d_U_slave,d_UPrim_master,d_UPrim_slave
#endif /* ((PP_NodeType==1) && defined(SPLIT_DG)) */
#endif /* USE_ACCEL */

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Flux(0:NRef,0:NRefZ,1:nSidesRef)
REAL                           :: Flux_nVar(PP_nVar,0:NRef,0:NRefZ,1:nSidesRef)
REAL                           :: Ut(PP_nVar,0:NRef,0:NRef,0:NRefZ,nElemsRef),Ut_ref(1,0:NRef,0:NRef,0:NRefZ,nElemsRef)
#if FV_ENABLED
REAL                           :: FV_Ut(PP_nVar,0:NRef,0:NRef,0:NRefZ,nElemsRef),FV_Ut_ref(1,0:NRef,0:NRef,0:NRefZ,nElemsRef)
#endif
! Therefore the data in CurvedSingleElementData.bin should be generate without formerly firstMasterSideID...
INTEGER                        :: i,j,k,l,nArgs
CHARACTER(LEN=*),PARAMETER     :: BinaryFluxString='SurfIntFlux.bin'
LOGICAL                        :: binaryExists,doGenerateReference=.FALSE.,equal,doGenerateFlux=.FALSE.
CHARACTER(LEN=255)             :: BinaryString,argument
!==================================================================================================================================
! Set file name for different node types
#if (PP_NodeType==1)
#ifdef SPLIT_DG
BinaryString='SurfInt_G3D_Split.bin'
#else
BinaryString='SurfInt_G3D.bin'
#endif /*SPLIT_DG*/
#elif (PP_NodeType==2)
#ifdef EXACT_MM
BinaryString='SurfInt_GL3D_EMM.bin'
#else
BinaryString='SurfInt_GL3D.bin'
#endif /*EXACT_MM*/
#endif /*PP_NodeType*/

WRITE(*,*) "Binary String: ",BinaryString

! Check for command line arguments to generate the reference solution
nArgs=COMMAND_ARGUMENT_COUNT()
IF (nArgs.GT.0) THEN
  CALL GET_COMMAND_ARGUMENT(1,argument)
  IF (argument.EQ.TRIM('--generate-reference')) THEN
    doGenerateReference = .TRUE.
  ELSE IF (argument.EQ.TRIM('--generate-flux')) THEN
    doGenerateFlux = .TRUE.
  ELSE
    WRITE(*,*) 'ERROR - Unknown command line argument.'
    STOP -1
  END IF
END IF

IF (doGenerateFlux) THEN
  ! Generate a random flux
  CALL RANDOM_NUMBER(Flux)
  ! Save the calculated flux to a binary file for later input
  OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(BinaryFluxString),FORM='unformatted')  ! replace an existing file or create a new one
  WRITE(10) Flux
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference flux to file ',BinaryFluxString
ELSE
  ! Read in the random flux
  OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryFluxString),FORM='unformatted')  ! open an existing file
  READ(10) Flux
  CLOSE(10) ! close the file
  ! Build flux on PP_nVar as expected by SurfInt
  DO i=1,PP_nVar
    Flux_nVar(i,:,:,:) = Flux(:,:,:)
  END DO
END IF

! Regardless of the formulation, we need Ut and master/slave fluxes 
#if (USE_ACCEL != ACCEL_OFF)
  CALL AllocateDeviceMemory(d_Flux_master, SIZE_C_DOUBLE, SIZE(Flux_nVar))
  CALL CopyToDevice(d_Flux_master, Flux_nVar, SIZE(Flux_nVar))

  CALL AllocateDeviceMemory(d_Ut,SIZE_C_DOUBLE, SIZE(Ut))
  CALL CopyToDevice(d_Ut, Ut, SIZE(Ut))
#endif

! Read in data from single curved element
CALL ReadInReferenceElementData()

! Initialize Ut
Ut = 0.

#if FV_ENABLED
FV_Ut = 0.
ALLOCATE(FV_w(0:NRef))        ! 1D width of FV-Subcells
ALLOCATE(FV_w_inv(0:NRef))
FV_w(:)     = 2./(9+1) ! equidistant widths of FV-Subcells
FV_w_inv(:) = 1./FV_w(:)
ALLOCATE(FV_Elems(1:1))
ALLOCATE(FV_Elems_master(1:nSidesRef))
ALLOCATE(FV_Elems_slave (1:nSidesRef))
FV_Elems = 0
FV_Elems_master = 0
FV_Elems_slave = 0
#endif
#if ((PP_NodeType==1) && defined(SPLIT_DG))
ALLOCATE(U(PP_nVar,0:NRef,0:NRef,0:NRefZ,nElems))
ALLOCATE(UPrim(PP_nVarPrim,0:NRef,0:NRef,0:NRefZ,nElems))
ALLOCATE(U_master(PP_nVar,0:NRef,0:NRefZ,nElems))
ALLOCATE(UPrim_master(PP_nVarPrim,0:NRef,0:NRefZ,nSidesRef))
U = 1.; UPrim=1.; U_master=1.; UPrim_master=1.
ALLOCATE(Metrics_fTilde(3,0:NRef,0:NRef,0:NRefZ,nElems,0:0))
ALLOCATE(Metrics_gTilde(3,0:NRef,0:NRef,0:NRefZ,nElems,0:0))
Metrics_fTilde=1.; Metrics_gTilde=1.
ALLOCATE(Metrics_hTilde(3,0:NRef,0:NRef,0:NRefZ,nElems,0:0))
Metrics_hTilde=1.
ALLOCATE( Ja_Face(3,3,0:NRef,0:NRefZ,1:nSidesRef)) ! temp
ALLOCATE(Ja_slave(3,3,0:NRef,0:NRefZ,1:nSidesRef)) ! temp
Ja_slave = 1.
Ja_Face  = 1.
CALL InitSplitDG(0)

! Allocate and copy device variables if needed
#if (USE_ACCEL != ACCEL_OFF)
CALL AllocateDeviceMemory(d_U, SIZE_C_DOUBLE, SIZE(U))
CALL CopyToDevice(d_U, U, SIZE(U))
CALL AllocateDeviceMemory(d_U_master, SIZE_C_DOUBLE, SIZE(U_master))
CALL CopyToDevice(d_U_master, U_master, SIZE(U_master))

CALL AllocateDeviceMemory(d_UPrim, SIZE_C_DOUBLE, SIZE(UPrim))
CALL CopyToDevice(d_UPrim, UPrim, SIZE(UPrim))
CALL AllocateDeviceMemory(d_UPrim_master, SIZE_C_DOUBLE, SIZE(UPrim_master))
CALL CopyToDevice(d_UPrim_master, UPrim_master, SIZE(UPrim_master))

CALL AllocateDeviceMemory(d_Metrics_fTilde, SIZE_C_DOUBLE, SIZE(Metrics_fTilde))
CALL CopyToDevice(d_Metrics_fTilde, Metrics_fTilde, SIZE(Metrics_fTilde))
CALL AllocateDeviceMemory(d_Metrics_gTilde, SIZE_C_DOUBLE, SIZE(Metrics_gTilde))
CALL CopyToDevice(d_Metrics_gTilde, Metrics_gTilde, SIZE(Metrics_fTilde))
CALL AllocateDeviceMemory(d_Metrics_hTilde, SIZE_C_DOUBLE, SIZE(Metrics_hTilde))
CALL CopyToDevice(d_Metrics_hTilde, Metrics_hTilde, SIZE(Metrics_hTilde))

CALL AllocateDeviceMemory(d_Ja_Face, SIZE_C_DOUBLE, SIZE(Ja_Face))
CALL CopyToDevice(d_Ja_Face, Ja_Face, SIZE(Ja_Face))
CALL AllocateDeviceMemory(d_Ja_slave, SIZE_C_DOUBLE, SIZE(Ja_slave))
CALL CopyToDevice(d_Ja_slave, Ja_slave, SIZE(Ja_slave))
#endif /* USE_ACCEL != ACCEL_OFF */
#endif /*((PP_NodeType==1) && defined(SPLIT_DG))*/

! Call SurfInt
CALL SurfIntCons(NRef,Flux_nVar,Flux_nVar,Ut,.FALSE.,L_HatMinus,L_HatPlus)
#if FV_ENABLED
FV_Elems = 1
FV_Elems_master = 1
FV_Elems_slave = 1
CALL SurfIntCons(NRef,Flux_nVar,Flux_nVar,FV_Ut,.FALSE.,L_HatMinus,L_HatPlus)
#endif

! If testing device code, get solution back from device
#if (USE_ACCEL != ACCEL_OFF)
CALL CopyFromDevice(Ut, d_Ut, SIZE(Ut))
#endif
#if ((PP_NodeType==1) && defined(SPLIT_DG))
DO i=2,PP_nVar
  Ut(i,:,:,:,:) = Ut(1,:,:,:,:)
END DO
#endif /*((PP_NodeType==1) && defined(SPLIT_DG))*/


IF (doGenerateReference) THEN
  ! Save the calculated solution to a binary file for later comparison
  OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(BinaryString),FORM='unformatted')  ! replace an existing file or create a new one
  WRITE(10) Ut(1,:,:,:,:)
#if FV_ENABLED
  WRITE(10) FV_Ut(1,:,:,:,:)
#endif
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference to file ',BinaryString
ELSE
  ! Check if binary results file exists
  INQUIRE(FILE=TRIM(BinaryString),EXIST=binaryExists)

  IF (binaryExists) THEN
    ! Read the reference solution
    OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryString),FORM='unformatted')  ! open an existing file
    READ(10) Ut_ref
#if FV_ENABLED
    READ(10) FV_Ut_ref
#endif
    CLOSE(10) ! close the file
    ! Check if the computed and the reference solutions are within a given tolerance
    equal =  .TRUE.
    DO i=1,PP_nVar; DO j=0,NRef; DO k=0,NRef; DO l=0,NRefZ
#if ((PP_NodeType==1) && defined(SPLIT_DG))
      equal = ALMOSTEQUALABSORREL(Ut(i,j,k,l,1),Ut_ref(1,j,k,l,1),200.*PP_RealTolerance)
#else
      equal = ALMOSTEQUALABSORREL(Ut(i,j,k,l,1),Ut_ref(1,j,k,l,1),100.*PP_RealTolerance)
#endif
      IF (.NOT. equal) THEN
        WRITE(*,*) 'ERROR - Calculated surface integral deviates from reference for variable: ', i, " at DOF: ", j,k,l
        WRITE(*,*) "Calculated solution is: ", Ut(i,j,k,l,1), " Should be: ", Ut_ref(1,j,k,l,1)
        STOP -1
      END IF
#if FV_ENABLED
      equal = ALMOSTEQUALABSORREL(FV_Ut(i,j,k,l,1),FV_Ut_ref(1,j,k,l,1),100.*PP_RealTolerance)
#endif
    END DO; END DO; END DO; END DO
    
    WRITE(*,*) 'Compared surface integral against stored data in ',TRIM(BinaryString), ' -- SUCCESSFUL.'
  ELSE
    WRITE(*,*) 'ERROR - No reference solution has been found.'
    STOP -1
  END IF
END IF

END PROGRAM SurfIntUnitTest
