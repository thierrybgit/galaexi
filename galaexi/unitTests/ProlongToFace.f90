#include "flexi.h"

!==================================================================================================================================
!> Unit test 'ProlongToFaceUnitTest'
!> Test the routine: 'ProlongToFace', from module: 'ProlongToFace'.
!> Compare against precomputed and stored values.
!==================================================================================================================================
PROGRAM ProlongToFaceUnitTest
! MODULES
USE ISO_C_BINDING,          ONLY: C_NULL_CHAR
USE MOD_Globals
USE MOD_PreProc
USE MOD_Unittest_Vars
USE MOD_Unittest,           ONLY: ReadInReferenceElementData
USE MOD_ProlongToFace
! Modules needed to read in reference element
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems,FV_Elems_master,FV_Elems_slave
#endif
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_DeviceMem
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Uvol(0:NRef,0:NRef,0:NRefZ,nElemsRef)
REAL                           :: Uvol_nVar(PP_nVar,0:NRef,0:NRef,0:NRefZ,nElemsRef)
REAL                           :: Uface_master(PP_nVar,0:NRef,0:NRefZ,1:nSidesRef),Uface_master_ref(0:NRef,0:NRefZ,1:nSidesRef)
REAL                           :: Uface_slave(PP_nVar,0:NRef,0:NRefZ,1:nSidesRef)!,Uface_slave_ref(0:9,0:9,7:6)
#if FV_ENABLED
REAL                           :: FV_Uface_master(PP_nVar,0:NRef,0:NRefZ,1:nSidesRef),FV_Uface_master_ref(0:NRef,0:NRefZ,1:nSidesRef)
REAL                           :: FV_Uface_slave(PP_nVar,0:NRef,0:NRefZ,1:nSidesRef)
#endif
INTEGER                        :: i,j,k,l,nArgs
CHARACTER(LEN=*),PARAMETER     :: BinaryUvolString='ProlongToFaceUvol.bin'
LOGICAL                        :: binaryExists,doGenerateReference=.FALSE.,equal,doGenerateUvol=.FALSE.
CHARACTER(LEN=255)             :: BinaryString,argument
INTEGER(C_INT)                    :: d_Uvol_nVar 
INTEGER(C_INT)                    :: d_Uface_master 
INTEGER(C_INT)                    :: d_Uface_slave 
#if FV_ENABLED
INTEGER(C_INT)                    :: d_FV_Uface_master 
INTEGER(C_INT)                    :: d_FV_Uface_slave 
#endif
!==================================================================================================================================
! Set binary file name for different node types
#if (PP_NodeType==1)
BinaryString='ProlongToFace_G3D.bin'
#elif (PP_NodeType==2)
BinaryString='ProlongToFace_GL3D.bin'
#endif

! Check for command line arguments to generate the reference solution
nArgs=COMMAND_ARGUMENT_COUNT()
IF (nArgs.GT.0) THEN
  CALL GET_COMMAND_ARGUMENT(1,argument)
  IF (argument.EQ.TRIM('--generate-reference')) THEN
    doGenerateReference = .TRUE.
  ELSE IF (argument.EQ.TRIM('--generate-uvol')) THEN
    doGenerateUvol = .TRUE.
  ELSE
    WRITE(*,*) 'ERROR - Unknown command line argument.'
    STOP -1
  END IF
END IF

IF (doGenerateUvol) THEN
  ! Generate a random volume solution
  CALL RANDOM_NUMBER(Uvol)
  ! Save the calculated volume solution to a binary file for later input
  OPEN(UNIT=10, STATUS='replace',FILE=TRIM(BinaryUvolString),FORM='unformatted')  ! replace an existing file or create a new one
  WRITE(10) Uvol
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference Uvol to file ',BinaryUvolString
ELSE
  ! Read in the random Uvol
  OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryUvolString),FORM='unformatted')  ! open an existing file
  READ(10) Uvol
  CLOSE(10) ! close the file
  ! Build Uvol on PP_nVar as expected by ProlongToFace
  DO i=1,PP_nVar
    Uvol_nVar(i,:,:,:,:) = Uvol(:,:,:,:)
  END DO
END IF

! Read in data from single curved element
CALL ReadInReferenceElementData()

#if FV_ENABLED
ALLOCATE(FV_Elems(1:1))
ALLOCATE(FV_Elems_master(1:nSidesRef))
ALLOCATE(FV_Elems_slave (1:nSidesRef))
FV_Elems = 0
FV_Elems_master = 0
FV_Elems_slave = 0
#endif /* FV_ENABLED */

! If this is a device accelerated test, move all information over to the GPU
#if (USE_ACCEL != ACCEL_OFF)
CALL AllocateDeviceMemory(d_Uvol_nVar,SIZE_C_DOUBLE,SIZE(Uvol_nVar))
CALL CopyToDevice(d_Uvol_nVar,Uvol_nVar,SIZE(Uvol_nVar))
! The sides arrays, we only need to allocate, there is no need to copy, 
! as they are set by the call to ProlongToFace
CALL AllocateDeviceMemory(d_Uface_master,SIZE_C_DOUBLE,SIZE(Uface_master))
CALL AllocateDeviceMemory(d_Uface_slave,SIZE_C_DOUBLE,SIZE(Uface_slave))
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
! The above comment is invalid for hybrid chips we have to call the copy routines
! so we associated the device pointer with the host memory
CALL CopyToDevice(d_Uface_master, Uface_master, SIZE(Uface_master))
CAll CopyToDevice(d_Uface_slave, Uface_slave, SIZE(Uface_slave))
#endif /* USE_HYBRID */

#if FV_ENABLED
! Don't need to allocate or copy FV_Elems arrays, as everything they are used for here is on the host

! This arrays are only allocated and not copied (same reason as d_Uface_... above)
CALL AllocateDeviceMemory(d_FV_Uface_master,SIZE_C_INT,SIZE(FV_Uface_master))
CALL AllocateDeviceMemory(d_FV_Uface_slave,SIZE_C_INT,SIZE(FV_Uface_slave))
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
! We have to call copies for hybrid chips (same reason as d_Uface_... above)
CALL CopyToDevice(d_FV_Uface_master, FV_Uface_master, SIZE(FV_Uface_master))
CAll CopyToDevice(d_FV_Uface_slave, FV_Uface_slave, SIZE(FV_Uface_slave))
#endif /* USE_HYBRID */
#endif /* FV_ENABLED */
#endif /* USE_ACCEL */

! Call ProlongToFace
CALL ProlongToFace(PP_nVar,NRef,Uvol_nVar,Uface_master,Uface_slave, &
                              d_Uvol_nVar,d_Uface_master,d_Uface_slave,L_Minus,L_Plus,.FALSE.)

#if FV_ENABLED
FV_Elems = 1
FV_Elems_master = 1
FV_Elems_slave = 1
CALL ProlongToFace(PP_nVar,NRef,Uvol_nVar,FV_Uface_master,FV_Uface_slave, &
                                        d_Uvol_nVar,d_FV_Uface_master,d_FV_Uface_slave,L_Minus,L_Plus,.FALSE.)
#endif /* FV_ENABLED */

! Get solution from the device, if needed
#if (USE_ACCEL != ACCEL_OFF)
CALL CopyFromDevice(Uface_master, d_Uface_master, SIZE(Uface_master))
! slave array isn't checked below, so don't bother copying back
#if FV_ENABLED
CALL CopyFromDevice(FV_Uface_master, d_FV_Uface_master, SIZE(FV_Uface_master))
#endif /* FV_ENABLED */
#endif /* USE_ACCEL */

IF (doGenerateReference) THEN
  ! Save the calculated solution to a binary file for later comparison
  OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(BinaryString),FORM='unformatted')  ! replace an existing file or create a new one
  WRITE(10) Uface_master(1,:,:,:)!,Uface_slave(1,:,:,:)
#if FV_ENABLED
  WRITE(10) FV_Uface_master(1,:,:,:)
#endif
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference to file ',BinaryString
ELSE
  ! Check if binary results file exists
  INQUIRE(FILE=TRIM(BinaryString),EXIST=binaryExists)

  IF (binaryExists) THEN
    ! Read the reference solution
    OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryString),FORM='unformatted')  ! open an existing file
    READ(10) Uface_master_ref!,Uface_slave_ref
#if FV_ENABLED
    READ(10) FV_Uface_master_ref
#endif
    CLOSE(10) ! close the file
    ! Check if the computed and the reference solutions are within a given tolerance
    equal =  .TRUE.
    DO i=1,PP_nVar; DO j=0,NRef; DO k=0,NRefZ; DO l=1,nSidesRef
      equal = ALMOSTEQUALABSORREL(Uface_master(i,j,k,l),Uface_master_ref(j,k,l),100.*PP_RealTolerance)
      IF (.NOT. equal) THEN
        WRITE(*,'(A,4(i2),A)') 'Calculated prolonged values deviate from reference for: (',i,j,k,l,' ).'
        WRITE(*,'(A,f12.6,A,f12.6)') "Value should be: ", Uface_master_ref(j,k,l), " is: ", Uface_master(i,j,k,l)
        ERROR STOP
      END IF
#if FV_ENABLED
      equal = ALMOSTEQUALABSORREL(FV_Uface_master(i,j,k,l),FV_Uface_master_ref(j,k,l),100.*PP_RealTolerance)
      IF (.NOT. equal) THEN
        WRITE(*,'(A,4(i2),A)') 'Calculated prolonged FV values deviate from reference for: (',i,j,k,l,' ).'
        WRITE(*,'(A,f12.6,A,f12.6)') "Value should be: ", FV_Uface_master_ref(j,k,l), " is: ", FV_Uface_master(i,j,k,l)
        ERROR STOP
      END IF
#endif
    END DO; END DO; END DO; END DO
    ! Plus sides not needed in single element case
    !DO i=1,PP_nVar; DO j=0,9; DO k=0,9; DO l=7,6
      !equal = ALMOSTEQUALABSORREL(Uface_slave(i,j,k,l),Uface_slave_ref(j,k,l),100.*PP_RealTolerance) .AND. equal
    !END DO; END DO; END DO; END DO
  ELSE
    WRITE(*,*) 'ERROR - No reference solution has been found.'
    ERROR STOP
  END IF
END IF

END PROGRAM ProlongToFaceUnitTest
