#include "flexi.h"

!==================================================================================================================================
!> Module containing general routines used in Unittests
!==================================================================================================================================
MODULE MOD_Unittest
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE GenerateUnittestReferenceData
  MODULE PROCEDURE GenerateUnittestReferenceData
END INTERFACE

INTERFACE ReadInReferenceElementData
  MODULE PROCEDURE ReadInReferenceElementData
END INTERFACE

PUBLIC::GenerateUnittestReferenceData,ReadInReferenceElementData
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Writes several arrays to a binary file, that are used for the unittest. Therewith the unittest do not have to compute
!> geometrical quantities, index arrays, ...
!==================================================================================================================================
SUBROUTINE GenerateUnittestReferenceData()
! MODULES
USE MOD_Mesh_Vars
USE MOD_Interpolation_Vars
USE MOD_DG_Vars
!----------------------------------------------------------------------------------------------------------------------------------
! Local Variables
CHARACTER(LEN=255)             :: Filename
!==================================================================================================================================
#if FV_ENABLED == 0
WRITE(*,*) 'To generate reference data for the unit tests, please compile FLEXI with FV activated.'
STOP
#endif

#if PP_dim == 3
Filename = "UnittestElementData3D.bin"
#else
Filename = "UnittestElementData2D.bin"
#endif
! Save the calculated solution to a binary file for later comparison
OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(Filename),FORM='unformatted')  ! replace an existing file or create a new one
WRITE(10) nElems,SideToElem,firstMPISide_YOUR,lastMPISide_MINE,nSides,S2V2,S2V,L_Minus,L_Plus,L_HatPlus,L_HatMinus,sJ
CLOSE(10) ! close the file
WRITE(*,*) "Generated Unittest reference data into: ", TRIM(Filename)
END SUBROUTINE GenerateUnittestReferenceData

!==================================================================================================================================
!> Read in the data for the curved reference element, allocate the necessary arrays beforehand.
!==================================================================================================================================
SUBROUTINE ReadInReferenceElementData()
! MODULES
Use Iso_C_Binding
USE MOD_Unittest_Vars
USE MOD_Mesh_Vars,             ONLY: SideToElem,ElemToSide,S2V2,nElems,nSides,firstMPISide_YOUR,lastMPISide_MINE,sJ,S2V
USE MOD_Interpolation_Vars,    ONLY: L_Minus,L_Plus
USE MOD_DG_Vars,               ONLY: L_HatPlus,L_HatMinus
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_Device
USE MOD_DeviceMem
USE MOD_DG_Vars      ,ONLY: d_L_HatMinus, d_L_HatPlus
USE MOD_Mesh_Vars    ,ONLY: d_S2V,  d_SideToElem
USE MOD_Mesh_Vars    ,ONLY: d_ElemToSide,d_S2V2,d_sJ
USE MOD_Interpolation_Vars, ONLY: d_L_Minus,d_L_Plus
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! Local Variables
CHARACTER(LEN=255)             :: Filename
INTEGER                        :: Flip_lower,Flip_upper,locSide_lower,locSide_upper
INTEGER                        :: side
!==================================================================================================================================
! Dimensions for mappings
  Flip_lower = 0
  Flip_upper = 4
  locSide_lower = 1
  locSide_upper = 6

  Filename = "UnittestElementData3D.bin"

  ! Read in data from single curved element
  ALLOCATE(SideToElem(1:5,1:nSidesRef))
  ALLOCATE(S2V2(1:2,0:NRef,0:NRefZ,Flip_lower:Flip_upper,locSide_lower:locSide_upper))
  ALLOCATE(S2V(1:3,0:NRef,0:NRef,0:NRefZ,Flip_lower:Flip_upper,locSide_lower:locSide_upper))
  ALLOCATE(L_Minus(0:NRef))
  ALLOCATE(L_Plus(0:NRef))
  ALLOCATE(L_HatMinus(0:NRef))
  ALLOCATE(L_HatPlus(0:NRef))
  ALLOCATE(sJ(0:NRef,0:NRef,0:NRef,nElemsRef,0:FV_SIZE))
  OPEN(UNIT = 10, STATUS='old',FILE=TRIM(Filename),FORM='unformatted')  ! open an existing file
  READ(10) nElems,SideToElem,firstMPISide_YOUR,lastMPISide_MINE,nSides,S2V2,S2V,L_Minus,L_Plus,L_HatPlus,L_HatMinus,sJ
  CLOSE(10) ! close the file

#if (USE_ACCEL != ACCEL_OFF)
  ! Need to calculate ElemToSide as well for SurfInt on device
  ALLOCATE( ElemToSide(2,6,nElemsRef) )
  DO side=1,6
    ElemToSide(E2S_SIDE_ID,side,1) = side
    ElemToSide(E2S_FLIP,side,1) = 0
  END DO ! LocSideID

  ! Allocate/copy example mesh to device if testing accelerated code
  CALL AllocateDeviceMemory(d_ElemToSide, SIZE_C_INT, SIZE(ElemToSide))
  CALL CopyToDevice(d_ElemToSide, C_Loc(ElemToSide), SIZE_C_INT, SIZE(ElemToSide))
  CALL AllocateDeviceMemory(d_SideToElem, SIZE_C_INT, SIZE(SideToElem))
  CALL CopyToDevice(d_SideToElem, C_Loc(SideToElem), SIZE_C_INT, SIZE(SideToElem))

  CALL AllocateDeviceMemory(d_S2V2, SIZE_C_INT, SIZE(S2V2))
  CALL CopyToDevice(d_S2V2, C_Loc(S2V2), SIZE_C_INT, SIZE(S2V2))
  CALL AllocateDeviceMemory(d_S2V, SIZE_C_INT, SIZE(S2V))
  CALL CopyToDevice(d_S2V, C_Loc(S2V), SIZE_C_INT, SIZE(S2V))

  CALL AllocateDeviceMemory(d_L_Minus, SIZE_C_DOUBLE, size(L_Minus))
  CALL CopyToDevice(d_L_Minus, C_Loc(L_Minus), SIZE_C_DOUBLE, size(L_Minus))
  CALL AllocateDeviceMemory(d_L_Plus, SIZE_C_DOUBLE, size(L_Plus))
  CALL CopyToDevice(d_L_Plus, C_Loc(L_Plus), SIZE_C_DOUBLE, size(L_Plus))

  CALL AllocateDeviceMemory(d_L_HatMinus,SIZE_C_DOUBLE, SIZE(L_HatMinus))
  CALL CopyToDevice(d_L_HatMinus, C_Loc(L_HatMinus), SIZE_C_DOUBLE, SIZE(L_HatMinus))
  CALL AllocateDeviceMemory(d_L_HatPlus,SIZE_C_DOUBLE, SIZE(L_HatPlus))
  CALL CopyToDevice(d_L_HatPlus, C_Loc(L_HatPlus), SIZE_C_DOUBLE, SIZE(L_HatPlus))

  CALL AllocateDeviceMemory(d_sJ,SIZE_C_DOUBLE, SIZE(sJ))
  CALL CopyToDevice(d_sJ, C_Loc(sJ), SIZE_C_DOUBLE, SIZE(sJ))
#endif

END SUBROUTINE ReadInReferenceElementData

END MODULE MOD_Unittest
