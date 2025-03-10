!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Contains routines to interpolate the interior solution to the boundary
!==================================================================================================================================
MODULE MOD_ProlongToFace
USE ISO_C_BINDING
IMPLICIT NONE
PRIVATE

INTERFACE ProlongToFace
  MODULE PROCEDURE ProlongToFaceEntropy
  MODULE PROCEDURE ProlongToFace
END INTERFACE

INTERFACE EvalElemFace
  MODULE PROCEDURE EvalElemFaceG_Host
  MODULE PROCEDURE EvalElemFaceGL_Host
END INTERFACE

! Device backend interface
#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
  SUBROUTINE ProlongToFace_Device(TP_nVar, Nloc, nElems, nSides, firstMPISide_YOUR, lastMPISide_MINE, &
                                  d_Uvol, d_Uface_master, d_Uface_slave, d_L_Minus, d_L_Plus,  &
                                  d_SideToElem, d_S2V2 &
#if FV_ENABLED
                                  ,d_isFV &
#endif
                                  ,streamID &
            ) BIND(C, NAME="ProlongToFace_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: TP_nVar
    INTEGER(C_INT),VALUE :: Nloc
    INTEGER(C_INT),VALUE :: nElems
    INTEGER(C_INT),VALUE :: nSides
    INTEGER(C_INT),VALUE :: firstMPISide_YOUR
    INTEGER(C_INT),VALUE :: lastMPISide_MINE
    INTEGER(C_INT),VALUE :: d_Uvol
    INTEGER(C_INT),VALUE :: d_Uface_master
    INTEGER(C_INT),VALUE :: d_Uface_slave
    INTEGER(C_INT),VALUE :: d_L_Minus
    INTEGER(C_INT),VALUE :: d_L_Plus
    INTEGER(C_INT),VALUE :: d_SideToElem
    INTEGER(C_INT),VALUE :: d_S2V2
#if FV_ENABLED
    INTEGER(C_INT),VALUE :: d_isFV
#endif
    INTEGER(C_INT),VALUE :: streamID
  END SUBROUTINE ProlongToFace_Device
END INTERFACE
#endif

PUBLIC :: ProlongToFace
PUBLIC :: EvalElemFace

CONTAINS

!----------------------------------------------------------------------------------------------------------------------------------
! ENTRY POINT METHODS
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
!> Entry point for ProlongToFace to hide backend split
!==================================================================================================================================
SUBROUTINE ProlongToFace(TP_nVar,Nloc,Uvol,Uface_master,Uface_slave,d_Uvol,d_Uface_Master,d_Uface_slave,L_Minus,L_Plus,doMPISides &
#if FV_ENABLED
                          ,pureDG,pureFV &
#endif
                          ,streamID &
                        )
! MODULES
USE ISO_C_BINDING,          ONLY: C_NULL_CHAR
USE MOD_Mesh_Vars,          ONLY: nElems,SideToElem, S2V2, firstMPISide_YOUR, lastMPISide_MINE, nSides
USE MOD_Device,             ONLY: STREAM_DEFAULT
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems
#endif
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_Interpolation_Vars, ONLY: d_L_Minus,d_L_Plus
USE MOD_Mesh_Vars,          ONLY: d_S2V2, d_SideToElem
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: TP_nVar
INTEGER,INTENT(IN)              :: Nloc
REAL,INTENT(IN)                 :: Uvol(TP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems)
REAL,INTENT(INOUT)              :: Uface_master(TP_nVar,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(INOUT)              :: Uface_slave( TP_nVar,0:Nloc,0:Nloc,1:nSides)
INTEGER(C_INT)                  :: d_Uvol          !< Key for device copy of Uvol array
INTEGER(C_INT)                  :: d_Uface_master  !< Key for device copy of Uface_master array
INTEGER(C_INT)                  :: d_Uface_slave   !< Key for device copy of Uface_slave array
REAL,INTENT(IN)                 :: L_Minus(0:Nloc),L_Plus(0:Nloc)
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE
#if FV_ENABLED
LOGICAL,INTENT(IN),OPTIONAL     :: pureDG      != .TRUE. prolongates all elements as DG elements
LOGICAL,INTENT(IN),OPTIONAL     :: pureFV      != .TRUE. prolongates all elements as FV elements
#endif
INTEGER,INTENT(IN),OPTIONAL     :: streamID    !< ID of device stream to run the device kernel in
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: firstSideID,lastSideID
INTEGER                         :: mystream
#if FV_ENABLED
INTEGER                         :: isFV(nElems)
INTEGER(C_INT)                  :: d_isFV 
#endif
!==================================================================================================================================

  mystream = STREAM_DEFAULT
  IF(PRESENT(streamID)) mystream = streamID

  ! Determining cell type and side ranges is backend independent and easier if we handle it here and pass the results to backends
#if FV_ENABLED
  isFV(:) = 0
  IF (.NOT.PRESENT(pureDG).AND..NOT.PRESENT(pureFV)) THEN
    isFV(:) = FV_Elems(:)
  ELSEIF (PRESENT(pureFV)) THEN
    IF(pureFV) isFV(:) = 1
  ELSE
    IF(.NOT.pureDG) isFV(:) = FV_Elems(:)
  END IF
#if (USE_ACCEL != ACCEL_OFF)
  CALL AllocateDeviceMemory(d_isFV, SIZE_C_INT, nElems)
  CALL CopyToDevice(d_isFV, C_Loc(isFV), nElems)
#endif
#endif /* FV_ENABLED */ 

  IF(doMPISides)THEN
    firstSideID = firstMPISide_YOUR
    lastSideID = nSides
  ELSE
    firstSideID = 1
    lastSideID = lastMPISide_MINE
  END IF

#if (USE_ACCEL == ACCEL_OFF)

  CALL ProlongToFace_Host(TP_nVar,Nloc,nElems,nSides,firstSideID,lastSideID, &
                          Uvol,Uface_master,Uface_slave,L_Minus,L_Plus,SideToElem,S2V2 &
#if FV_ENABLED
                          ,isFV &
#endif
                         )
#else /* USE_ACCEL */
  CALL ProlongToFace_Device(TP_nVar,Nloc,nElems,nSides,firstSideID,lastSideID, &
                            d_Uvol,d_Uface_master,d_Uface_slave,d_L_Minus,d_L_Plus,d_SideToElem,d_S2V2 &
#if FV_ENABLED
                            ,d_isFV &
#endif
                            ,mystream &
                           )
#endif /* USE_ACCEL */

END SUBROUTINE ProlongToFace

!==================================================================================================================================
!> Prolongs the conservatives variables to the sides
!> In the case of Gauss disc2, we project the entropy variables and then transform back to conservative variables
!> This is method isn't an entry point in the strictest definition, all it does it call the entry points for other kernels
!==================================================================================================================================
SUBROUTINE ProlongToFaceEntropy(TP_nVar,Nloc,Vvol,Vface_master,Vface_slave,d_Vvol,d_Vface_Master,d_Vface_slave, &
                                        Uface_master,Uface_slave,d_Uface_master,d_Uface_slave,L_Minus,L_Plus,doMPISides)
! MODULES
USE MOD_EOS,                ONLY: EntropyToCons
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR, lastMPISide_MINE, nSides
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: TP_nVar
INTEGER,INTENT(IN)              :: Nloc
REAL,INTENT(IN)                 :: Vvol(TP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems)
REAL,INTENT(INOUT)              :: Vface_master(TP_nVar,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(INOUT)              :: Vface_slave( TP_nVar,0:Nloc,0:Nloc,1:nSides)
INTEGER(C_INT),INTENT(IN)       :: d_Vvol
INTEGER(C_INT),INTENT(IN)       :: d_Vface_master
INTEGER(C_INT),INTENT(IN)       :: d_Vface_slave
REAL,INTENT(INOUT)              :: Uface_master(TP_nVar,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(INOUT)              :: Uface_slave( TP_nVar,0:Nloc,0:Nloc,1:nSides)
INTEGER(C_INT),INTENT(IN)       :: d_Uface_master
INTEGER(C_INT),INTENT(IN)       :: d_Uface_slave
REAL,INTENT(IN)                 :: L_Minus(0:Nloc),L_Plus(0:Nloc)
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: firstSideID,lastSideID
INTEGER                         :: ElemID,nbElemID,SideID
!==================================================================================================================================

  ! Prolong the entropy variables
  CALL ProlongToFace(TP_nVar,Nloc,Vvol,Vface_master,Vface_slave,d_Vvol,d_Vface_Master,d_Vface_slave,L_Minus,L_Plus,doMPISides)

  ! Transform back to conservative variables
  IF(doMPISides)THEN
    firstSideID = firstMPISide_YOUR
    lastSideID = nSides
  ELSE
    firstSideID = 1
    lastSideID = lastMPISide_MINE
  END IF

  DO SideID=firstSideID,lastSideID
    ElemID    = SideToElem(S2E_ELEM_ID,SideID)
    IF(ElemID.GT.0)   CALL EntropyToCons(Nloc,SideID,1,Vface_master(:,:,:,SideID),Uface_master(:,:,:,SideID),d_Vface_master,d_Uface_master)

    nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID)
    IF(nbElemID.GT.0) CALL EntropyToCons(Nloc,SideID,1,Vface_slave(:,:,:,SideID) ,Uface_slave(:,:,:,SideID),d_Vface_slave,d_Uface_slave)
  END DO

END SUBROUTINE ProlongToFaceEntropy


#if (USE_ACCEL == ACCEL_OFF)
!----------------------------------------------------------------------------------------------------------------------------------
! HOST BACKEND METHODS
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
!> Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
!> integration points, using fast 1D Interpolation and store in global side structure
!==================================================================================================================================
PPURE SUBROUTINE ProlongToFace_Host(TP_nVar,Nloc,nElems,nSides,firstSideID,lastSideID, &
                                    Uvol,Uface_master,Uface_slave,L_Minus,L_Plus,SideToElem,S2V2 &
#if FV_ENABLED
                                    ,isFV &
#endif
                                    )
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: TP_nVar
INTEGER,INTENT(IN)              :: Nloc
INTEGER,INTENT(IN)              :: nElems
INTEGER,INTENT(IN)              :: nSides
INTEGER,INTENT(IN)              :: firstSideID
INTEGER,INTENT(IN)              :: lastSideID
REAL,INTENT(IN)                 :: Uvol(TP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems)
REAL,INTENT(INOUT)              :: Uface_master(TP_nVar,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(INOUT)              :: Uface_slave( TP_nVar,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(IN)                 :: L_Minus(0:Nloc),L_Plus(0:Nloc)
INTEGER,INTENT(IN)              :: SideToElem(5,nSides)
INTEGER,INTENT(IN)              :: S2V2(2,0:Nloc,0:Nloc,0:4,1:6)
#if FV_ENABLED
INTEGER,INTENT(IN)              :: isFV(nElems)
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: p,q
INTEGER                         :: ElemID,nbElemID,locSide,nblocSide,SideID,flip
REAL                            :: Uface(TP_nVar,0:Nloc,0:Nloc)
!==================================================================================================================================

  DO SideID=firstSideID,lastSideID
    ElemID    = SideToElem(S2E_ELEM_ID,SideID)
    nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID)

    !master sides
    IF(ElemID.GT.0)THEN
      locSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
      flip    = 0

#if FV_ENABLED
      IF(PP_NodeType.EQ.1 .AND. isFV(ElemID).EQ.0)THEN
#else
      IF(PP_NodeType.EQ.1)THEN
#endif
        CALL EvalElemFaceG_Host(TP_nVar,Nloc,UVol(:,:,:,:,ElemID),Uface,L_Minus,L_Plus,locSide)
      ELSE
        CALL EvalElemFaceGL_Host(TP_nVar,Nloc,UVol(:,:,:,:,ElemID),Uface,locSide)
      END IF

      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        Uface_master(:,p,q,SideID)=Uface(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
      END DO; END DO
    END IF

    !slave side (ElemID,locSide and flip =-1 if not existing)
    IF(nbElemID.GT.0)THEN
      nblocSide = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
      flip      = SideToElem(S2E_FLIP,SideID)

#if FV_ENABLED
      IF(PP_NodeType.EQ.1 .AND. isFV(nbElemID).EQ.0)THEN
#else
      IF(PP_NodeType.EQ.1)THEN
#endif
        CALL EvalElemFaceG_Host(TP_nVar,Nloc,UVol(:,:,:,:,nbElemID),Uface,L_Minus,L_Plus,nblocSide)
      ELSE
        CALL EvalElemFaceGL_Host(TP_nVar,Nloc,UVol(:,:,:,:,nbElemID),Uface,nblocSide)
      END IF

      DO q=0,Nloc; DO p=0,Nloc
        Uface_slave( :,p,q,SideID)=Uface(:,S2V2(1,p,q,flip,nblocSide),S2V2(2,p,q,flip,nblocSide))
      END DO; END DO
    END IF
  END DO

END SUBROUTINE ProlongToFace_Host
#endif /* USE_ACCEL == ACCEL_OFF */


!==================================================================================================================================
!> Interpolates the element volume data stored at Gauss points
!> There is no entry point for this method. It and its device side counterpart are only ever called within other kernels
!> While this method is technically a HOST backend, it is called directly in at least one I/O routine that is called for all builds
!==================================================================================================================================
PPURE SUBROUTINE EvalElemFaceG_Host(TP_nVar,Nloc,Uvol,Uface,L_Minus,L_Plus,locSide)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: TP_nVar
INTEGER,INTENT(IN)              :: Nloc
INTEGER,INTENT(IN)              :: locSide
REAL,INTENT(IN)                 :: L_Minus(0:Nloc),L_Plus(0:Nloc)
REAL,INTENT(IN)                 :: Uvol( TP_nVar,0:Nloc,0:Nloc,0:Nloc)
REAL,INTENT(OUT)                :: Uface(TP_nVar,0:Nloc,0:Nloc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: l
!==================================================================================================================================

  SELECT CASE(locSide)
  CASE(XI_MINUS)
    Uface=Uvol(:,0,:,:)*L_Minus(0)
    DO l=1,Nloc
      Uface=Uface+Uvol(:,l,:,:)*L_Minus(l)
    END DO ! l
  CASE(ETA_MINUS)
    Uface=Uvol(:,:,0,:)*L_Minus(0)
    DO l=1,Nloc
      Uface=Uface+Uvol(:,:,l,:)*L_Minus(l)
    END DO ! l
  CASE(ZETA_MINUS)
    Uface=Uvol(:,:,:,0)*L_Minus(0)
    DO l=1,Nloc
      Uface=Uface+Uvol(:,:,:,l)*L_Minus(l)
    END DO ! l
  CASE(XI_PLUS)
    Uface=Uvol(:,0,:,:)*L_Plus(0)
    DO l=1,Nloc
      Uface=Uface+Uvol(:,l,:,:)*L_Plus(l)
    END DO ! l
  CASE(ETA_PLUS)
    Uface=Uvol(:,:,0,:)*L_Plus(0)
    DO l=1,Nloc
      Uface=Uface+Uvol(:,:,l,:)*L_Plus(l)
    END DO ! l
  CASE(ZETA_PLUS)
    Uface=Uvol(:,:,:,0)*L_Plus(0)
    DO l=1,Nloc
      Uface=Uface+Uvol(:,:,:,l)*L_Plus(l)
    END DO ! l
  END SELECT

END SUBROUTINE EvalElemFaceG_Host


!==================================================================================================================================
!> Interpolates the element volume data stored at Gauss-Lobatto points.
!> There is no entry point for this method. It and its device side counterpart are only ever called within other kernels
!> While this method is technically a HOST backend, it is called directly in at least one I/O routine that is called for all builds
!==================================================================================================================================
PPURE SUBROUTINE EvalElemFaceGL_Host(TP_nVar,Nloc,Uvol,Uface,locSide)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: TP_nVar
INTEGER,INTENT(IN)              :: Nloc
INTEGER,INTENT(IN)              :: locSide
REAL,INTENT(IN)                 :: Uvol( TP_nVar,0:Nloc,0:Nloc,0:Nloc)
REAL,INTENT(OUT)                :: Uface(TP_nVar,0:Nloc,0:Nloc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
  
  SELECT CASE(locSide)
  CASE(XI_MINUS)
    Uface=Uvol(:,0,:,:)
  CASE(ETA_MINUS)
    Uface=Uvol(:,:,0,:)
  CASE(ZETA_MINUS)
    Uface=Uvol(:,:,:,0)
  CASE(XI_PLUS)
    Uface=Uvol(:,Nloc,:,:)
  CASE(ETA_PLUS)
    Uface=Uvol(:,:,Nloc,:)
  CASE(ZETA_PLUS)
    Uface=Uvol(:,:,:,Nloc)
  END SELECT

END SUBROUTINE EvalElemFaceGL_Host

END MODULE MOD_ProlongToFace

