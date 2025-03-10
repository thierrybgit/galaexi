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

!==================================================================================================================================
!> \brief Prepares the integrand for the flux integrals, i.e. computes the common inteface fluxes.
!> This module prepares the computation of the fluxes over the sides by filling the two parts (advective and viscous)
!> of the common interface flux by calling the associated numerical flux functions.
!> We distinguish between inner sides, sides with boundary conditions and mpi sides.
!==================================================================================================================================
MODULE MOD_FillFlux
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE FillFlux
  MODULE PROCEDURE FillFlux
END INTERFACE

#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
  SUBROUTINE FillFlux_Device(Nloc,nSides,t,firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE,firstBCSide,nBCSides,doMPISides, &
                              d_Flux_master,d_Flux_slave,d_U_master,d_U_slave,d_UPrim_master,d_UPrim_slave, &
                              d_NormVec,d_TangVec1,d_TangVec2,d_SurfElem,d_Face_xGP &
                              ,streamID &
                            ) BIND(C, NAME="FillFlux_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: Nloc
    INTEGER(C_INT),VALUE :: nSides
    REAL(C_DOUBLE),VALUE :: t
    INTEGER(C_INT),VALUE :: firstInnerSide
    INTEGER(C_INT),VALUE :: lastInnerSide
    INTEGER(C_INT),VALUE :: firstMPISide_MINE
    INTEGER(C_INT),VALUE :: lastMPISide_MINE
    INTEGER(C_INT),VALUE :: firstBCSide
    INTEGER(C_INT),VALUE :: nBCSides
    LOGICAL,VALUE        :: doMPISides
    INTEGER(C_INT),VALUE :: d_Flux_master
    INTEGER(C_INT),VALUE :: d_Flux_slave
    INTEGER(C_INT),VALUE :: d_U_master
    INTEGER(C_INT),VALUE :: d_U_slave
    INTEGER(C_INT),VALUE :: d_UPrim_master
    INTEGER(C_INT),VALUE :: d_UPrim_slave
    INTEGER(C_INT),VALUE :: d_NormVec
    INTEGER(C_INT),VALUE :: d_TangVec1
    INTEGER(C_INT),VALUE :: d_TangVec2
    INTEGER(C_INT),VALUE :: d_SurfElem
    INTEGER(C_INT),VALUE :: d_Face_xGP
    INTEGER(C_INT),VALUE :: streamID
  END SUBROUTINE FillFlux_Device
END INTERFACE
#endif

PUBLIC::FillFlux
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Entry point to switch between backends for FillFlux
!==================================================================================================================================
SUBROUTINE FillFlux(t,Flux_master,Flux_slave,U_master,U_slave,UPrim_master,UPrim_slave,doMPISides &
#if FV_ENABLED
                    ,pureFV &
#endif
                    ,streamID)
USE MOD_PreProc
USE MOD_Device
USE MOD_Globals
USE MOD_Mesh_Vars,       ONLY: NormVec, TangVec1, TangVec2, SurfElem, Face_xGP
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE
USE MOD_Mesh_Vars,       ONLY: nSides,firstBCSide,nBCSides
#if PARABOLIC
USE MOD_Lifting_Vars,    ONLY: gradUx_master,gradUy_master,gradUz_master,gradUx_slave,gradUy_slave,gradUz_slave
#endif /*PARABOLIC*/
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,   ONLY: muSGS_master,muSGS_slave
#endif
#if FV_ENABLED
USE MOD_FV_Vars,         ONLY: FV_Elems_master,FV_Elems_slave,FV_Elems_Sum,FV_sVdm
#endif
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_DG_Vars, ONLY: d_Flux_master, d_Flux_slave, d_U_master, d_U_slave, d_UPrim_master, d_UPrim_slave
USE MOD_Mesh_Vars, ONLY: d_NormVec, d_TangVec1, d_TangVec2, d_SurfElem, d_Face_xGP
#if PARABOLIC
#endif
#if EDDYVISCOSITY
#endif
#if FV_ENABLED
#endif
#endif /* USE_ACCEL */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  !< = .TRUE. only MINE (where the proc is master)  MPISides are filled, =.FALSE. InnerSides
REAL,INTENT(IN)    :: t           !< physical time required for BC state evaluation in case of time dependent BCs
REAL,INTENT(OUT)   :: Flux_master(1:PP_nVar,0:PP_N,0:PP_N,1:nSides)      !< sum of advection and diffusion fluxes across the boundary
REAL,INTENT(OUT)   :: Flux_slave (1:PP_nVar,0:PP_N,0:PP_N,1:nSides)      !< sum of advection and diffusion fluxes across the boundary
REAL,INTENT(INOUT) :: U_master(    PP_nVar,0:PP_N, 0:PP_N,1:nSides)      !< solution on master sides
REAL,INTENT(INOUT) :: U_slave(     PP_nVar,0:PP_N, 0:PP_N,1:nSides)      !< solution on slave sides
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N, 0:PP_N, 1:nSides) !< primitive solution on master sides
REAL,INTENT(IN)    :: UPrim_slave( PP_nVarPrim,0:PP_N, 0:PP_N, 1:nSides) !< primitive solution on slave sides
#if FV_ENABLED
LOGICAL,INTENT(IN),OPTIONAL :: pureFV      != .TRUE. prolongates all elements as FV elements
#endif
INTEGER,INTENT(IN),OPTIONAL :: streamID
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: mystream
!==================================================================================================================================

  mystream = STREAM_DEFAULT
  IF (PRESENT(streamID)) mystream = streamID

#if (USE_ACCEL == ACCEL_OFF)

  CALL FillFlux_Host(PP_N,nSides,t,firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE,firstBCSide,nBCSides,doMPISides, &
                      Flux_master,Flux_slave,U_master,U_slave,UPrim_master,UPrim_slave, &
                      NormVec,TangVec1,TangVec2,SurfElem,Face_xGP &
#if PARABOLIC
                      ,gradUx_master,gradUy_master,gradUz_master,gradUx_slave,gradUy_slave,gradUz_slave &
#endif
#if EDDYVISCOSITY
                      ,muSGS_master,muSGS_slave &
#endif
#if FV_ENABLED
                      ,FV_Elems_master,FV_Elems_slave,FV_Elems_Sum,FV_sVdm,pureFV &
#endif
                    )

#else /* USE_ACCEL */

  Call FillFlux_Device(PP_N,nSides,t,firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE,firstBCSide,nBCSides,doMPISides, &
                        d_Flux_master,d_Flux_slave,d_U_master,d_U_slave,d_UPrim_master,d_UPrim_slave, &
                        d_NormVec,d_TangVec1,d_TangVec2,d_SurfElem,d_Face_xGP &
                        ,mystream &
                      )

#endif /* USE_ACCEL */

END SUBROUTINE FillFlux


! #if (USE_ACCEL == ACCEL_OFF)
!==================================================================================================================================
!> HOST BACKEND. Computes the fluxes for inner sides, MPI sides where the local proc is "master"  and boundary conditions.
!> The flux computation is performed separately for advection and diffusion fluxes in case
!> parabolic terms are considered.
!==================================================================================================================================
SUBROUTINE FillFlux_Host(Nloc,nSides,t,firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE,firstBCSide,nBCSides,doMPISides, &
                         Flux_master,Flux_slave,U_master,U_slave,UPrim_master,UPrim_slave, &
                         NormVec,TangVec1,TangVec2,SurfElem,Face_xGP &
#if PARABOLIC
                         ,gradUx_master,gradUy_master,gradUz_master,gradUx_slave,gradUy_slave,gradUz_slave &
#endif
#if EDDYVISCOSITY
                         ,muSGS_master,muSGS_slave &
#endif
#if FV_ENABLED
                         ,FV_Elems_master,FV_Elems_slave,FV_Elems_Sum,FV_sVdm,pureFV &
#endif
                        )
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_ChangeBasisByDim,ONLY: ChangeBasisSurf
USE MOD_Riemann,         ONLY: Riemann
USE MOD_GetBoundaryFlux, ONLY: GetBoundaryFlux
#if PARABOLIC
USE MOD_Riemann,         ONLY: ViscousFlux
#endif /*PARABOLIC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                !< Local polynomial order
INTEGER,INTENT(IN) :: nSides              !< Number of sides in the domain
REAL,INTENT(IN)    :: t                   !< physical time required for BC state evaluation in case of time dependent BCs
INTEGER,INTENT(IN) :: firstInnerSide      
INTEGER,INTENT(IN) :: lastInnerSide
INTEGER,INTENT(IN) :: firstMPISide_MINE
INTEGER,INTENT(IN) :: lastMPISide_MINE
INTEGER,INTENT(IN) :: firstBCSide
INTEGER,INTENT(IN) :: nBCSides
LOGICAL,INTENT(IN) :: doMPISides          !< = .TRUE. only MINE (where the proc is master)  MPISides are filled, =.FALSE. InnerSides
REAL,INTENT(OUT)   :: Flux_master(1:PP_nVar,0:Nloc,0:Nloc,1:nSides)      !< sum of advection and diffusion fluxes across the boundary
REAL,INTENT(OUT)   :: Flux_slave (1:PP_nVar,0:Nloc,0:Nloc,1:nSides)      !< sum of advection and diffusion fluxes across the boundary
REAL,INTENT(INOUT) :: U_master(    PP_nVar,0:Nloc, 0:Nloc,1:nSides)      !< solution on master sides
REAL,INTENT(INOUT) :: U_slave(     PP_nVar,0:Nloc, 0:Nloc,1:nSides)      !< solution on slave sides
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:Nloc, 0:Nloc, 1:nSides) !< primitive solution on master sides
REAL,INTENT(IN)    :: UPrim_slave( PP_nVarPrim,0:Nloc, 0:Nloc, 1:nSides) !< primitive solution on slave sides
REAL,INTENT(IN)    :: NormVec(3,0:Nloc,0:Nloc,1:nSides,0:FV_SIZE)
REAL,INTENT(IN)    :: TangVec1(3,0:Nloc,0:Nloc,1:nSides,0:FV_SIZE)
REAL,INTENT(IN)    :: TangVec2(3,0:Nloc,0:Nloc,1:nSides,0:FV_SIZE)
REAL,INTENT(IN)    :: SurfElem(0:Nloc,0:Nloc,1:nSides,0:FV_SIZE)
REAL,INTENT(IN)    :: Face_xGP(3,0:Nloc,0:Nloc,0:FV_SIZE,1:nSides)
#if PARABOLIC
REAL,INTENT(IN)    :: gradUx_master(PP_nVarLifting,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(IN)    :: gradUy_master(PP_nVarLifting,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(IN)    :: gradUz_master(PP_nVarLifting,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(IN)    :: gradUx_slave(PP_nVarLifting,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(IN)    :: gradUy_slave(PP_nVarLifting,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(IN)    :: gradUz_slave(PP_nVarLifting,0:Nloc,0:Nloc,1:nSides)
#endif
#if EDDYVISCOSITY
REAL,INTENT(IN)    :: muSGS_master(1,0:Nloc,0:Nloc,nSides)
REAL,INTENT(IN)    :: muSGS_slave(1,0:Nloc,0:Nloc,nSides)
#endif
#if FV_ENABLED
INTEGER,INTENT(IN) :: FV_Elems_master(nSides)
INTEGER,INTENT(IN) :: FV_Elems_slave(nSides)
INTEGER,INTENT(IN) :: FV_Elems_Sum(nSides)
INTEGER,INTENT(IN) :: FV_sVdm(0:Nloc,0:Nloc)
LOGICAL,INTENT(IN),OPTIONAL :: pureFV      != .TRUE. prolongates all elements as FV elements
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: SideID,p,q,firstSideID_wo_BC,firstSideID,lastSideID,FVEM
#if PARABOLIC
REAL    :: FluxV_loc(PP_nVar,0:Nloc, 0:Nloc)
#endif
INTEGER :: FV_Elems_Max(nSides) ! 0 if both sides DG, 1 else
!==================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver for advection and viscous terms
! Set the side range according to MPI or no MPI
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides (where the local proc is master)
  firstSideID_wo_BC = firstMPISide_MINE
  firstSideID = firstMPISide_MINE
   lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides that do not need communication
  firstSideID_wo_BC = firstInnerSide ! for fluxes
  firstSideID = firstBCSide    ! include BCs for master sides
   lastSideID = lastInnerSide
END IF

FV_Elems_Max = 0.0
#if FV_ENABLED
IF(PRESENT(pureFV))THEN
  IF(pureFV) FV_Elems_Max(firstSideID:lastSideID) = 1.
ELSE
  DO SideID=1,nBCSides
    FV_Elems_Max(SideID) = FV_Elems_master(SideID)
  END DO
  DO SideID=firstSideID_wo_BC,lastSideID
    FV_Elems_Max(SideID) = MAX(FV_Elems_master(SideID),FV_Elems_slave(SideID))
  END DO
END IF
#endif

! =============================
! Workflow:
!
!  1.  compute flux for non-BC sides
!  1.1) advective flux
!  1.2) viscous flux
!  1.3) add up viscous flux to Flux_master
!  2.  compute flux for BC sides
!  3.  multiply by SurfElem
!  4.  copy flux from Flux_master to Flux_slave
!  5.  convert FV flux to DG flux at mixed interfaces
!==============================

! 1. compute flux for non-BC sides
DO SideID=firstSideID_wo_BC,lastSideID
  ! 1.1) advective part of flux
  CALL Riemann(Nloc,Flux_master(:,:,:,SideID),&
      U_master    (:,:,:,SideID),U_slave    (:,:,:,SideID),       &
      UPrim_master(:,:,:,SideID),UPrim_slave(:,:,:,SideID),       &
      NormVec (:,:,:,SideID,FV_Elems_Max(SideID)), &
      TangVec1(:,:,:,SideID,FV_Elems_Max(SideID)), &
      TangVec2(:,:,:,SideID,FV_Elems_Max(SideID)),doBC=.FALSE.)

#if PARABOLIC
  ! 1.2) Fill viscous flux for non-BC sides
  CALL ViscousFlux(Nloc,FluxV_loc, UPrim_master(:,:,:,SideID), UPrim_slave  (:,:,:,SideID), &
      gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID), gradUz_master(:,:,:,SideID),&
      gradUx_slave (:,:,:,SideID),gradUy_slave (:,:,:,SideID), gradUz_slave (:,:,:,SideID),&
      NormVec(:,:,:,SideID,FV_Elems_Max(SideID))&
#if EDDYVISCOSITY
      ,muSGS_master(:,:,:,SideID),muSGS_slave(:,:,:,SideID)&
#endif
  )
  ! 1.3) add up viscous flux
  Flux_master(:,:,:,SideID) = Flux_master(:,:,:,SideID) + FluxV_loc(:,:,:)
#endif /*PARABOLIC*/
END DO ! SideID

! 2. Compute the fluxes at the boundary conditions: 1..nBCSides
IF(.NOT.doMPISides)THEN
  DO SideID=1,nBCSides
#if FV_ENABLED
    FVEM = FV_Elems_Max(SideID)
#endif
    CALL GetBoundaryFlux(SideID,t,Nloc,&
       Flux_master(  :,:,:,     SideID),&
       UPrim_master( :,:,:,     SideID),&
#if PARABOLIC
       gradUx_master(:,:,:,     SideID),&
       gradUy_master(:,:,:,     SideID),&
       gradUz_master(:,:,:,     SideID),&
#endif
       NormVec(      :,:,:,SideID,FVEM),&
       TangVec1(     :,:,:,SideID,FVEM),&
       TangVec2(     :,:,:,SideID,FVEM),&
       Face_xGP(     :,:,:,SideID,FVEM))
  END DO
END IF ! .NOT. MPISIDES


! 3. multiply by SurfElem
DO SideID=firstSideID,lastSideID
  ! multiply with SurfElem
  DO q=0,Nloc; DO p=0,Nloc
    Flux_master(:,p,q,SideID) = Flux_master(:,p,q,SideID) * SurfElem(p,q,SideID,FV_Elems_Max(SideID))
  END DO; END DO
END DO ! SideID

! 4. copy flux from master side to slave side
Flux_slave(:,:,:,firstSideID:lastSideID) = Flux_master(:,:,:,firstSideID:lastSideID)

#if FV_ENABLED == 1
! 5. convert flux on FV points to DG points for all DG faces at mixed interfaces
! only inner sides can be mixed (BC do not require a change basis)
DO SideID=firstSideID_wo_BC,lastSideID
  IF (FV_Elems_Sum(SideID).EQ.2) THEN
    CALL ChangeBasisSurf(PP_nVar,Nloc,Nloc,FV_sVdm,Flux_master(:,:,:,SideID))
  ELSE IF (FV_Elems_Sum(SideID).EQ.1) THEN
    CALL ChangeBasisSurf(PP_nVar,Nloc,Nloc,FV_sVdm,Flux_slave (:,:,:,SideID))
  END IF
END DO
#endif /* FV_ENABLED == 1 */

END SUBROUTINE FillFlux_Host
! #endif /* USE_ACCEL == ACCEL_OFF */


END MODULE MOD_FillFlux
