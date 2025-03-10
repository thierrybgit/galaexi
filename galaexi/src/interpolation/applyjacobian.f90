!==================================================================================================================================
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
!==================================================================================================================================
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> \brief Module containing routines that changes scalar solution representation between physical and reference coordinates
!==================================================================================================================================
MODULE MOD_ApplyJacobian
USE ISO_C_BINDING
IMPLICIT NONE
PRIVATE

PUBLIC :: ApplyJacobian

! Device backend interface
#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
   SUBROUTINE ApplyJacobian_Device(TP_nVar,Nloc,nElems,d_U,toPhysical,d_sJ,FVE &
#if FV_ENABLED
                                   ,FV_Elems,chooseFV &
#endif                                  
                                   ,streamID &
                                  ) BIND(C, NAME="ApplyJacobian_Device")
   USE ISO_C_BINDING
   IMPLICIT NONE
      INTEGER(C_INT),VALUE :: TP_nVar
      INTEGER(C_INT),VALUE :: Nloc
      INTEGER(C_INT),VALUE :: nElems
      INTEGER(C_INT),VALUE :: d_U
      LOGICAL,VALUE        :: toPhysical
      INTEGER(C_INT),VALUE :: d_sJ
      INTEGER(C_INT),VALUE :: FVE
#if FV_ENABLED
      INTEGER(C_INT)       :: FV_Elems
      LOGICAL,VALUE        :: chooseFV
#endif
      INTEGER(C_INT),VALUE :: streamID
   END SUBROUTINE ApplyJacobian_Device
END INTERFACE
#endif

CONTAINS

   !==================================================================================================================================
   !> Entry point method to choose backend
   !==================================================================================================================================
   SUBROUTINE ApplyJacobian(TP_nVar,Nloc,U,d_U,toPhysical,FVE,streamID)
   ! MODULES
   USE MOD_PreProc
   USE MOD_Mesh_Vars ,ONLY: sJ,nElems
#if FV_ENABLED
   USE MOD_FV_Vars   ,ONLY: FV_Elems
#endif
   USE MOD_Device
#if (USE_ACCEL != ACCEL_OFF)
   USE MOD_Mesh_Vars ,ONLY: d_sJ
#endif
   IMPLICIT NONE
   !----------------------------------------------------------------------------------------------------------------------------------
   ! INPUT/OUTPUT VARIABLES
   INTEGER,INTENT(IN)            :: TP_nVar                                   !< Number of variables in 1st dim of array
   INTEGER,INTENT(IN)            :: Nloc                                      !< Local polynomial order
   REAL,INTENT(INOUT)            :: U(TP_nVar,0:Nloc,0:Nloc,0:Nloc,nElems)    !< Input/Output: Solution to be transformed
   INTEGER(C_INT),INTENT(IN)     :: d_U                                       !< Key for device copy of U
   LOGICAL,INTENT(IN)            :: toPhysical                                !< Switch for physical<-->reference transformation
   INTEGER,INTENT(IN),OPTIONAL   :: FVE                                       !< FV element switch
   INTEGER,INTENT(IN),OPTIONAL   :: streamID                                  !< ID of device stream to run backend with
   !----------------------------------------------------------------------------------------------------------------------------------
   ! LOCAL VARIABLES
   INTEGER :: locFVE
   LOGICAL :: chooseFV
   INTEGER :: mystream
   !==================================================================================================================================

      locFVE = 0
      chooseFV = .FALSE.
      IF(PRESENT(FVE)) THEN
         locFVE = FVE
         chooseFV = .TRUE.
      END IF

      mystream = STREAM_DEFAULT
      IF(PRESENT(streamID)) mystream = streamID

#if (USE_ACCEL == ACCEL_OFF)
      CALL ApplyJacobian_Host(TP_nVar,Nloc,nElems,U,toPhysical,sJ,locFVE &
#if FV_ENABLED
                              ,FV_Elems,chooseFV &
#endif
                              )
#else /* USE_ACCEL */
      CALL ApplyJacobian_Device(TP_nVar,Nloc,nElems,d_U,toPhysical,d_sJ,locFVE &
#if FV_ENABLED
                                ,FV_Elems,chooseFV &
#endif                                  
                                ,mystream &
                               )
#endif /* USE_ACCEL */

   END SUBROUTINE ApplyJacobian


#if (USE_ACCEL == ACCEL_OFF)
   !==================================================================================================================================
   !> Convert solution between physical <-> reference space, input will be overwritten with transformed solution
   !==================================================================================================================================
   PPURE SUBROUTINE ApplyJacobian_Host(TP_nVar,Nloc,nElems,U,toPhysical,sJ,FVE &
#if FV_ENABLED
                                       ,FV_Elems,chooseFV &
#endif
                                       )
   ! MODULES
   IMPLICIT NONE
   !----------------------------------------------------------------------------------------------------------------------------------
   ! INPUT/OUTPUT VARIABLES
   INTEGER,INTENT(IN) :: TP_nVar
   INTEGER,INTENT(IN) :: Nloc
   INTEGER,INTENT(IN) :: nElems
   REAL,INTENT(INOUT) :: U(TP_nVar,0:Nloc,0:Nloc,0:Nloc,nElems)      !< Input/Output: Solution to be transformed
   LOGICAL,INTENT(IN) :: toPhysical                                  !< Switch for physical<-->reference transformation
   REAL,INTENT(IN)    :: sJ(0:Nloc,0:Nloc,0:Nloc,nElems,0:FV_SIZE)   !< Jacobian
   INTEGER,INTENT(IN) :: FVE                                         !< Switch to handle FV elements (0 if no FV)
#if FV_ENABLED
   INTEGER,INTENT(IN) :: FV_Elems(nElems)                            !< Indicators for DG/FV type of element
   LOGICAL,INTENT(IN) :: chooseFV                                    !< Flag for whether or not to choose elements based on type
#endif
   !----------------------------------------------------------------------------------------------------------------------------------
   ! LOCAL VARIABLES
   INTEGER :: i,j,k,iElem
   INTEGER :: FVidx
   !==================================================================================================================================
   
      ! IF FV isn't on, then we will always use FVidx = 0, so set that as the default
      FVidx = 0

      IF(toPhysical)THEN
         DO iElem=1,nElems
#if FV_ENABLED
            IF (chooseFV) THEN
               IF (FV_Elems(iElem) == FVE) THEN
                  FVidx = FV_Elems(iElem)
               ELSE
                  CYCLE
               ENDIF               
            ELSE
               FVidx = FV_Elems(iElem)
            END IF
#endif
            DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
               U(:,i,j,k,iElem)=U(:,i,j,k,iElem)*sJ(i,j,k,iElem,FVidx)
            END DO; END DO; END DO
         END DO
      ELSE
         DO iElem=1,nElems
#if FV_ENABLED
            IF (chooseFV) THEN
               IF (FV_Elems(iElem) == FVE) THEN
                  FVidx = FV_Elems(iElem)
               ELSE
                  CYCLE
               ENDIF               
            ELSE
               FVidx = FV_Elems(iElem)
            END IF
#endif
            DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
               U(:,i,j,k,iElem)=U(:,i,j,k,iElem)/sJ(i,j,k,iElem,FVidx)
            END DO; END DO; END DO
         END DO
      END IF
   
   END SUBROUTINE ApplyJacobian_Host
#endif

END MODULE MOD_ApplyJacobian

