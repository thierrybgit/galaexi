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
!> Multiples volume data within an DG element linewise with a metric
!==================================================================================================================================
MODULE MOD_ApplyDMatrix
IMPLICIT NONE
PRIVATE

! Interfaces to C++ device backends
!----------------------------------------------------------------------------------------------------------------------------------
#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
    SUBROUTINE ApplyDMatrix_Device(Nloc,nElems,d_Ut,d_F,d_G,d_H,d_D_Hat_T,streamID,doOverwrite) &
                            BIND(C, NAME="ApplyDMatrix_Device")
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(C_INT), VALUE :: Nloc
        INTEGER(C_INT), VALUE :: nElems
        INTEGER(C_INT),VALUE  :: d_Ut
        INTEGER(C_INT),VALUE  :: d_F
        INTEGER(C_INT),VALUE  :: d_G
        INTEGER(C_INT),VALUE  :: d_H
        INTEGER(C_INT),VALUE  :: d_D_Hat_T
        INTEGER(C_INT), VALUE :: streamID
        LOGICAL, VALUE        :: doOverwrite
    END SUBROUTINE ApplyDMatrix_Device
END INTERFACE
#endif

PUBLIC::ApplyDMatrix

CONTAINS

!==================================================================================================================================
!> Entry point method for ApplyDMatrix. Splits to appropriate backend based on build.
!==================================================================================================================================
SUBROUTINE ApplyDMatrix(Nloc,nElems, Ut, F, G, H, D_Hat_T, streamID, doOverwrite)
USE MOD_PreProc
USE MOD_Device, ONLY: STREAM_DEFAULT
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_DG_Vars, ONLY: d_Ut, d_F, d_G, d_H, d_D_Hat_T
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc     !< Local polynomial order
INTEGER,INTENT(IN) :: nElems   !< Number of elements
REAL,INTENT(INOUT) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< time update
REAL,INTENT(IN)    :: f( PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< flux in x
REAL,INTENT(IN)    :: g( PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< flux in y
REAL,INTENT(IN)    :: h( PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< flux in z
REAL,INTENT(IN)    :: D_Hat_T(   0:PP_N,0:PP_N)                  !< Dervative matrix
INTEGER,OPTIONAL,INTENT(IN) :: streamID                                   !< Index of the stream to use on the device
LOGICAL,OPTIONAL,INTENT(IN) :: doOverwrite                                !< Toggle whether to overwrite or sum into Ut array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: mystream
LOGICAL :: doOverwrite_In
!==================================================================================================================================

    mystream = STREAM_DEFAULT
    IF (PRESENT(streamID)) mystream = streamID

    doOverwrite_In = .TRUE.
    IF (PRESENT(doOverwrite)) doOverwrite_In = doOverwrite

#if (USE_ACCEL == ACCEL_OFF)
    CALL ApplyDMatrix_Host(Nloc, nElems, Ut, F, G, H, D_Hat_T, doOverwrite_In)
#else
    CALL ApplyDMatrix_Device(Nloc, nElems, d_Ut, d_F, d_G, d_H, d_D_Hat_T, mystream, doOverwrite_In)
#endif

END SUBROUTINE ApplyDMatrix


#if (USE_ACCEL == ACCEL_OFF)
!==================================================================================================================================
!> Host backend for applying the derivative matrix linewise to the fluxes and adding to time derivative.
!==================================================================================================================================
SUBROUTINE ApplyDMatrix_Host(Nloc,nElems,Ut,F,G,H,D_Hat_T,doOverwrite)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc     !< Element local polynomial order
INTEGER,INTENT(IN) :: nElems   !< Number of elements
REAL,INTENT(INOUT) :: Ut(PP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems)  !< time update
REAL,INTENT(IN)    :: F( PP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems)  !< flux in x
REAL,INTENT(IN)    :: G( PP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems)  !< flux in y
REAL,INTENT(IN)    :: H( PP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems)  !< flux in z
REAL,INTENT(IN)    :: D_Hat_T(   0:Nloc,0:Nloc)                  !< Dervative matrix
LOGICAL,INTENT(IN) :: doOverwrite                                !< Toggle whether to overwrite or sum into Ut array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER   :: i,j,k,l,iElem
REAL      :: Ut_loc(PP_nVar)
!==================================================================================================================================
    
    DO iElem=1,nElems
        DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
            Ut_loc(:) = 0.
            DO l=0,Nloc
                ! Update the time derivative with the spatial derivatives of the transformed fluxes
                Ut_loc(:) = Ut_loc(:) + D_Hat_T(l,i)*F(:,l,j,k,iElem) &
                                      + D_Hat_T(l,j)*G(:,i,l,k,iElem) &
                                      + D_Hat_T(l,k)*H(:,i,j,l,iElem)
            END DO ! l
            ! Write temporary array to global array
            IF (doOverwrite) THEN
                Ut(:,i,j,k,iElem) =                     Ut_loc(:)
            ELSE
                Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Ut_loc(:)
            END IF
        END DO;END DO;END DO ! i,j,k
    END DO ! iElem

END SUBROUTINE ApplyDMatrix_Host
#endif /* ACCEL_OFF */

END MODULE MOD_ApplyDMatrix