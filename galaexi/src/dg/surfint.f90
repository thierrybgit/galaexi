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
!> Contains the different Surface integral formulations
!> Computes the Surface integral for all faces using U and updates Ut
!> Computes only inner surface integrals!
!> Surface integrals are separated for each direction
!==================================================================================================================================
MODULE MOD_SurfInt
IMPLICIT NONE
PRIVATE

INTEGER,PARAMETER :: TP_nVar = 1

INTERFACE SurfInt
  MODULE PROCEDURE SurfInt
END INTERFACE

#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
    SUBROUTINE SurfInt_Device(&
#if FV_ENABLED
                              d_FV_Elems_master, d_FV_Elems_slave, &
#endif
                              Nloc, nSides, nElems, firstMPISide_YOUR, lastMPISide_MINE, doMPISides, streamID, &
                              d_Flux_master, d_Flux_slave, d_Ut, &
                              d_L_HatMinus, d_L_HatPlus, &
                              d_ElemToSide, d_SideToElem, d_S2V2 &
#if (PP_NodeType == 1 && defined(SPLIT_DG))
                              ,d_U, d_UPrim, &
                              d_U_master, d_UPrim_master, &
                              d_U_slave, d_UPrim_slave, &
                              d_L_minus, d_L_plus, &
                              d_S2V, d_Metrics_fTilde,  &
                              d_Metrics_gTilde, d_Metrics_hTilde, d_Ja_Face, d_Ja_slave &
#endif
                            ) BIND(C, NAME="SurfInt_Device")
    USE ISO_C_BINDING 
    IMPLICIT NONE
#if FV_ENABLED
      INTEGER(C_INT),VALUE :: d_FV_Elems_master
      INTEGER(C_INT),VALUE :: d_FV_Elems_slave
#endif
      INTEGER, VALUE :: Nloc
      INTEGER, VALUE :: nSides
      INTEGER, VALUE :: nElems
      INTEGER, VALUE :: firstMPISide_YOUR
      INTEGER, VALUE :: lastMPISide_MINE
      LOGICAL, VALUE :: doMPISides
      INTEGER, VALUE :: streamID
      INTEGER(C_INT),VALUE :: d_Flux_master
      INTEGER(C_INT),VALUE :: d_Flux_slave
      INTEGER(C_INT),VALUE :: d_Ut
      INTEGER(C_INT),VALUE :: d_L_HatMinus
      INTEGER(C_INT),VALUE :: d_L_HatPlus
      INTEGER(C_INT),VALUE :: d_ElemToSide
      INTEGER(C_INT),VALUE :: d_SideToElem
      INTEGER(C_INT),VALUE :: d_S2V2
#if (PP_NodeType == 1 && defined(SPLIT_DG))
      INTEGER(C_INT),VALUE :: d_U
      INTEGER(C_INT),VALUE :: d_UPrim
      INTEGER(C_INT),VALUE :: d_U_master
      INTEGER(C_INT),VALUE :: d_U_slave
      INTEGER(C_INT),VALUE :: d_UPrim_master
      INTEGER(C_INT),VALUE :: d_UPrim_slave
      INTEGER(C_INT),VALUE :: d_L_minus
      INTEGER(C_INT),VALUE :: d_L_plus
      INTEGER(C_INT),VALUE :: d_S2V
      INTEGER(C_INT),VALUE :: d_Metrics_fTilde
      INTEGER(C_INT),VALUE :: d_Metrics_gTilde
      INTEGER(C_INT),VALUE :: d_Metrics_hTilde
      INTEGER(C_INT),VALUE :: d_Ja_Face
      INTEGER(C_INT),VALUE :: d_Ja_slave
#endif
    END SUBROUTINE SurfInt_Device
END INTERFACE 
#endif /* (USE_ACCEL != ACCEL_OFF) */

INTERFACE DoSurfInt
  MODULE PROCEDURE DoSurfInt_Host
END INTERFACE

PUBLIC :: DoSurfInt
PUBLIC :: SurfInt

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfInt

!==================================================================================================================================
!> Contains the surface integral for conservative quantities
!==================================================================================================================================
MODULE MOD_SurfIntCons
IMPLICIT NONE
PRIVATE

INTEGER,PARAMETER :: TP_nVar = PP_nVar

INTERFACE SurfIntCons
  MODULE PROCEDURE SurfInt
#if ((FV_ENABLED==2) && (PP_NodeType==1))
  MODULE PROCEDURE SurfIntBlend
#endif
END INTERFACE

#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
    SUBROUTINE SurfInt_Device(&
#if FV_ENABLED
                              d_FV_Elems_master, d_FV_Elems_slave, &
#endif
                              Nloc, nSides, nElems, firstMPISide_YOUR, lastMPISide_MINE, doMPISides, streamID, &
                              d_Flux_master, d_Flux_slave, d_Ut, &
                              d_L_HatMinus, d_L_HatPlus, &
                              d_ElemToSide, d_SideToElem, d_S2V2 &
#if (PP_NodeType == 1 && defined(SPLIT_DG))
                              ,d_U, d_UPrim, &
                              d_U_master, d_UPrim_master, &
                              d_U_slave, d_UPrim_slave, &
                              d_L_minus, d_L_plus, &
                              d_S2V, d_Metrics_fTilde,  &
                              d_Metrics_gTilde, d_Metrics_hTilde, d_Ja_Face, d_Ja_slave &
#endif
                            ) BIND(C, NAME="SurfInt_Device")
    USE ISO_C_BINDING 
    IMPLICIT NONE
#if FV_ENABLED
      INTEGER(C_INT),VALUE :: d_FV_Elems_master
      INTEGER(C_INT),VALUE :: d_FV_Elems_slave
#endif
      INTEGER, VALUE :: Nloc
      INTEGER, VALUE :: nSides
      INTEGER, VALUE :: nElems
      INTEGER, VALUE :: firstMPISide_YOUR
      INTEGER, VALUE :: lastMPISide_MINE
      LOGICAL, VALUE :: doMPISides
      INTEGER, VALUE :: streamID
      INTEGER(C_INT),VALUE :: d_Flux_master
      INTEGER(C_INT),VALUE :: d_Flux_slave
      INTEGER(C_INT),VALUE :: d_Ut
      INTEGER(C_INT),VALUE :: d_L_HatMinus
      INTEGER(C_INT),VALUE :: d_L_HatPlus
      INTEGER(C_INT),VALUE :: d_ElemToSide
      INTEGER(C_INT),VALUE :: d_SideToElem
      INTEGER(C_INT),VALUE :: d_S2V2
#if (PP_NodeType == 1 && defined(SPLIT_DG))
      INTEGER(C_INT),VALUE :: d_U
      INTEGER(C_INT),VALUE :: d_UPrim
      INTEGER(C_INT),VALUE :: d_U_master
      INTEGER(C_INT),VALUE :: d_U_slave
      INTEGER(C_INT),VALUE :: d_UPrim_master
      INTEGER(C_INT),VALUE :: d_UPrim_slave
      INTEGER(C_INT),VALUE :: d_L_minus
      INTEGER(C_INT),VALUE :: d_L_plus
      INTEGER(C_INT),VALUE :: d_S2V
      INTEGER(C_INT),VALUE :: d_Metrics_fTilde
      INTEGER(C_INT),VALUE :: d_Metrics_gTilde
      INTEGER(C_INT),VALUE :: d_Metrics_hTilde
      INTEGER(C_INT),VALUE :: d_Ja_Face
      INTEGER(C_INT),VALUE :: d_Ja_slave
#endif
    END SUBROUTINE SurfInt_Device
END INTERFACE

! INTERFACE
! #if ((FV_ENABLED==2) && (PP_NodeType==1))
!     SUBROUTINE SurfIntBlend_Device() BIND(C, NAME="SurfIntBlend_Device")
!       USE ISO_C_BINDING
!       IMPLICIT NONE
!     END SUBROUTINE SurfIntBlend_Device
! #endif
! END INTERFACE
#endif /* (USE_ACCEL != ACCEL_OFF) */

INTERFACE DoSurfIntCons
  MODULE PROCEDURE DoSurfInt_Host
END INTERFACE


PUBLIC :: DoSurfIntCons
PUBLIC :: SurfIntCons

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfIntCons

!==================================================================================================================================
!> Contains the surface integral for primitive quantities
!==================================================================================================================================
MODULE MOD_SurfIntLifting
IMPLICIT NONE
PRIVATE

INTEGER,PARAMETER :: TP_nVar = PP_nVarLifting

INTERFACE SurfIntLifting
  MODULE PROCEDURE SurfInt
END INTERFACE

#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
    SUBROUTINE SurfInt_Device(&
#if FV_ENABLED
                              d_FV_Elems_master, d_FV_Elems_slave, &
#endif
                              Nloc, nSides, nElems, firstMPISide_YOUR, lastMPISide_MINE, doMPISides, streamID, &
                              d_Flux_master, d_Flux_slave, d_Ut, &
                              d_L_HatMinus, d_L_HatPlus, &
                              d_ElemToSide, d_SideToElem, d_S2V2 &
#if (PP_NodeType == 1 && defined(SPLIT_DG))
                              ,d_U, d_UPrim, &
                              d_U_master, d_UPrim_master, &
                              d_U_slave, d_UPrim_slave, &
                              d_L_minus, d_L_plus, &
                              d_S2V, d_Metrics_fTilde,  &
                              d_Metrics_gTilde, d_Metrics_hTilde, d_Ja_Face, d_Ja_slave &
#endif
                            ) BIND(C, NAME="SurfInt_Device")
    USE ISO_C_BINDING 
    IMPLICIT NONE
#if FV_ENABLED
      INTEGER(C_INT),VALUE :: d_FV_Elems_master
      INTEGER(C_INT),VALUE :: d_FV_Elems_slave
#endif
      INTEGER, VALUE :: Nloc
      INTEGER, VALUE :: nSides
      INTEGER, VALUE :: nElems
      INTEGER, VALUE :: firstMPISide_YOUR
      INTEGER, VALUE :: lastMPISide_MINE
      LOGICAL, VALUE :: doMPISides
      INTEGER, VALUE :: streamID
      INTEGER(C_INT),VALUE :: d_Flux_master
      INTEGER(C_INT),VALUE :: d_Flux_slave
      INTEGER(C_INT),VALUE :: d_Ut
      INTEGER(C_INT),VALUE :: d_L_HatMinus
      INTEGER(C_INT),VALUE :: d_L_HatPlus
      INTEGER(C_INT),VALUE :: d_ElemToSide
      INTEGER(C_INT),VALUE :: d_SideToElem
      INTEGER(C_INT),VALUE :: d_S2V2
#if (PP_NodeType == 1 && defined(SPLIT_DG))
      INTEGER(C_INT),VALUE :: d_U
      INTEGER(C_INT),VALUE :: d_UPrim
      INTEGER(C_INT),VALUE :: d_U_master
      INTEGER(C_INT),VALUE :: d_U_slave
      INTEGER(C_INT),VALUE :: d_UPrim_master
      INTEGER(C_INT),VALUE :: d_UPrim_slave
      INTEGER(C_INT),VALUE :: d_L_minus
      INTEGER(C_INT),VALUE :: d_L_plus
      INTEGER(C_INT),VALUE :: d_S2V
      INTEGER(C_INT),VALUE :: d_Metrics_fTilde
      INTEGER(C_INT),VALUE :: d_Metrics_gTilde
      INTEGER(C_INT),VALUE :: d_Metrics_hTilde
      INTEGER(C_INT),VALUE :: d_Ja_Face
      INTEGER(C_INT),VALUE :: d_Ja_slave
#endif
    END SUBROUTINE SurfInt_Device
END INTERFACE
#endif /* (USE_ACCEL != ACCEL_OFF) */

INTERFACE DoSurfIntLifting
  MODULE PROCEDURE DoSurfInt_Host
END INTERFACE

PUBLIC :: DoSurfIntLifting
PUBLIC :: SurfIntLifting

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfIntLifting

!==================================================================================================================================
!> Contains the surface integral for primitive quantities
!==================================================================================================================================
MODULE MOD_SurfIntPrim
IMPLICIT NONE
PRIVATE

INTEGER,PARAMETER :: TP_nVar = PP_nVarPrim

INTERFACE SurfIntPrim
  MODULE PROCEDURE SurfInt
END INTERFACE

#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
    SUBROUTINE SurfInt_Device(&
#if FV_ENABLED
                              d_FV_Elems_master, d_FV_Elems_slave, &
#endif
                              Nloc, nSides, nElems, firstMPISide_YOUR, lastMPISide_MINE, doMPISides, streamID, &
                              d_Flux_master, d_Flux_slave, d_Ut, &
                              d_L_HatMinus, d_L_HatPlus, &
                              d_ElemToSide, d_SideToElem, d_S2V2 &
#if (PP_NodeType == 1 && defined(SPLIT_DG))
                              ,d_U, d_UPrim, &
                              d_U_master, d_UPrim_master, &
                              d_U_slave, d_UPrim_slave, &
                              d_L_minus, d_L_plus, &
                              d_S2V, d_Metrics_fTilde,  &
                              d_Metrics_gTilde, d_Metrics_hTilde, d_Ja_Face, d_Ja_slave &
#endif
                            ) BIND(C, NAME="SurfInt_Device")
    USE ISO_C_BINDING 
    IMPLICIT NONE
#if FV_ENABLED
      INTEGER(C_INT),VALUE :: d_FV_Elems_master
      INTEGER(C_INT),VALUE :: d_FV_Elems_slave
#endif
      INTEGER, VALUE :: Nloc
      INTEGER, VALUE :: nSides
      INTEGER, VALUE :: nElems
      INTEGER, VALUE :: firstMPISide_YOUR
      INTEGER, VALUE :: lastMPISide_MINE
      LOGICAL, VALUE :: doMPISides
      INTEGER, VALUE :: streamID
      INTEGER(C_INT),VALUE :: d_Flux_master
      INTEGER(C_INT),VALUE :: d_Flux_slave
      INTEGER(C_INT),VALUE :: d_Ut
      INTEGER(C_INT),VALUE :: d_L_HatMinus
      INTEGER(C_INT),VALUE :: d_L_HatPlus
      INTEGER(C_INT),VALUE :: d_ElemToSide
      INTEGER(C_INT),VALUE :: d_SideToElem
      INTEGER(C_INT),VALUE :: d_S2V2
#if (PP_NodeType == 1 && defined(SPLIT_DG))
      INTEGER(C_INT),VALUE :: d_U
      INTEGER(C_INT),VALUE :: d_UPrim
      INTEGER(C_INT),VALUE :: d_U_master
      INTEGER(C_INT),VALUE :: d_U_slave
      INTEGER(C_INT),VALUE :: d_UPrim_master
      INTEGER(C_INT),VALUE :: d_UPrim_slave
      INTEGER(C_INT),VALUE :: d_L_minus
      INTEGER(C_INT),VALUE :: d_L_plus
      INTEGER(C_INT),VALUE :: d_S2V
      INTEGER(C_INT),VALUE :: d_Metrics_fTilde
      INTEGER(C_INT),VALUE :: d_Metrics_gTilde
      INTEGER(C_INT),VALUE :: d_Metrics_hTilde
      INTEGER(C_INT),VALUE :: d_Ja_Face
      INTEGER(C_INT),VALUE :: d_Ja_slave
#endif
    END SUBROUTINE SurfInt_Device
END INTERFACE
#endif /* (USE_ACCEL != ACCEL_OFF) */

INTERFACE DoSurfIntPrim
  MODULE PROCEDURE DoSurfInt_Host
END INTERFACE

PUBLIC :: DoSurfIntPrim
PUBLIC :: SurfIntPrim

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfIntPrim
