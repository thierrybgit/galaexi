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
!> Routines to map simple operations (copy,scalar multiply) from multi-dimensional arrays onto 1D arrays for better performance.
!==================================================================================================================================
MODULE MOD_Vector
! MODULES
USE ISO_C_BINDING
IMPLICIT NONE

PRIVATE

PUBLIC :: VCopy,VNullify,VSetConst
PUBLIC :: VAXPB_2STEP,VAXPB_IN,VAXPB_OUT
PUBLIC :: VAXPB_IN_OUT,VAXPB_ADD,VAXPB_CONST

!----------------------------------------------------------------------------------------------------------------------------------
! VAXPB Device Interfaces
!----------------------------------------------------------------------------------------------------------------------------------
#if (USE_ACCEL != ACCEL_OFF)

INTERFACE
  SUBROUTINE VAXPB_2STEP_Device(nTotal,d_VecOut,d_VecOut2,d_VecIn,Const1,Const2) &
    BIND(C, NAME="VAXPB_2STEP_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: nTotal
    INTEGER(C_INT),VALUE :: d_VecOut
    INTEGER(C_INT),VALUE :: d_VecOut2
    INTEGER(C_INT),VALUE :: d_VecIn
    REAL(C_DOUBLE),VALUE :: Const1
    REAL(C_DOUBLE),VALUE :: Const2
  END SUBROUTINE VAXPB_2STEP_Device
END INTERFACE


INTERFACE
  SUBROUTINE VAXPB_IN_OUT_Device(nTotal,d_VecOut,d_VecIn,ConstOut,ConstIn) &
    BIND(C, NAME="VAXPB_IN_OUT_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: nTotal
    INTEGER(C_INT),VALUE :: d_VecOut
    INTEGER(C_INT),VALUE :: d_VecIn
    REAL(C_DOUBLE),VALUE :: ConstOut
    REAL(C_DOUBLE),VALUE :: ConstIn
  END SUBROUTINE VAXPB_IN_OUT_Device
END INTERFACE


INTERFACE
  SUBROUTINE VAXPB_OUT_Device(nTotal,d_VecOut,d_VecIn,Const) &
    BIND(C, NAME="VAXPB_OUT_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: nTotal
    INTEGER(C_INT),VALUE :: d_VecOut
    INTEGER(C_INT),VALUE :: d_VecIn
    REAL(C_DOUBLE),VALUE :: Const
  END SUBROUTINE VAXPB_OUT_Device
END INTERFACE


INTERFACE
  SUBROUTINE VAXPB_IN_Device(nTotal,d_VecOut,d_VecIn,Const) &
    BIND(C, NAME="VAXPB_IN_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: nTotal
    INTEGER(C_INT),VALUE :: d_VecOut
    INTEGER(C_INT),VALUE :: d_VecIn
    REAL(C_DOUBLE),VALUE :: Const
  END SUBROUTINE VAXPB_IN_Device
END INTERFACE


INTERFACE
  SUBROUTINE VAXPB_ADD_Device(nTotal,d_VecOut,d_VecIn) &
    BIND(C, NAME="VAXPB_ADD_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: nTotal
    INTEGER(C_INT),VALUE :: d_VecOut
    INTEGER(C_INT),VALUE :: d_VecIn
  END SUBROUTINE VAXPB_ADD_Device
END INTERFACE


INTERFACE
  SUBROUTINE VAXPB_CONST_Device(nTotal,d_VecOut,Const) &
    BIND(C, NAME="VAXPB_CONST_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: nTotal
    INTEGER(C_INT),VALUE :: d_VecOut
    REAL(C_DOUBLE),VALUE :: Const
  END SUBROUTINE VAXPB_CONST_Device
END INTERFACE

#endif /* USE_ACCEL */

CONTAINS

!==================================================================================================================================
!> Y=0
!==================================================================================================================================
PURE SUBROUTINE VNullify(nTotal,Vec)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal                                               !< vector length
REAL,INTENT(OUT)      :: Vec(nTotal)                                          !< input vector
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
Vec=0.
END SUBROUTINE VNullify


!==================================================================================================================================
!> Y=a
!==================================================================================================================================
PURE SUBROUTINE VSetConst(nTotal,Vec,Const)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal                                               !< vector length
REAL,INTENT(OUT)      :: Vec(nTotal)                                          !< input vector
REAL,INTENT(IN)       :: Const                                                !< constant to set
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
Vec=Const
END SUBROUTINE VSetConst


!==================================================================================================================================
!> Y=X
!==================================================================================================================================
PURE SUBROUTINE VCopy(nTotal,VecOut,VecIn)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal                                               !< vector length
REAL,INTENT(IN)       :: VecIn(nTotal)                                        !< input vector
REAL,INTENT(OUT)      :: VecOut(nTotal)                                       !< output vector
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
VecOut=VecIn
END SUBROUTINE VCopy


!----------------------------------------------------------------------------------------------------------------------------------
! VAXPB Entry Point Methods
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
!> Two fused operations:
!> Y=a*Y+X
!> Z=b*Y+Z
!==================================================================================================================================
SUBROUTINE VAXPB_2STEP(nTotal,VecOut,VecOut2,VecIn,d_Out,d_Out2,d_In,Const1,Const2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal            !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)    !< output vector
REAL,INTENT(INOUT)    :: VecOut2(nTotal)   !< second output vector for 2 step variant
REAL,INTENT(IN)       :: VecIn(nTotal)     !< input vector
INTEGER(C_INT)        :: d_Out,d_In,d_Out2 !< device variable keys
REAL,INTENT(IN)       :: Const1            !< constant to multiply with in 1st op
REAL,INTENT(IN)       :: Const2            !< constant to multiply with in 2nd op
!==================================================================================================================================

#if (USE_ACCEL == ACCEL_OFF)
  CALL VAXPB_2STEP_Host(nTotal,VecOut,VecOut2,VecIn,Const1,Const2)
#else
  CALL VAXPB_2STEP_Device(nTotal,d_Out,d_Out2,d_In,Const1,Const2)
#endif

END SUBROUTINE VAXPB_2STEP
  
!==================================================================================================================================
!> Y = Y+a*X
!==================================================================================================================================
SUBROUTINE VAXPB_IN(nTotal,VecOut,VecIn,d_Out,d_In,Const)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal           !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)   !< output vector
REAL,INTENT(IN)       :: VecIn(nTotal)    !< input vector
INTEGER(C_INT)        :: d_Out,d_In       !< device variable keys
REAL,INTENT(IN)       :: Const            !< constant to multiply with
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================

#if (USE_ACCEL == ACCEL_OFF)
  CALL VAXPB_IN_Host(nTotal,VecOut,VecIn,Const)
#else
  CALL VAXPB_IN_Device(nTotal,d_Out,d_In,Const)
#endif

END SUBROUTINE VAXPB_IN

!==================================================================================================================================
!> Y = a*Y+X
!==================================================================================================================================
SUBROUTINE VAXPB_OUT(nTotal,VecOut,VecIn,d_Out,d_In,Const)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal           !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)   !< output vector
REAL,INTENT(IN)       :: VecIn(nTotal)    !< input vector
INTEGER(C_INT)        :: d_Out,d_In       !< device variable keys
REAL,INTENT(IN)       :: Const            !< constant to multiply with
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================

#if (USE_ACCEL == ACCEL_OFF)
  CALL VAXPB_OUT_Host(nTotal,VecOut,VecIn,Const)
#else
  CALL VAXPB_OUT_Device(nTotal,d_Out,d_In,Const)
#endif

END SUBROUTINE VAXPB_OUT

!==================================================================================================================================
!> Y = a*Y+b*X
!==================================================================================================================================
SUBROUTINE VAXPB_IN_OUT(nTotal,VecOut,VecIn,d_Out,d_In,ConstOut,ConstIn)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal            !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)    !< output vector
REAL,INTENT(IN)       :: VecIn(nTotal)     !< input vector
INTEGER(C_INT)        :: d_Out,d_In        !< device variable keys
REAL,INTENT(IN)       :: ConstIn           !< constant to multiply with input vec
REAL,INTENT(IN)       :: ConstOut          !< constant to multiply with output vec
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================

#if (USE_ACCEL == ACCEL_OFF)
  CALL VAXPB_IN_OUT_Host(nTotal,VecOut,VecIn,ConstOut,ConstIn)
#else
  CALL VAXPB_IN_OUT_Device(nTotal,d_Out,d_In,ConstOut,ConstIn)
#endif

END SUBROUTINE VAXPB_IN_OUT

!==================================================================================================================================
!> Y = Y+X
!==================================================================================================================================
SUBROUTINE VAXPB_ADD(nTotal,VecOut,VecIn,d_Out,d_In)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal           !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)   !< output vector
REAL,INTENT(IN)       :: VecIn(nTotal)    !< input vector
INTEGER(C_INT)        :: d_Out,d_In       !< device variable keys
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================

#if (USE_ACCEL == ACCEL_OFF)
  CALL VAXPB_ADD_Host(nTotal,VecOut,VecIn)
#else
  CALL VAXPB_ADD_Device(nTotal,d_Out,d_In)
#endif

END SUBROUTINE VAXPB_ADD

!==================================================================================================================================
!> Y = a*Y
!==================================================================================================================================
SUBROUTINE VAXPB_CONST(nTotal,VecOut,d_Out,Const)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal           !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)   !< output vector
INTEGER(C_INT)        :: d_Out            !< device variable keys
REAL,INTENT(IN)       :: Const            !< constant to multiply with
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================

#if (USE_ACCEL == ACCEL_OFF)
  CALL VAXPB_CONST_Host(nTotal,VecOut,Const)
#else
  CALL VAXPB_CONST_Device(nTotal,d_Out,Const)
#endif

END SUBROUTINE VAXPB_CONST


!----------------------------------------------------------------------------------------------------------------------------------
! VAXPB Host Backend Methods
!----------------------------------------------------------------------------------------------------------------------------------
#if (USE_ACCEL == ACCEL_OFF)
SUBROUTINE VAXPB_2STEP_Host(nTotal,VecOut,VecOut2,VecIn,Const1,Const2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal           !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)   !< output vector
REAL,INTENT(INOUT)    :: VecOut2(nTotal)  !< Optional second output vector for 2 step variant
REAL,INTENT(IN)       :: VecIn(nTotal)    !< input vector
REAL,INTENT(IN)       :: Const1           !< constant to multiply with in 1st op
REAL,INTENT(IN)       :: Const2           !< constant to multiply with in 2nd op
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================

  DO i=1,nTotal
    VecOut(i)  = VecOut(i)*Const1 + VecIn(i)
    VecOut2(i) = VecOut2(i) + VecOut(i)*Const2
  END DO

END SUBROUTINE VAXPB_2STEP_Host


SUBROUTINE VAXPB_IN_Host(nTotal,VecOut,VecIn,Const)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal           !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)   !< output vector
REAL,INTENT(IN)       :: VecIn(nTotal)    !< input vector
REAL,INTENT(IN)       :: Const            !< constant to multiply with
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================

  DO i=1,nTotal
    VecOut(i) = VecOut(i) + VecIn(i)*Const
  END DO

END SUBROUTINE VAXPB_IN_Host


SUBROUTINE VAXPB_OUT_Host(nTotal,VecOut,VecIn,Const)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal           !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)   !< output vector
REAL,INTENT(IN)       :: VecIn(nTotal)    !< input vector
REAL,INTENT(IN)       :: Const            !< constant to multiply with
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================

  DO i=1,nTotal
    VecOut(i) = VecOut(i)*Const + VecIn(i)
  END DO

END SUBROUTINE VAXPB_OUT_Host


SUBROUTINE VAXPB_IN_OUT_Host(nTotal,VecOut,VecIn,ConstOut,ConstIn)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal           !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)   !< output vector
REAL,INTENT(IN)       :: VecIn(nTotal)    !< input vector
REAL,INTENT(IN)       :: ConstIn          !< constant to multiply with input vec
REAL,INTENT(IN)       :: ConstOut         !< constant to multiply with output vec
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================

  DO i=1,nTotal
    VecOut(i) = VecOut(i)*ConstOut + VecIn(i)*ConstIn
  END DO

END SUBROUTINE VAXPB_IN_OUT_Host


SUBROUTINE VAXPB_ADD_Host(nTotal,VecOut,VecIn)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal           !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)   !< output vector
REAL,INTENT(IN)       :: VecIn(nTotal)    !< input vector
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================

  DO i=1,nTotal
    VecOut(i)=VecOut(i)+VecIn(i)
  END DO

END SUBROUTINE VAXPB_ADD_Host


SUBROUTINE VAXPB_CONST_Host(nTotal,VecOut,Const)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal           !< vector length
REAL,INTENT(INOUT)    :: VecOut(nTotal)   !< output vector
REAL,INTENT(IN)       :: Const            !< constant to multiply with
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================

  DO i=1,nTotal
    VecOut(i)=VecOut(i)*Const
  END DO

END SUBROUTINE VAXPB_CONST_Host

#endif /* (USE_ACCEL == ACCEL_OFF) */

END MODULE MOD_Vector
