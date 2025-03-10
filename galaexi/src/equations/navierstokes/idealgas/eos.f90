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
!> Subroutines  that provide gas properties and conversion between primitive and conservative variables for a ideal gas.
!==================================================================================================================================
MODULE MOD_EOS
! MODULES
USE ISO_C_BINDING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitEOS
  MODULE PROCEDURE InitEOS
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim_Point
  MODULE PROCEDURE ConsToPrim_Side
  MODULE PROCEDURE ConsToPrim_Sides
  MODULE PROCEDURE ConsToPrim_Elem
  MODULE PROCEDURE ConsToPrim_Volume
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
  MODULE PROCEDURE PrimToCons_Side
  MODULE PROCEDURE PrimToCons_Elem
  MODULE PROCEDURE PrimToCons_Volume
END INTERFACE

INTERFACE ConsToEntropy
  MODULE PROCEDURE ConsToEntropy_Volume
END INTERFACE

INTERFACE EntropyToCons
  MODULE PROCEDURE EntropyToCons_Side
END INTERFACE

INTERFACE PRESSURE_RIEMANN
  MODULE PROCEDURE PRESSURE_RIEMANN
END INTERFACE

! Interfaces to C++ device backends
!----------------------------------------------------------------------------------------------------------------------------------
#if (USE_ACCEL != ACCEL_OFF)
! ConsToPrim
INTERFACE
    SUBROUTINE ConsToPrim_Volume_Device(primKey, consKey, Nloc, KappaM1, R, nElems, streamID) &
                      BIND(C, name="ConsToPrim_Volume_Device")
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT),VALUE   :: primKey
      INTEGER(C_INT),VALUE   :: consKey
      INTEGER, VALUE :: Nloc
      REAL, VALUE    :: KappaM1
      REAL, VALUE    :: R
      INTEGER, VALUE :: nElems
      INTEGER, VALUE :: streamID
    END SUBROUTINE ConsToPrim_Volume_Device
END INTERFACE

INTERFACE
  SUBROUTINE ConsToPrim_Elem_Device(primKey, consKey, Nloc, KappaM1, R, elemIdx) &
                      BIND(C, name="ConsToPrim_Elem_Device")
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(C_INT),VALUE   :: primKey
    INTEGER(C_INT),VALUE   :: consKey
    INTEGER, VALUE :: Nloc
    REAL, VALUE    :: KappaM1
    REAL, VALUE    :: R
    INTEGER, VALUE :: elemIdx
  END SUBROUTINE ConsToPrim_Elem_Device
END INTERFACE

INTERFACE
  SUBROUTINE ConsToPrim_Side_Device(primKey, consKey, Nloc, KappaM1, R, elemIdx, sideIdx) &
                      BIND(C, name="ConsToPrim_Side_Device")
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(C_INT),VALUE   :: primKey
    INTEGER(C_INT),VALUE   :: consKey
    INTEGER, VALUE :: Nloc
    REAL, VALUE    :: KappaM1
    REAL, VALUE    :: R
    INTEGER, VALUE :: elemIdx
    INTEGER, VALUE :: sideIdx
  END SUBROUTINE ConsToPrim_Side_Device
END INTERFACE

INTERFACE
  SUBROUTINE ConsToPrim_Sides_Device(primKey, consKey, Nloc, KappaM1, R, elemIdx, sideIdx, nSides, streamID) &
                      BIND(C, name="ConsToPrim_Sides_Device")
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(C_INT),VALUE   :: primKey
    INTEGER(C_INT),VALUE   :: consKey
    INTEGER, VALUE :: Nloc
    REAL, VALUE    :: KappaM1
    REAL, VALUE    :: R
    INTEGER, VALUE :: elemIdx
    INTEGER, VALUE :: sideIdx
    INTEGER, VALUE :: nSides
    INTEGER, VALUE :: streamID
  END SUBROUTINE ConsToPrim_Sides_Device
END INTERFACE

INTERFACE
  SUBROUTINE ConsToPrim_Point_Device(primKey, consKey, Nloc, KappaM1, R, elemIdx, i, j, k) &
                      BIND(C, name="ConsToPrim_Point_Device")
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(C_INT),VALUE   :: primKey
    INTEGER(C_INT),VALUE   :: consKey
    INTEGER, VALUE :: Nloc
    REAL, VALUE    :: KappaM1
    REAL, VALUE    :: R
    INTEGER, VALUE :: elemIdx
    INTEGER, VALUE :: i,j,k
  END SUBROUTINE ConsToPrim_Point_Device
END INTERFACE

! ConsToEntropy
INTERFACE
  SUBROUTINE ConsToEntropy_Volume_Device(Nloc, nElems, d_entropy, d_cons, kappa, kappaM1, sKappaM1, streamID) &
            BIND(C, NAME="ConsToEntropy_Volume_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: Nloc
    INTEGER(C_INT),VALUE :: nElems
    INTEGER(C_INT),VALUE         :: d_entropy
    INTEGER(C_INT),VALUE         :: d_cons
    REAL(C_DOUBLE),VALUE :: kappa
    REAL(C_DOUBLE),VALUE :: kappaM1
    REAL(C_DOUBLE),VALUE :: sKappaM1
    INTEGER(C_INT),VALUE :: streamID
  END SUBROUTINE ConsToEntropy_Volume_Device
END INTERFACE

! EntropyToCons
INTERFACE
  SUBROUTINE EntropyToCons_Side_Device(Nloc, sideID, elemID, d_entropy, d_cons, kappa, kappaM1, sKappaM1) &
            BIND(C, NAME="EntropyToCons_Side_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: Nloc
    INTEGER(C_INT),VALUE :: sideID
    INTEGER(C_INT),VALUE :: elemID
    INTEGER(C_INT),VALUE         :: d_entropy
    INTEGER(C_INT),VALUE         :: d_cons
    REAL(C_DOUBLE),VALUE :: kappa
    REAL(C_DOUBLE),VALUE :: kappaM1
    REAL(C_DOUBLE),VALUE :: sKappaM1
  END SUBROUTINE EntropyToCons_Side_Device
END INTERFACE

! Init method
INTERFACE
  SUBROUTINE InitEOS_Device(EOS_Vars) BIND(C, NAME="InitEOS_Device")
    IMPLICIT NONE
    REAL :: EOS_Vars(PP_nVarEOS)
  END SUBROUTINE InitEOS_Device
END INTERFACE
#endif /* USE_ACCEL != ACCEL_OFF */

PUBLIC::InitEos
PUBLIC::ConsToPrim
PUBLIC::PrimToCons
PUBLIC::ConsToEntropy
PUBLIC::EntropyToCons
PUBLIC::PRESSURE_RIEMANN
PUBLIC::DefineParametersEos
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for the used eos
!==================================================================================================================================
SUBROUTINE DefineParametersEos()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("Equation of State")
CALL prms%CreateLogicalOption('UseNonDimensionalEqn',"Set true to compute R and mu from bulk Mach Reynolds (nondimensional form.",&
                                                '.FALSE.')
CALL prms%CreateRealOption(     'BulkMach',     "Bulk Mach     (UseNonDimensionEqn=T)")
CALL prms%CreateRealOption(     'BulkReynolds', "Bulk Reynolds (UseNonDimensionEqn=T)")
CALL prms%CreateRealOption(     'kappa',        "Heat capacity ratio / isentropic exponent"        , '1.4')
CALL prms%CreateRealOption(     'R',            "Specific gas constant"                            , '287.058')
CALL prms%CreateRealOption(     'Pr',           "Prandtl number"                                   , '0.72')
CALL prms%CreateRealOption(     'mu0',          "Dynamic Viscosity"                                , '0.')
CALL prms%CreateRealOption(     'Ts',           "Sutherland's law for variable viscosity: Ts"      , '110.4')
CALL prms%CreateRealOption(     'Tref',         "Sutherland's law for variable viscosity: Tref "   , '273.15')
CALL prms%CreateRealOption(     'ExpoSuth',     "Sutherland's law for variable viscosity: Exponent", '1.5')

END SUBROUTINE DefineParametersEos

!==================================================================================================================================
!> Initialize variables needed by the ideal gas equation of state.
!==================================================================================================================================
SUBROUTINE InitEos()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1,KappaP1,cp,cv
USE MOD_EOS_Vars      ,ONLY: R,sKappaM1,sKappaP1
#if PARABOLIC
USE MOD_EOS_Vars      ,ONLY: mu0,Pr,KappaSpr
#if PP_VISC == 1
USE MOD_EOS_Vars      ,ONLY: Ts,cSuth
#endif
#if (PP_VISC == 1) || (PP_VISC == 2)
USE MOD_EOS_Vars      ,ONLY: Tref,ExpoSuth
#endif
#endif /*PARABOLIC*/
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_EOS_Vars       ,ONLY: EOS_Vars
#endif

! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: BulkMach,BulkReynolds
LOGICAL :: UseNonDimensionalEqn=.FALSE.
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT IDEAL GAS...'


UseNonDimensionalEqn=GETLOGICAL('UseNonDimensionalEqn')
IF(UseNonDimensionalEqn)THEN
  BulkMach     = GETREAL('BulkMach')
  BulkReynolds = GETREAL('BulkReynolds')
END IF

! Gas constants
Kappa    = GETREAL('kappa')
KappaM1  = Kappa-1.
sKappaM1 = 1./KappaM1
KappaP1  = Kappa+1.
sKappaP1 = 1./(KappaP1)
IF(.NOT. UseNonDimensionalEqn)THEN
  R=GETREAL('R')
ELSE
  R=1./(Kappa*BulkMach*BulkMach)
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                              R | Set to 1/(Kappa*BulkMach**2)=',R
END IF

cp=R*kappa/(kappa-1.)
cv=R/(kappa-1.)

#if PARABOLIC
Pr       =GETREAL('Pr')
KappasPr =Kappa/Pr

! Viscosity
#if   PP_VISC == 0
! Constant viscosity
IF(.NOT. UseNonDimensionalEqn)THEN
  mu0=GETREAL('mu0')
ELSE
  mu0=1./BulkReynolds
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                            mu0 | Set to 1/BulkReynolds=',mu0
END IF
#elif PP_VISC == 1
! mu-Sutherland
! Coefficients from White, F. M., Viscous fluid flow, McGraw-Hill, 2006
mu0     =GETREAL('mu0')
Ts      =GETREAL('Ts')
Tref    =1./GETREAL('Tref')
ExpoSuth=GETREAL('ExpoSuth')
Ts      =Ts*Tref
cSuth   =Ts**ExpoSuth*(1+Ts)/(2*Ts*Ts)
#elif PP_VISC == 2
! mu power-law
IF(.NOT. UseNonDimensionalEqn)THEN
  mu0=GETREAL('mu0')
  Tref    =GETREAL('Tref')
ELSE
  mu0=1./BulkReynolds
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                            mu0 | Set to 1/BulkReynolds=',mu0
  Tref=1.
  SWRITE(UNIT_stdOut,'(A)')       ' |                           Tref | Set to 1.'
END IF
ExpoSuth=GETREAL('ExpoSuth')
mu0     =mu0/Tref**ExpoSuth
#endif
#endif /*PARABOLIC*/

! Do any initialization on the device
#if (USE_ACCEL != ACCEL_OFF)
! Fill HOST EOS array
EOS_Vars(EOS_KAPPA    ) = Kappa
EOS_Vars(EOS_KAPPAM1) = KappaM1
EOS_Vars(EOS_SKAPPAM1) = sKappaM1
EOS_Vars(EOS_KAPPAP1) = KappaP1
EOS_Vars(EOS_SKAPPAP1) = sKappaP1
EOS_Vars(EOS_R        ) = R
#if PARABOLIC
EOS_Vars(EOS_PR       ) = Pr
EOS_Vars(EOS_KAPPASPR ) = KappasPr
EOS_Vars(EOS_MU0      ) = mu0
#if PP_VISC==1
EOS_Vars(EOS_TS       ) = Ts
EOS_Vars(EOS_CSUTH    ) = cSuth
#endif
#if (PP_VISC==1) || (PP_VISC==2)
EOS_Vars(EOS_TREF     ) = Tref
EOS_Vars(EOS_EXPOSUTH ) = ExpoSuth
#endif
#endif /* PARABOLIC */

CALL InitEOS_Device(EOS_Vars)
#endif /* USE_ACCEL != ACCEL_OFF */

SWRITE(UNIT_stdOut,'(A)')' INIT IDEAL-GAS DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitEos

!----------------------------------------------------------------------------------------------------------------------------------
! ConstToPrim ENTRY POINT METHODS
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
!> Transformation from conservative variables to primitive variables in the whole volume
!==================================================================================================================================
SUBROUTINE ConsToPrim_Volume(Nloc,prim,cons,primKey,consKey,streamID)
! MODULES
USE MOD_EOS_Vars,ONLY:KappaM1,R
USE MOD_Mesh_Vars,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)       :: Nloc                                                  !< local polynomial degree of solution representation
REAL,INTENT(OUT)         :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:Nloc,1:nElems)       !< vector of primitive variables
REAL,INTENT(IN)          :: cons(PP_nVar    ,0:Nloc,0:Nloc,0:Nloc,1:nElems)       !< vector of conservative variables
INTEGER(C_INT),VALUE,INTENT(IN)  :: primKey                                               !< Key for primitive variable array device copy
INTEGER(C_INT),VALUE,INTENT(IN)  :: consKey                                               !< Key for convervative variable array device copy
INTEGER, INTENT(IN)      :: streamID                                              !< ID for device to stream to use on this call

!==================================================================================================================================

! Choose which backend to call based on configuration
#if (USE_ACCEL == ACCEL_OFF)
  call ConsToPrim_Volume_Host(prim, cons, Nloc, KappaM1, R, nElems, streamID)
#else
  call ConsToPrim_Volume_Device(primKey, consKey, Nloc, KappaM1, R, nElems, streamID)
#endif

END SUBROUTINE ConsToPrim_Volume

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables for a single element
!==================================================================================================================================
SUBROUTINE ConsToPrim_Elem(Nloc,prim,cons,primKey,consKey,elemIdx)
! MODULES
USE MOD_EOS_Vars,ONLY:KappaM1,R
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)       :: Nloc                                          !< local polynomial degree of solution representation
REAL,INTENT(OUT)         :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:Nloc)        !< vector of primitive variables
REAL,INTENT(IN)          :: cons(PP_nVar    ,0:Nloc,0:Nloc,0:Nloc)        !< vector of conservative variables
INTEGER(C_INT),VALUE,INTENT(IN)  :: primKey                                       !< Key for primitive variable array device copy
INTEGER(C_INT),VALUE,INTENT(IN)  :: consKey                                       !< Key for convervative variable array device copy
INTEGER, INTENT(IN)      :: elemIdx                                       !< Index of the element that is being converted
!==================================================================================================================================

! Choose which backend to call based on configuration
#if (USE_ACCEL == ACCEL_OFF)
  call ConsToPrim_Elem_Host(prim, cons, Nloc, KappaM1, R, elemIdx)
#else
  call ConsToPrim_Elem_Device(primKey, consKey, Nloc, KappaM1, R, elemIdx)
#endif

END SUBROUTINE ConsToPrim_Elem

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables on a single side
!==================================================================================================================================
SUBROUTINE ConsToPrim_Side(Nloc,prim,cons,primKey,consKey,elemIdx,sideIdx)
! MODULES
USE MOD_EOS_Vars,ONLY:KappaM1,R
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                         !< local polynomial degree of solution representation
REAL,INTENT(OUT)         :: prim(PP_nVarPrim,0:Nloc,0:Nloc)        !< vector of primitive variables
REAL,INTENT(IN)          :: cons(PP_nVar    ,0:Nloc,0:Nloc)        !< vector of conservative variables
INTEGER(C_INT),VALUE,INTENT(IN)  :: primKey                                !< Key for primitive variable array device copy
INTEGER(C_INT),VALUE,INTENT(IN)  :: consKey                                !< Key for convervative variable array device copy
INTEGER, INTENT(IN)      :: elemIdx                                !< Index of the element that is being converted
INTEGER, INTENT(IN)      :: sideIdx                                !< Index of the desired side within the element

!==================================================================================================================================

! Choose which backend to call based on configuration
#if (USE_ACCEL == ACCEL_OFF)
  call ConsToPrim_Side_Host(prim, cons, Nloc, KappaM1, R, elemIdx, sideIdx)
#else
  call ConsToPrim_Side_Device(primKey, consKey, Nloc, KappaM1, R, elemIdx, sideIdx)
#endif

END SUBROUTINE ConsToPrim_Side

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables on a multiple contiguous sides
!==================================================================================================================================
SUBROUTINE ConsToPrim_Sides(Nloc,prim,cons,primKey,consKey,elemIdx,sideIdx,nSides,streamID)
! MODULES
USE MOD_EOS_Vars,ONLY:KappaM1,R
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)       :: Nloc                                   !< local polynomial degree of solution representation
REAL,INTENT(OUT)         :: prim(PP_nVarPrim,0:Nloc,0:Nloc,nSides) !< vector of primitive variables
REAL,INTENT(IN)          :: cons(PP_nVar    ,0:Nloc,0:Nloc,nSides) !< vector of conservative variables
INTEGER(C_INT),VALUE,INTENT(IN)  :: primKey                                !< Key for primitive variable array device copy
INTEGER(C_INT),VALUE,INTENT(IN)  :: consKey                                !< Key for convervative variable array device copy
INTEGER, INTENT(IN)      :: elemIdx                                !< Index of the element that is being converted (0 for nSides only arrays)
INTEGER, INTENT(IN)      :: sideIdx                                !< Index of the starting side
INTEGER, INTENT(IN)      :: nSides                                 !< Number of sides to convert within the element
INTEGER, INTENT(IN)      :: streamID                               !< Device stream to run backend on (accelerated only)
!==================================================================================================================================

! Choose which backend to call based on configuration
#if (USE_ACCEL == ACCEL_OFF)
  call ConsToPrim_Sides_Host(prim, cons, Nloc, KappaM1, R, elemIdx, sideIdx, nSides)
#else
  call ConsToPrim_Sides_Device(primKey, consKey, Nloc, KappaM1, R, elemIdx, sideIdx, nSides, streamID)
#endif

END SUBROUTINE ConsToPrim_Sides


!==================================================================================================================================
!> Transformation from conservative variables to primitive variables for a single DOF
!==================================================================================================================================
SUBROUTINE ConsToPrim_Point(Nloc,prim,cons,primKey,consKey,elemIdx,i,j,k,runOnHost)
! MODULES
USE MOD_EOS_Vars,ONLY:KappaM1,R
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                           !< local polynomial degree of solution representation
REAL,INTENT(OUT)         :: prim(PP_nVarPrim)        !< vector of primitive variables
REAL,INTENT(IN)          :: cons(PP_nVar    )        !< vector of conservative variables
INTEGER(C_INT),VALUE,INTENT(IN)  :: primKey                  !< Key for primitive variable array device copy
INTEGER(C_INT),VALUE,INTENT(IN)  :: consKey                  !< Key for convervative variable array device copy
INTEGER, INTENT(IN)      :: elemIdx                  !< Index of the element that is being converted
INTEGER, INTENT(IN)      :: i,j,k                    !< Indices of the desired DOF within the element
LOGICAL,OPTIONAL,INTENT(IN) :: runOnHost             !< If .TRUE., runs the host backend. Does nothing when accel is off.
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: runOnHostLoc
!==================================================================================================================================

! Choose which backend to call based on configuration
#if (USE_ACCEL == ACCEL_OFF)
  call ConsToPrim_Point_Host(prim, cons, Nloc, KappaM1, R, elemIdx, i, j, k)
#else
  runOnHostLoc = .FALSE.
  IF(PRESENT(runOnHost)) runOnHostLoc = runOnHost

  IF (runOnHostLoc) THEN
    call ConsToPrim_Point_Host(prim, cons, Nloc, KappaM1, R, elemIdx, i, j, k)
  ELSE
    call ConsToPrim_Point_Device(primKey, consKey, Nloc, KappaM1, R, elemIdx, i, j, k)
  END IF
#endif

END SUBROUTINE ConsToPrim_Point


!----------------------------------------------------------------------------------------------------------------------------------
! HOST ConsToPrim METHODS
!----------------------------------------------------------------------------------------------------------------------------------

!==================================================================================================================================
!> CPU backend for transformation for core variable transformation for ConsToPrim
!> See note on ConsToPrim_Point_Host below for why this is compile regardless of status
!==================================================================================================================================
SUBROUTINE ConsToPrimCore_Host(prim,cons,KappaM1,R)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT) :: prim(:) !< vector of primitive variables (density,velocities,temperature,pressure)
REAL,INTENT(IN)  :: cons(:) !< vector of conservative variables
REAL,INTENT(IN)  :: kappaM1 !< = \f$\kappa - 1\f$
REAL,INTENT(IN)  :: R
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho    ! 1/Rho
!==================================================================================================================================

  sRho=1.0/cons(DENS)
  ! density
  prim(DENS)=cons(DENS)
  ! velocity
  prim(VEL1:VEL3)=cons(MOM1:MOM3)*sRho
  ! pressure
  prim(PRES)=kappaM1*(cons(ENER)-0.5*SUM(cons(MOMV)*prim(VELV)))
  ! temperature
  prim(TEMP) = prim(PRES)*sRho / R

END SUBROUTINE ConsToPrimCore_Host

#if (USE_ACCEL == ACCEL_OFF)
!==================================================================================================================================
!> CPU backend of transformation from conservative variables to primitive variables in the whole volume
!==================================================================================================================================
SUBROUTINE ConsToPrim_Volume_Host(prim, cons, Nloc, KappaM1, R, nElems, streamID)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(OUT)   :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:Nloc,1:nElems)       !< vector of primitive variables
REAL, INTENT(IN)    :: cons(PP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems    )       !< vector of conservative variables
INTEGER, INTENT(IN) :: Nloc                  !< Local polynomial degree
REAL, INTENT(IN)    :: KappaM1               !< ratio of specific heats - 1.0
REAL, INTENT(IN)    :: R                     !< specific gas constant
INTEGER, INTENT(IN) :: nElems                !< Total elements in the proc local mesh
INTEGER, INTENT(IN) :: streamID              !< ID for device to stream to use on this call

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
!==================================================================================================================================

    DO iElem=1,nElems
        DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
            CALL ConsToPrimCore_Host(prim(:,i,j,k,iElem),cons(:,i,j,k,iElem),KappaM1,R)
        END DO; END DO; END DO! i,j,k=0,Nloc
    END DO ! iElem  

END SUBROUTINE ConsToPrim_Volume_Host

!==================================================================================================================================
!> CPU backend of transformation from conservative variables to primitive variables in the whole volume
!==================================================================================================================================
SUBROUTINE ConsToPrim_Elem_Host(prim, cons, Nloc, KappaM1, R, elemIdx)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(OUT)   :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:Nloc)       !< vector of primitive variables
REAL, INTENT(IN)    :: cons(PP_nVar    ,0:Nloc,0:Nloc,0:Nloc)       !< vector of conservative variables
INTEGER, INTENT(IN) :: Nloc                  !< Local polynomial degree
REAL, INTENT(IN)    :: KappaM1               !< ratio of specific heats - 1.0
REAL, INTENT(IN)    :: R                     !< specific gas constant
INTEGER, INTENT(IN) :: elemIdx               !< Index of the element that is being converted

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k
!==================================================================================================================================

  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    CALL ConsToPrimCore_Host(prim(:,i,j,k),cons(:,i,j,k),KappaM1,R)
  END DO; END DO; END DO! i,j,k=0,Nloc

END SUBROUTINE ConsToPrim_Elem_Host

!==================================================================================================================================
!> CPU backend of transformation from conservative variables to primitive variables on a single side
!==================================================================================================================================
SUBROUTINE ConsToPrim_Side_Host(prim, cons, Nloc, KappaM1, R, elemIdx, sideIdx)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(OUT)   :: prim(PP_nVarPrim,0:Nloc,0:Nloc)           !< vector of primitive variables
REAL, INTENT(IN)    :: cons(PP_nVar    ,0:Nloc,0:Nloc)           !< vector of conservative variables
INTEGER, INTENT(IN) :: Nloc                  !< Local polynomial degree
REAL, INTENT(IN)    :: KappaM1               !< ratio of specific heats - 1.0
REAL, INTENT(IN)    :: R                     !< specific gas constant
INTEGER, INTENT(IN) :: elemIdx               !< Index of the element that is being converted.  Zero for "sides only" arrays
INTEGER, INTENT(IN) :: sideIdx               !< Index of the desired side

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q
!==================================================================================================================================

  DO q=0,Nloc; DO p=0,Nloc
    CALL ConsToPrimCore_Host(prim(:,p,q),cons(:,p,q),KappaM1,R)
  END DO; END DO

END SUBROUTINE ConsToPrim_Side_Host

!==================================================================================================================================
!> CPU backend of transformation from conservative variables to primitive variables on a multiple contiguous sides
!==================================================================================================================================
SUBROUTINE ConsToPrim_Sides_Host(prim, cons, Nloc, KappaM1, R, elemIdx, sideIdx, nSides)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(OUT)   :: prim(PP_nVarPrim,0:Nloc,0:Nloc,nSides)         !< vector of primitive variables
REAL, INTENT(IN)    :: cons(PP_nVar    ,0:Nloc,0:Nloc,nSides)         !< vector of conservative variables
INTEGER, INTENT(IN) :: Nloc                  !< Local polynomial degree
REAL, INTENT(IN)    :: KappaM1               !< ratio of specific heats - 1.0
REAL, INTENT(IN)    :: R                     !< specific gas constant
INTEGER, INTENT(IN) :: elemIdx               !< Index of the element that is being converted.  Zero for "sides only" arrays
INTEGER, INTENT(IN) :: sideIdx               !< Index of the starting side
INTEGER, INTENT(IN) :: nSides                !< Number of sides to convert

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,iSide
!==================================================================================================================================

  DO iSide = 1,nSides
    DO q=0,Nloc; DO p=0,Nloc
      CALL ConsToPrimCore_Host(prim(:,p,q,iSide),cons(:,p,q,iSide),KappaM1,R)
    END DO; END DO
  END DO

END SUBROUTINE ConsToPrim_Sides_Host


#endif /* USE_ACCEL == ACCEL_OFF */

!==================================================================================================================================
!> CPU backend of transformation from conservative variables to primitive variables for a single DOF
!> We build this regardless of acceleration, because sometime we want to be able to this on the host
!> Regardless if acceleration is turned on, simply because it isn't that much work.
!==================================================================================================================================
SUBROUTINE ConsToPrim_Point_Host(prim, cons, Nloc, KappaM1, R, elemIdx, i, j, k)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(OUT)   :: prim(PP_nVarPrim)     !< vector of primitive variables
REAL, INTENT(IN)    :: cons(PP_nVar)         !< vector of conservative variables
INTEGER, INTENT(IN) :: Nloc                  !< Local polynomial degree
REAL, INTENT(IN)    :: KappaM1               !< ratio of specific heats - 1.0
REAL, INTENT(IN)    :: R                     !< specific gas constant
INTEGER, INTENT(IN) :: elemIdx               !< Index of the element that is being converted
INTEGER, INTENT(IN) :: i,j,k                 !< Indices of the desired DOF within the element

!==================================================================================================================================

  CALL ConsToPrimCore_Host(prim(:),cons(:),KappaM1,R)

END SUBROUTINE ConsToPrim_Point_Host


!----------------------------------------------------------------------------------------------------------------------------------
! PrimToCons Methods
!----------------------------------------------------------------------------------------------------------------------------------

!==================================================================================================================================
!> Transformation from primitive to conservative variables for a single state
!==================================================================================================================================
PPURE SUBROUTINE PrimToCons(prim,cons)
! MODULES
USE MOD_EOS_Vars,ONLY:sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: prim(PP_nVarPrim) !< vector of primitive variables
REAL,INTENT(OUT) :: cons(PP_nVar)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! density
cons(DENS)=prim(DENS)
! momentum
cons(MOM1:MOM2)=prim(VEL1:VEL2)*prim(DENS)
#if (PP_dim==3)
cons(MOM3)=prim(VEL3)*prim(DENS)
#else
cons(MOM3)=0.
#endif
! energy
cons(ENER)=sKappaM1*prim(PRES)+0.5*SUM(cons(MOMV)*prim(VELV))
END SUBROUTINE PrimToCons

!==================================================================================================================================
!> Transformation from primitive to conservative variables on a single side
!==================================================================================================================================
PPURE SUBROUTINE PrimToCons_Side(Nloc,prim,cons)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                  !< local polynomial degree of solution representation
REAL,INTENT(IN)    :: prim(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)) !< vector of primitive variables
REAL,INTENT(OUT)   :: cons(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q
!==================================================================================================================================
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  CALL PrimToCons(prim(:,p,q),cons(:,p,q))
END DO; END DO ! p,q=0,Nloc
END SUBROUTINE PrimToCons_Side

!==================================================================================================================================
!> Transformation from primitive to conservative variables in the whole volume
!==================================================================================================================================
PPURE SUBROUTINE PrimToCons_Elem(Nloc,prim,cons)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                             !< local polynomial degree of solution representation
REAL,INTENT(IN)    :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc))     !< vector of primitive variables
REAL,INTENT(OUT)   :: cons(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc))     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k
!==================================================================================================================================
DO k=0,ZDIM(Nloc); DO j=0,Nloc; DO i=0,Nloc
  CALL PrimToCons(prim(:,i,j,k),cons(:,i,j,k))
END DO; END DO; END DO
END SUBROUTINE PrimToCons_Elem

!==================================================================================================================================
!> Transformation from primitive to conservative variables in the whole volume
!==================================================================================================================================
PPURE SUBROUTINE PrimToCons_Volume(Nloc,prim,cons)
! MODULES
USE MOD_Mesh_Vars,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                                  !< local polynomial degree of solution representation
REAL,INTENT(IN)    :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)     !< vector of primitive variables
REAL,INTENT(OUT)   :: cons(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iElem
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,ZDIM(Nloc); DO j=0,Nloc; DO i=0,Nloc
    CALL PrimToCons(prim(:,i,j,k,iElem),cons(:,i,j,k,iElem))
  END DO; END DO; END DO
END DO
END SUBROUTINE PrimToCons_Volume


!----------------------------------------------------------------------------------------------------------------------------------
! ConsToEntropy Methods
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
!> Entry point for ConsToEntropy_Volume
!==================================================================================================================================
SUBROUTINE ConsToEntropy_Volume(Nloc,entropy,cons,d_entropy,d_cons,streamID)
! MODULES
USE MOD_Mesh_Vars, ONLY: nElems
USE MOD_EOS_Vars, ONLY: KappaM1,Kappa,sKappaM1
USE MOD_Device, ONLY: STREAM_DEFAULT
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: Nloc                                             !< Local polynomial order
REAL,INTENT(OUT)            :: entropy(PP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems)   !< vector of entropy variables
REAL,INTENT(IN)             :: cons(PP_nVar   ,0:Nloc,0:Nloc,0:Nloc,1:nElems)   !< vector of conservative variables
INTEGER(C_INT),VALUE,INTENT(IN)     :: d_entropy                                        !< Device variable key for entropy array
INTEGER(C_INT),VALUE,INTENT(IN)     :: d_cons                                           !< Device variable key for cons array
INTEGER,INTENT(IN),OPTIONAL :: streamID                                         !< Index for stream to run kernel on device with
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: mystream
!==================================================================================================================================

  mystream = STREAM_DEFAULT
  IF(PRESENT(streamID)) mystream = streamID

#if (USE_ACCEL == ACCEL_OFF)
  CALL ConsToEntropy_Volume_Host(Nloc,nElems,entropy,cons,Kappa,KappaM1,sKappaM1)
#else /* USE_ACCEL */
  CALL ConsToEntropy_Volume_Device(Nloc,nElems,d_entropy,d_cons,Kappa,KappaM1,sKappaM1,mystream)
#endif /* USE_ACCEL */

END SUBROUTINE ConsToEntropy_Volume

!==================================================================================================================================
!> Transformation from conservative variables U to entropy vector, dS/dU, S = -rho*s/(kappa-1), s=ln(p)-kappa*ln(rho)
!==================================================================================================================================
PPURE SUBROUTINE ConsToEntropy_Core_Host(entropy,cons,Kappa,KappaM1,sKappaM1)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: cons    !< vector of conservative variables
REAL,INTENT(IN)                     :: Kappa,KappaM1,sKappaM1
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: entropy !< vector of entropy variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: vel(3),s,p,rho_p
!==================================================================================================================================
vel(:) = cons(MOMV)/cons(DENS)
p      = KappaM1*(cons(ENER)-0.5*SUM(cons(MOMV)*vel(:)))
! entropy: log(p) - eq.γ * log(ρ)
s      = LOG(p) - kappa*LOG(cons(DENS))
rho_p  = cons(DENS)/p

! Convert to entropy variables
entropy(DENS)      = (kappa-s)*skappaM1 - rho_p * 0.5*SUM(vel**2)  ! (γ - s) / (γ - 1) - (ρu^2 + ρv^2 + ρw^2) / ρ / 2p,
entropy(MOM1:MOM3) = rho_p*vel(1:3)        ! ρu / p
entropy(ENER)      = - rho_p          ! -ρ / p

END SUBROUTINE ConsToEntropy_Core_Host

!==================================================================================================================================
!> Transformation from primitive to conservative variables in the whole volume
!==================================================================================================================================
PPURE SUBROUTINE ConsToEntropy_Volume_Host(Nloc,nElems,entropy,cons,Kappa,KappaM1,sKappaM1)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                            !< local polynomial degree of solution representation
INTEGER,INTENT(IN) :: nElems
REAL,INTENT(OUT)   :: entropy(PP_nVar,0:Nloc,0:Nloc,0:Nloc,1:nElems)  !< vector of entropy variables
REAL,INTENT(IN)    :: cons(PP_nVar   ,0:Nloc,0:Nloc,0:Nloc,1:nElems)  !< vector of conservative variables
REAL,INTENT(IN)    :: Kappa,KappaM1,sKappaM1
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iElem
!==================================================================================================================================

  DO iElem=1,nElems
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      CALL ConsToEntropy_Core_Host(entropy(:,i,j,k,iElem),cons(:,i,j,k,iElem),Kappa,KappaM1,sKappaM1)
    END DO; END DO; END DO
  END DO

END SUBROUTINE ConsToEntropy_Volume_Host


!----------------------------------------------------------------------------------------------------------------------------------
! EntropyToCons Methods
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
!> Entry point for EntropyToCons_Sides
!==================================================================================================================================
SUBROUTINE EntropyToCons_Side(Nloc,sideID,elemID,entropy,cons,d_entropy,d_cons)
! MODULES
USE MOD_EOS_Vars,  ONLY: KappaM1,Kappa,sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)      :: Nloc                            !< local polynomial degree of solution representation
INTEGER,INTENT(IN)      :: sideID                          !< Index of the side to convert
INTEGER,INTENT(IN)      :: elemID                          !< Index of the element this side belongs to. Zero for "sides only" arrays
REAL,INTENT(IN)         :: entropy(PP_nVar,0:Nloc,0:Nloc)  !< vector of entropy variables
REAL,INTENT(OUT)        :: cons(PP_nVar   ,0:Nloc,0:Nloc)  !< vector of conservative variables
INTEGER(C_INT),VALUE,INTENT(IN) :: d_entropy                       !< Device variable key for entropy array
INTEGER(C_INT),VALUE,INTENT(IN) :: d_cons                          !< Device variable key for cons array
!----------------------------------------------------------------------------------------------------------------------------------

#if (USE_ACCEL == ACCEL_OFF)
  CALL EntropyToCons_Side_Host(Nloc,sideID,elemID,entropy,cons,Kappa,KappaM1,sKappaM1)
#else /* USE_ACCEL */
  CALL EntropyToCons_Side_Device(Nloc,sideID,elemID,d_entropy,d_cons,Kappa,KappaM1,sKappaM1)
#endif /* USE_ACCEL */

END SUBROUTINE EntropyToCons_Side

!==================================================================================================================================
!> Transformation from entropy to conservative variables U, dS/dU, S = -rho*s/(kappa-1), s=ln(p)-kappa*ln(rho)
!==================================================================================================================================
PPURE SUBROUTINE EntropyToCons_Core_Host(entropy,cons,Kappa,KappaM1,sKappaM1)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)   :: entropy !< vector of entropy variables
REAL,INTENT(IN)                      :: KappaM1,kappa,sKappaM1
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT)  :: cons    !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: s,entropy2(PP_nVar),rhoe
!==================================================================================================================================
entropy2 = entropy*kappaM1
s        = kappa - entropy2(DENS) + 0.5 * SUM(entropy2(MOMV)**2) / entropy2(ENER)
rhoe     = (kappaM1 / ((-entropy2(ENER))**kappa))**(skappaM1) * EXP(-s*skappaM1)

cons(DENS) = - rhoe * entropy2(ENER) ! ρ = -p * W[5]
cons(MOMV) = rhoe * entropy2(MOMV)
cons(ENER) = rhoe * (1 - SUM(entropy2(MOMV)**2) * 0.5/ entropy2(ENER)) !sKappaM1*p+0.5*SUM(cons(MOMV)*vel)

END SUBROUTINE EntropyToCons_Core_Host

!> Transformation from primitive to conservative variables on a single side
!==================================================================================================================================
PPURE SUBROUTINE EntropyToCons_Side_Host(Nloc,sideID,elemID,entropy,cons,Kappa,KappaM1,sKappaM1)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                            !< local polynomial degree of solution representation
INTEGER,INTENT(IN) :: sideID                          !< (only passed for consistency with device API)
INTEGER,INTENT(IN) :: elemID                          !< (only passed for consistency with device API)
REAL,INTENT(IN)    :: entropy(PP_nVar,0:Nloc,0:Nloc)  !< vector of entropy variables
REAL,INTENT(OUT)   :: cons(PP_nVar   ,0:Nloc,0:Nloc)  !< vector of conservative variables
REAL,INTENT(IN)    :: KappaM1,Kappa,sKappaM1
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,iSide
!==================================================================================================================================

  DO q=0,Nloc; DO p=0,Nloc
    CALL EntropyToCons_Core_Host(entropy(:,p,q),cons(:,p,q),Kappa,KappaM1,sKappaM1)
  END DO; END DO ! p,q=0,Nloc

END SUBROUTINE EntropyToCons_Side_Host


!==================================================================================================================================
!> Riemann solver function to get pressure at BCs
!==================================================================================================================================
PPURE FUNCTION PRESSURE_RIEMANN(U_Prim)
!==================================================================================================================================
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1,sKappaM1,sKappaP1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN) :: U_Prim(PP_nVarPrim) !< vector of primitive variables
REAL            :: PRESSURE_RIEMANN    !< pressure as the return value of the Riemann problem
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: kappaFac,ar,br,P_RP
!==================================================================================================================================
kappaFac=2.*Kappa*sKappaM1
IF(U_Prim(VEL1) .LE. 0.)THEN ! rarefaction
  P_RP=U_Prim(PRES) * MAX(0.0001,(1.+0.5*KappaM1*U_Prim(VEL1)/SQRT(Kappa*U_Prim(PRES)/U_Prim(DENS))))**kappaFac
ELSE ! shock
  ar=2.*sKappaP1/U_Prim(DENS)
  br=KappaM1*sKappaP1*U_Prim(PRES)
  P_RP=U_Prim(PRES)+U_Prim(VEL1)/ar*0.5*(U_Prim(VEL1)+SQRT(U_Prim(VEL1)*U_Prim(VEL1)+4.*ar*(U_Prim(PRES)+br)))
END IF
PRESSURE_RIEMANN=P_RP
END FUNCTION PRESSURE_RIEMANN


END MODULE MOD_EOS
