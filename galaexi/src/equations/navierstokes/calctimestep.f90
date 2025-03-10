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
!> This module contains the routines to calculate the equation system specific allowable timestep.
!==================================================================================================================================
MODULE MOD_CalcTimeStep
! MODULES
USE ISO_C_BINDING, ONLY: C_NULL_CHAR, C_INT, C_Loc
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! Device side method interfaces
#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
  SUBROUTINE InitCalcTimestep_Device(Nloc,nElems,d_sJ,d_Mf,d_Mg,d_Mh,d_MetAdv &
#if PARABOLIC
                                    ,d_MetVisc,KappasPr &
#endif /*PARABOLIC*/
                                    ) BIND(C, NAME="InitCalcTimestep_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: Nloc
    INTEGER(C_INT),VALUE :: nElems
    INTEGER(C_INT),VALUE :: d_sJ
    INTEGER(C_INT),VALUE :: d_Mf
    INTEGER(C_INT),VALUE :: d_Mg
    INTEGER(C_INT),VALUE :: d_Mh
    INTEGER(C_INT),VALUE :: d_MetAdv
#if PARABOLIC
    INTEGER(C_INT),VALUE :: d_MetVisc
    REAL(C_DOUBLE),VALUE :: KappasPr
#endif
  END SUBROUTINE InitCalcTimestep_Device
END INTERFACE

INTERFACE
  SUBROUTINE CalcTimeStep_Device(Nloc,nElems,d_errType,d_TSconv,d_TSvisc,d_U,d_sJ,d_Mf,d_Mg,d_Mh, & 
                                d_CFLScale,d_MetAdv,Kappa,KappaM1,R &
#if PARABOLIC
                                ,d_MetVisc,d_DFLScale &
#endif /*PARABOLIC*/
#if FV_ENABLED
                                ,d_FV_Elems &
#if FV_ENABLED == 2
                                ,d_FV_alpha,FV_alpha_min &
#endif /* FV_ENABLED == 2 */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
                                ,d_muSGS &
#endif /* EDDYVISCOSITY */
                                ) BIND(C, NAME="CalcTimeStep_Device")
  USE ISO_C_BINDING
  IMPLICIT NONE
    INTEGER(C_INT),VALUE :: Nloc
    INTEGER(C_INT),VALUE :: nElems
    INTEGER(C_INT),VALUE :: d_errType
    INTEGER(C_INT),VALUE :: d_TSconv
    INTEGER(C_INT),VALUE :: d_TSvisc
    INTEGER(C_INT),VALUE :: d_U
    INTEGER(C_INT),VALUE :: d_sJ
    INTEGER(C_INT),VALUE :: d_Mf
    INTEGER(C_INT),VALUE :: d_Mg
    INTEGER(C_INT),VALUE :: d_Mh
    INTEGER(C_INT),VALUE :: d_CFLScale
    INTEGER(C_INT),VALUE :: d_MetAdv
    REAL(C_DOUBLE),VALUE :: Kappa
    REAL(C_DOUBLE),VALUE :: KappaM1
    REAL(C_DOUBLE),VALUE :: R
#if PARABOLIC
    INTEGER(C_INT),VALUE :: d_MetVisc
    INTEGER(C_INT),VALUE :: d_DFLScale
#endif /*PARABOLIC*/
#if FV_ENABLED
    INTEGER(C_INT),VALUE :: d_FV_Elems
#if FV_ENABLED == 2
    INTEGER(C_INT),VALUE :: d_FV_alpha
    REAL(C_DOUBLE),VALUE :: FV_alpha_min
#endif /* FV_ENABLED == 2 */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
    INTEGER(C_INT),VALUE :: d_muSGS
#endif /* EDDYVISCOSITY */
    
  END SUBROUTINE CalcTimeStep_Device
END INTERFACE

#endif

PUBLIC :: InitCalctimestep,CalcTimeStep,FinalizeCalctimestep
!==================================================================================================================================

REAL,ALLOCATABLE :: MetricsAdv(:,:,:,:,:,:)  !< support variable: NORM2(Metricsfgh)/J
INTEGER(C_INT)   :: d_MetricsAdv 
REAL             :: TimeStepConv
REAL             :: TimeStepVisc
INTEGER          :: TimeStepErrType
INTEGER(C_INT)   :: d_TimeStepErrType 
INTEGER(C_INT)   :: d_TimeStepConv 
INTEGER(C_INT)   :: d_TimeStepVisc 
#if PARABOLIC
REAL,ALLOCATABLE :: MetricsVisc(:,:,:,:,:,:) !< support variable: kappa/Pr*(SUM((Metricsfgh/J)**2))
INTEGER(C_INT)   :: d_MetricsVisc 
#endif

CONTAINS

!==================================================================================================================================
!> Precompute some metric support variables
!==================================================================================================================================
SUBROUTINE InitCalctimestep()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
#if PARABOLIC
USE MOD_EOS_Vars ,ONLY:KappasPr
#endif /*PARABOLIC*/
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_DeviceMem
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem,FVE
#if PARABOLIC
REAL                         :: KappasPr_max
#endif /*PARABOLIC*/
!==================================================================================================================================

  ALLOCATE(MetricsAdv(3,0:PP_N,0:PP_N,0:PP_N,nElems,0:FV_SIZE))
#if PARABOLIC
  ALLOCATE(MetricsVisc(3,0:PP_N,0:PP_N,0:PP_N,nElems,0:FV_SIZE))
#endif /*PARABOLIC*/

#if PARABOLIC
  KappasPr_max=KAPPASPR_MAX_TIMESTEP_H()
#endif /*PARABOLIC*/

  DO FVE=0,FV_SIZE
    DO iElem=1,nElems
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        MetricsAdv(1,i,j,k,iElem,FVE)=sJ(i,j,k,iElem,FVE)*NORM2(Metrics_fTilde(:,i,j,k,iElem,FVE))
        MetricsAdv(2,i,j,k,iElem,FVE)=sJ(i,j,k,iElem,FVE)*NORM2(Metrics_gTilde(:,i,j,k,iElem,FVE))
        MetricsAdv(3,i,j,k,iElem,FVE)=sJ(i,j,k,iElem,FVE)*NORM2(Metrics_hTilde(:,i,j,k,iElem,FVE))
#if PARABOLIC
        MetricsVisc(1,i,j,k,iElem,FVE)=KappasPR_max*(SUM((Metrics_fTilde(:,i,j,k,iElem,FVE)*sJ(i,j,k,iElem,FVE))**2))
        MetricsVisc(2,i,j,k,iElem,FVE)=KappasPR_max*(SUM((Metrics_gTilde(:,i,j,k,iElem,FVE)*sJ(i,j,k,iElem,FVE))**2))
        MetricsVisc(3,i,j,k,iElem,FVE)=KappasPR_max*(SUM((Metrics_hTilde(:,i,j,k,iElem,FVE)*sJ(i,j,k,iElem,FVE))**2))
#endif /*PARABOLIC*/
      END DO; END DO; END DO
    END DO
  END DO

#if (USE_ACCEL != ACCEL_OFF)
  CALL AllocateDeviceMemory(d_MetricsAdv, SIZE_C_DOUBLE, SIZE(MetricsAdv))
  CALL CopyToDevice(d_MetricsAdv, C_Loc(MetricsAdv), SIZE_C_DOUBLE, SIZE(MetricsAdv))

  CALL AllocateDeviceMemory(d_TimeStepConv, SIZE_C_DOUBLE, 1)
  CALL AllocateDeviceMemory(d_TimeStepVisc, SIZE_C_DOUBLE, 1)
  CALL AllocateDeviceMemory(d_TimeStepErrType, SIZE_C_INT, 1)
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
  ! On hybrid architectures, we have to call the copy methods to associate
  ! the device pointer with the host memory
  CALL CopyToDevice(d_TimeStepConv, C_Loc(TimeStepConv), SIZE_C_DOUBLE, 1)
  CALL CopyToDevice(d_TimeStepVisc, C_Loc(TimeStepVisc), SIZE_C_DOUBLE, 1)
  CALL CopyToDevice(d_TimeStepErrType, C_Loc(TimeStepErrType), SIZE_C_INT, 1)
#endif /* USE_HYBRID */
#if PARABOLIC
  CALL AllocateDeviceMemory(d_MetricsVisc, SIZE_C_DOUBLE, SIZE(MetricsVisc))
  CALL CopyToDevice(d_MetricsVisc, C_Loc(MetricsVisc), SIZE_C_DOUBLE, SIZE(MetricsVisc))
#endif /*PARABOLIC*/
#endif /* USE_ACCEL */

END SUBROUTINE InitCalctimestep

!==================================================================================================================================
!> Entry point for the CalcTimeStep backends
!==================================================================================================================================
SUBROUTINE CalcTimeStep(minTimeStep,errType,Nloc)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY:U
USE MOD_EOS_Vars     ,ONLY:Kappa,KappaM1,R
USE MOD_Mesh_Vars    ,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
USE MOD_TimeDisc_Vars,ONLY:CFLScale,ViscousTimeStep,dtElem
#ifndef GNU
USE, INTRINSIC :: IEEE_ARITHMETIC,ONLY:IEEE_IS_NAN
#endif
#if PARABOLIC
USE MOD_TimeDisc_Vars,ONLY:DFLScale
#endif /*PARABOLIC*/
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems
#if FV_ENABLED == 2
USE MOD_FV_Vars      ,ONLY: FV_alpha,FV_alpha_min
#endif
#endif
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars, ONLY: muSGS
#endif
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_DeviceMem
USE MOD_DG_Vars      ,ONLY: d_U
USE MOD_EOS_Vars     ,ONLY: Kappa,KappaM1,R
USE MOD_Mesh_Vars    ,ONLY: d_sJ,d_Metrics_fTilde,d_Metrics_gTilde,d_Metrics_hTilde
USE MOD_TimeDisc_Vars,ONLY: d_CFLScale
#if PARABOLIC
USE MOD_TimeDisc_Vars,ONLY: d_DFLScale
#endif /*PARABOLIC*/
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: d_FV_Elems
#if FV_ENABLED == 2
USE MOD_FV_Vars      ,ONLY: d_FV_alpha
#endif
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars, ONLY: d_muSGS
#endif
#endif /* USE_ACCEL */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)             :: minTimeStep
INTEGER, INTENT(OUT)         :: errType
INTEGER,INTENT(IN),OPTIONAL  :: Nloc      !< Local polynomail order. Only passed for unit test. Is never passed in main code.
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: locN
REAL    :: TimeStep(3)
!==================================================================================================================================

  locN = PP_N
  IF(PRESENT(Nloc)) locN = Nloc

  TimeStepErrType = 0

#if (USE_ACCEL == ACCEL_OFF)

  CALL CalcTimeStep_Host(locN,nElems,TimeStepErrType,TimeStepConv,TimeStepVisc,U,sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde, &
                          CFLScale,dtElem,MetricsAdv,Kappa,KappaM1,R &
#if PARABOLIC
                          ,MetricsVisc,DFLScale &
#endif /*PARABOLIC*/
#if FV_ENABLED
                          ,FV_Elems &
#if FV_ENABLED == 2
                          ,FV_alpha,FV_alpha_min &
#endif /* FV_ENABLED == 2 */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
                          ,muSGS &
#endif /* EDDYVISCOSITY */
                        )

#else /* USE_ACCEL */

  CALL CalcTimeStep_Device(locN,nElems,d_TimeStepErrType,d_TimeStepConv,d_TimeStepVisc,d_U,d_sJ, &
                          d_Metrics_fTilde,d_Metrics_gTilde,d_Metrics_hTilde,d_CFLScale,d_MetricsAdv,Kappa,KappaM1,R &
#if PARABOLIC
                          ,d_MetricsVisc,d_DFLScale &
#endif /*PARABOLIC*/
#if FV_ENABLED
                          ,d_FV_Elems &
#if FV_ENABLED == 2
                          ,d_FV_alpha,FV_alpha_min &
#endif /* FV_ENABLED == 2 */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
                          ,d_muSGS &
#endif /* EDDYVISCOSITY */
                        )
  CALL CopyFromDevice(C_Loc(TimeStepConv), d_TimeStepConv, SIZE_C_DOUBLE, 1)
  CALL CopyFromDevice(C_Loc(TimeStepVisc), d_TimeStepVisc, SIZE_C_DOUBLE, 1)
  CALL CopyFromDevice(C_Loc(TimeStepErrType), d_TimeStepErrType, SIZE_C_INT, 1)
#endif /* USE_ACCEL */

  ! Check for NaN timesteps
  IF(TimeStepErrType /= 0) THEN
    TimeStepErrType = 1
  ELSE
    IF(IEEE_IS_NAN(TimeStepConv))THEN
      TimeStepErrType = 2
    END IF
    IF(IEEE_IS_NAN(TimeStepVisc))THEN
      TimeStepErrType = 3
    END IF
  END IF

  TimeStep(1) = TimeStepConv
  TimeStep(2) = TimeStepVisc
#if USE_MPI
  TimeStep(3) = -TimeStepErrType ! reduce with timestep, minus due to MPI_MIN
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,TimeStep,3,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_FLEXI,iError)
  TimeStepErrType = INT(-TimeStep(3))
#endif /*USE_MPI*/
  ViscousTimeStep = (TimeStep(2) .LT. TimeStep(1))
  errType = TimeStepErrType
  minTimeStep = MINVAL(TimeStep(1:2))

END SUBROUTINE CalcTimeStep

#if (USE_ACCEL == ACCEL_OFF)
!==================================================================================================================================
!> Compute the time step for the current update of U for the Navier-Stokes-Equations
!==================================================================================================================================
SUBROUTINE CalcTimeStep_Host(Nloc,nElems,errType,TSconv,TSvisc,U,sJ,Mf,Mg,Mh,CFLScale,dtElem,MetAdv, &
                              Kappa,KappaM1,R &
#if PARABOLIC
                              ,MetVisc,DFLScale &
#endif
#if FV_ENABLED
                              ,FV_Elems &
#if FV_ENABLED == 2
                              ,FV_alpha,FV_alpha_min &
#endif
#endif
#if EDDYVISCOSITY
                              ,muSGS &
#endif
                            )
! MODULES
#ifndef GNU
USE, INTRINSIC :: IEEE_ARITHMETIC,ONLY:IEEE_IS_NAN
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)           :: Nloc
INTEGER,INTENT(IN)           :: nElems
INTEGER,INTENT(INOUT)        :: errType
REAL,INTENT(OUT)             :: TSconv
REAL,INTENT(OUT)             :: TSvisc
REAL,INTENT(IN)              :: U(PP_nVar,0:Nloc,0:Nloc,0:Nloc,nElems)
REAL,INTENT(IN)              :: sJ(0:Nloc,0:Nloc,0:Nloc,nElems,0:FV_SIZE)
REAL,DIMENSION(3,0:Nloc,0:Nloc,0:Nloc,nElems,0:FV_SIZE),INTENT(IN)  :: Mf, Mg, Mh
REAL,INTENT(IN)              :: CFLScale(0:FV_SIZE)
REAL,INTENT(INOUT)           :: dtElem(nElems)
REAL,INTENT(IN)              :: MetAdv(3,0:Nloc,0:Nloc,0:Nloc,nElems,0:FV_SIZE)
REAL,INTENT(IN)              :: Kappa,KappaM1,R
#if PARABOLIC
REAL,INTENT(IN)              :: MetVisc(3,0:Nloc,0:Nloc,0:Nloc,nElems,0:FV_SIZE)
REAL,INTENT(IN)              :: DFLScale(0:FV_SIZE)
#endif /* PARABOLIC */
#if FV_ENABLED
INTEGER,INTENT(IN)           :: FV_Elems(nElems)
#if FV_ENABLED == 2
REAL,INTENT(IN)              :: FV_alpha(1,0:Nloc,0:Nloc,0:Nloc,nElems)
REAL,INTENT(IN)              :: FV_alpha_min
#endif /* FV_ENABLED == 2 */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
REAL,INTENT(IN)              :: muSGS(1,0:Nloc,0:Nloc,0:Nloc,nElems)
#endif /* EDDYVISCOSITY */
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem
REAL,DIMENSION(PP_2Var)      :: UE
REAL                         :: Max_Lambda(3),c,vsJ(3)
#if EDDYVISCOSITY
REAL                         :: muSGSmax
#endif
#if PARABOLIC
REAL                         :: Max_Lambda_v(3),mu,prim(PP_nVarPrim)
#endif /*PARABOLIC*/
INTEGER                      :: FVE
!==================================================================================================================================

  errType = 0

  TSconv = HUGE(1.0)
  TSvisc = HUGE(1.0)

  DO iElem=1,nElems
#if FV_ENABLED
    FVE = FV_Elems(iElem)
#else
    FVE = 0
#endif
    Max_Lambda=0.
#if PARABOLIC
    Max_Lambda_v=0.
#endif /*PARABOLIC*/
#if EDDYVISCOSITY
    muSGSMax = MAXVAL(muSGS(1,:,:,:,iElem))
#endif

    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      ! TODO: ATTENTION: Temperature of UE not filled!!!
      UE(EXT_CONS)=U(:,i,j,k,iElem)
      UE(EXT_SRHO)=1./UE(EXT_DENS)
      UE(EXT_VELV)=VELOCITY_HE(UE)
      UE(EXT_PRES)=PRESSURE_HE(UE)
      UE(EXT_TEMP)=TEMPERATURE_HE(UE)

      if (IEEE_IS_NAN(U(DENS,i,j,k,iElem))) errType = 1

      c=SPEEDOFSOUND_HE(UE)
      vsJ=UE(EXT_VELV)*sJ(i,j,k,iElem,FVE)
      Max_Lambda(1)=MAX(Max_Lambda(1),ABS(SUM(Mf(:,i,j,k,iElem,FVE)*vsJ)) + &
                                                c*MetAdv(1,i,j,k,iElem,FVE))
      Max_Lambda(2)=MAX(Max_Lambda(2),ABS(SUM(Mg(:,i,j,k,iElem,FVE)*vsJ)) + &
                                                c*MetAdv(2,i,j,k,iElem,FVE))
      Max_Lambda(3)=MAX(Max_Lambda(3),ABS(SUM(Mh(:,i,j,k,iElem,FVE)*vsJ)) + &
                                                c*MetAdv(3,i,j,k,iElem,FVE))
#if PARABOLIC
      ! Viscous Eigenvalues
      prim = UE(EXT_PRIM)
      mu=VISCOSITY_PRIM(prim)
#if EDDYVISCOSITY
      mu = mu+muSGSMax
#endif
      Max_Lambda_v=MAX(Max_Lambda_v,mu*UE(EXT_SRHO)*MetVisc(:,i,j,k,iElem,FVE))
#endif /* PARABOLIC*/
    END DO; END DO; END DO ! i,j,k

#if FV_ENABLED == 2
    dtElem(iElem)=MERGE(CFLScale(0),CFLScale(1),FV_alpha(iElem).LE.FV_alpha_min)*2./SUM(Max_Lambda)
#elif FV_ENABLED == 3
    dtElem(iElem)=MERGE(CFLScale(0),CFLScale(1),ALL(FV_alpha(:,:,:,:,iElem).LE.FV_alpha_min))*2./SUM(Max_Lambda)
#else
    dtElem(iElem)=CFLScale(FVE)*2./SUM(Max_Lambda)
#endif
    TSconv=MIN(TSconv,dtElem(iElem))

#if PARABOLIC
    IF(SUM(Max_Lambda_v).GT.0.)THEN
#if FV_ENABLED == 2 || FV_ENABLED == 3
      dtElem(iElem)=MIN(dtElem(iElem),MINVAL(DFLScale(:))*4./SUM(Max_Lambda_v))
      TSvisc= MIN(TSvisc, MINVAL(DFLScale(:))*4./SUM(Max_Lambda_v))
#else
      dtElem(iElem)=MIN(dtElem(iElem),DFLScale(FVE)*4./SUM(Max_Lambda_v))
      TSvisc= MIN(TSvisc, DFLScale(FVE)*4./SUM(Max_Lambda_v))
#endif
    END IF
#endif /* PARABOLIC*/

  END DO ! iElem=1,nElems

END SUBROUTINE CalcTimeStep_Host
#endif

!==================================================================================================================================
!> Deallocate CalcTimeStep arrays
!==================================================================================================================================
SUBROUTINE FinalizeCalctimestep()
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(MetricsAdv)
#if PARABOLIC
SDEALLOCATE(MetricsVisc)
#endif
END SUBROUTINE

END MODULE MOD_CalcTimeStep
