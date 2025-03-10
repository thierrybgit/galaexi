#define COPYRIGHT "=============================================================================================================================\
! Copyright (c) 2010-2022  Prof. Claus-Dieter Munz                                                                                                \
! Copyright (c) 2022-2024 Prof. Andrea Beck                                                                                                       \
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.                   \
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/                                                    \
!                                                                                                                                                 \
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License                              \
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.                           \
!                                                                                                                                                 \
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty                            \
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.                               \
!                                                                                                                                                 \
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.                 \
!================================================================================================================================               "

#ifdef SUN /*Abbrevations*/
#  define __DATE__ '__TIME__ and __DATE__ not'
#  define __TIME__ 'available for SUN COMPILER'
#  define IEEE_ISNAN
#elif SX
#  define __DATE__ '__TIME__ and __DATE__ not'
#  define __TIME__ 'available for SX COMPILER'
#elif PGI
#  define NO_ISNAN
#endif
#ifndef __FILENAME__
#define __FILENAME__ __FILE__
#endif
#define __STAMP__ __FILENAME__,__LINE__,__DATE__,__TIME__

#ifdef GNU
#  define IEEE_IS_NAN ISNAN
#endif

#ifdef NVHPC /*Calculate compiler version*/
#define NVHPC_VERSION (NVHPC_MAJOR * 10000 + NVHPC_MINOR * 100 + NVHPC_PATCH)
#if GCC_VERSION < 230000
#define NORM2(x) SQRT(SUM(x**2))
#endif
#endif

#define SIZEOF_F(x) (STORAGE_SIZE(x)/8)
#define NO_OP(x)    ASSOCIATE( x => x ); END ASSOCIATE

#ifdef GNU
#define CHECKSAFEINT(x,k)  IF(x>HUGE(INT( 1,KIND=k)).OR.x<-HUGE(INT( 1,KIND=k))) CALL Abort(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(REAL(1,KIND=k)).OR.x<-HUGE(REAL(1,KIND=k))) CALL Abort(__STAMP__,'Real conversion failed: out of range!')
#elif CRAY
#define CHECKSAFEINT(x,k)
#define CHECKSAFEREAL(x,k)
#else
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  ## k).OR.x<-HUGE(1_  ## k)) CALL Abort(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ ## k).OR.x<-HUGE(1._ ## k)) CALL Abort(__STAMP__,'Real conversion failed: out of range!')
#endif

#define DT_NVAR       3 /*Time Step Minimum: dt_Min*/
#define DT_MIN        1
#define DT_ANALYZE    2
#define DT_END        3


#define ALMOSTEQUALABSOLUTE(x,y,tol)  (ABS((x)-(y)).LE.(tol)) /*Test for equality: read description in mathtools.f90 for further info*/
#define ALMOSTEQUALRELATIVE(x,y,tol)  (ABS((x)-(y)).LE.MAX(ABS(x),ABS(y))*(tol))
#define ALMOSTEQUALABSORREL(x,y,tol)  (ALMOSTEQUALABSOLUTE(x,y,tol) .OR.  ALMOSTEQUALRELATIVE(x,y,tol))
#define ALMOSTEQUALABSANDREL(x,y,tol) (ALMOSTEQUALABSOLUTE(x,y,tol) .AND. ALMOSTEQUALRELATIVE(x,y,tol))

#if USE_MPI /*Define MPI specific write shortcuts*/
#  define SWRITE IF(MPIRoot) WRITE
#  define IPWRITE(a,b) WRITE(a,b)myRank,
#  define GETTIME(a) a=MPI_WTIME()
#else
#  define SWRITE WRITE
#  define IPWRITE WRITE
#  define GETTIME(a) CALL CPU_TIME(a)
#endif
#define ERRWRITE(a,b) CALL CreateErrFile(); IF(ErrorFiles) WRITE(UNIT_errOut,b)
#define LOGWRITE(a,b) IF(Logging) WRITE(UNIT_logOut,b)
#define SDEALLOCATE(A) IF(ALLOCATED(A))  DEALLOCATE(A)
#define ADEALLOCATE(A) IF(ASSOCIATED(A)) DEALLOCATE(A)

#if USE_OPENMP /*Define OpenMP specific shortcuts*/
#  define OMP_FLEXITIME() OMP_GET_WTIME()
#else
#  define OMP_FLEXITIME() FLEXITIME()
#endif

#define PP_IJK     i,j,k /*Loop variables*/
#define PP_ij      i,j

#define XI_MINUS   5 /*Predefined "PARAMETER-like" variables*/
#define XI_PLUS    3
#define ETA_MINUS  2
#define ETA_PLUS   4
#define ZETA_MINUS 1
#define ZETA_PLUS  6

#define S2E_ELEM_ID        1 /*Entry position in SideToElem*/
#define S2E_NB_ELEM_ID     2
#define S2E_LOC_SIDE_ID    3
#define S2E_NB_LOC_SIDE_ID 4
#define S2E_FLIP           5

#define E2S_SIDE_ID 1 /*Entry position in ElemToSide*/
#define E2S_FLIP    2

#define MI_SIDEID 1 /*Entry position in BC*/
#define MI_FLIP   2

#define BC_TYPE  1 /*Entry position in BC*/
#define BC_STATE 2
#define BC_ALPHA 3

#define SEND 1 /*Entry position in BC*/
#define RECV 2

#if FV_ENABLED
#define FV_SIZE 1
#else
#define FV_SIZE 0
#endif

#if !(FV_ENABLED)
#define FV_Elems(x) 0
#define FV_Elems_master(x) 0
#define FV_Elems_slave(x) 0
#define FV_Elems_Sum(x) 0
#endif

#if (FV_ENABLED>=2) || !PARABOLIC /*Compute viscous contributions in volume integral*/
#define VOLINT_VISC 0
#else
#define VOLINT_VISC 1
#endif

#define KILL(x) SWRITE(*,*) __FILE__,__LINE__,x; stop

#define OVERINTEGRATIONTYPE_NONE       0 /*overintegration*/
#define OVERINTEGRATIONTYPE_CUTOFF     1
#define OVERINTEGRATIONTYPE_CONSCUTOFF 2

#define FILTERTYPE_NONE   0 /*filter*/
#define FILTERTYPE_CUTOFF 1
#define FILTERTYPE_MODAL  2
#define FILTERTYPE_LAF    3

#define PRM_SPLITDG_SD 0 /*Split DG formulation options*/
#define PRM_SPLITDG_MO 1
#define PRM_SPLITDG_DU 2
#define PRM_SPLITDG_KG 3
#define PRM_SPLITDG_PI 4
#define PRM_SPLITDG_CH 5

#define PRM_RIEMANN_SAME         -1
#define PRM_RIEMANN_LF            1
#define PRM_RIEMANN_HLLC          2
#define PRM_RIEMANN_ROE           3
#define PRM_RIEMANN_ROEL2         32
#define PRM_RIEMANN_ROEENTROPYFIX 33
#define PRM_RIEMANN_HLL           4
#define PRM_RIEMANN_HLLE          5
#define PRM_RIEMANN_HLLEM         6
#ifdef SPLIT_DG
#define PRM_RIEMANN_CH            7
#define PRM_RIEMANN_Average       0
#endif

#if DEBUG /*PURE debug switch*/
#define PPURE
#else
#define PPURE PURE
#endif

#define ZDIM(a) a
#define PP_NZ   PP_N
#define DIMV    1:3

#define ACCEL_OFF 0 /*GPU options*/
#define ACCEL_CUDA 1
#define ACCEL_HIP 2

#undef COPYRIGHT
