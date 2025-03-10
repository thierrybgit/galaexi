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
!> \brief Contains the definitions of the physical fluxes of the equation system.
!>
!> The routine EvalFlux3D will compute the advection (Euler) part only, and can be called either for a single point or for
!> a volume cell. The fluxes are computed in three spatial dimension - for 2D computations, the fluxes in the third dimension
!> will always be set to 0.
!> EvalDiffFlux3D will do the same thing, but compute only the diffusive part of the fluxes. Additionally, a routine to compute
!> the fluxes on a single side is provided (used in the riemann routines).
!> The EvalEulerFlux1D routines are used in the Riemann solver, where only a flux in one spatial dimension is needed.
!>
!> The flux definitions are only done once in the single point routines, all other (side, volume) routines will simply wrap
!> to this definition.
!==================================================================================================================================
MODULE MOD_Flux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE EvalEulerFlux1D_fast
  MODULE PROCEDURE EvalEulerFlux1D_Host
END INTERFACE

#if PARABOLIC
INTERFACE EvalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux3D_Point
  MODULE PROCEDURE EvalDiffFlux3D_Surface
  MODULE PROCEDURE EvalDiffFlux3D_Volume
#if FV_ENABLED
  MODULE PROCEDURE EvalDiffFlux3D_Volume_FV
#endif /*FV_ENABLED*/
END INTERFACE
#endif /*PARABOLIC*/

!----------------------------------------------------------------------------------------------------------------------------------
! DEVICE METHOD INTERFACES
!----------------------------------------------------------------------------------------------------------------------------------
#if (USE_ACCEL != ACCEL_OFF)
INTERFACE
    SUBROUTINE EvalTransformedFlux3D_Device(Nloc, nElems, d_U, d_UPrim &
#if PARABOLIC
                                          ,d_gradUx,d_gradUy,d_gradUz &
#endif
                                          ,d_F,d_G,d_H,d_Metrics_fTilde,d_Metrics_gTilde,d_Metrics_hTilde, myStream) &
              BIND(C, NAME="EvalTransformedFlux3D_Device")
    USE ISO_C_BINDING
    IMPLICIT NONE
      INTEGER(C_INT), VALUE :: Nloc
      INTEGER(C_INT), VALUE :: nElems
      INTEGER(C_INT),VALUE  :: d_U
      INTEGER(C_INT),VALUE  :: d_Uprim
#if PARABOLIC
      INTEGER(C_INT),VALUE  :: d_gradUx,d_gradUy,d_gradUz
#endif
      INTEGER(C_INT),VALUE  :: d_F, d_G, d_H
      INTEGER(C_INT),VALUE  :: d_Metrics_fTilde,d_Metrics_gTilde,d_Metrics_hTilde
      INTEGER(C_INT), VALUE :: myStream
    END SUBROUTINE EvalTransformedFlux3D_Device
END INTERFACE
#endif

PUBLIC:: EvalEulerFlux1D_fast, EvalTransformedFlux3D
#if PARABOLIC
PUBLIC::EvalDiffFlux3D
#endif /*PARABOLIC*/
!==================================================================================================================================

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
! EvalTransformedFlux3D METHODS
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
!> Entry point for EvalTransformedFlux3D. Chooses which backend to call for the computation.
!==================================================================================================================================
SUBROUTINE EvalTransformedFlux3D(Nloc,nElems, U,UPrim &
#if PARABOLIC
                                  ,gradUx,gradUy,gradUz &
#endif
                                  ,f,g,h,Mf,Mg,Mh,streamID)
! MODULES
USE MOD_DEVICE, ONLY: STREAM_DEFAULT
#if (USE_ACCEL != ACCEL_OFF)
USE MOD_DG_Vars, ONLY: d_U, d_Uprim, d_F, d_G, d_H
USE MOD_Mesh_Vars, ONLY: d_Metrics_fTilde, d_Metrics_gTilde, d_Metrics_hTilde
#if PARABOLIC
USE MOD_Lifting_Vars, ONLY: d_gradUx, d_gradUy, d_gradUz
#endif
#endif /* USE_ACCEL */
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER                              ,INTENT(IN)  :: Nloc                   !< polynomial order
INTEGER                              ,INTENT(IN)  :: nElems          !< Number of elements in the current block being passed to the kernel
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:Nloc,nElems),INTENT(IN)  :: U                      !< Conservative solution
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:Nloc,0:Nloc,nElems),INTENT(IN)  :: UPrim                  !< Primitive solution
REAL,DIMENSION(3          ,0:Nloc,0:Nloc,0:Nloc,nElems),INTENT(IN)  :: Mf,Mg,Mh               !< Metrics in x,y,z
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:Nloc,nElems),INTENT(OUT) :: f,g,h                  !> Physical fluxes in x,y,z
#if PARABOLIC
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:Nloc,0:Nloc,nElems),INTENT(IN)  :: gradUx,gradUy,gradUz   !> gradients in x,y,z
#endif
INTEGER,OPTIONAL,INTENT(IN) :: streamID                                     !> Index for the device stream to perform the calc on
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: mystream
!==================================================================================================================================
  mystream=STREAM_DEFAULT
  IF (PRESENT(streamID)) mystream=streamID

#if (USE_ACCEL == ACCEL_OFF)

  ! Call the CPU backend
  CALL EvalTransformedFlux3D_Host((Nloc+1)**3*nElems,U,UPrim &
#if PARABOLIC
                                  ,gradUx,gradUy,gradUz &
#endif
                                  ,F,G,H,Mf,Mg,Mh)

#else /* (USE_ACCEL != ACCEL_OFF) */
  ! Call the device backend
  CALL EvalTransformedFlux3D_Device(Nloc,nElems, d_U, d_UPrim &
#if PARABOLIC
                                    ,d_gradUx,d_gradUy,d_gradUz &
#endif
                                    ,d_F,d_G,d_H,d_Metrics_fTilde,d_Metrics_gTilde,d_Metrics_hTilde, mystream)
#endif /* USE_ACCEL */

END SUBROUTINE EvalTransformedFlux3D

#if (USE_ACCEL == ACCEL_OFF)
!==================================================================================================================================
!> Computes the advection part of the Navier-Stokes fluxes on the CPU
!==================================================================================================================================
SUBROUTINE EvalTransformedFlux3D_Host(nDOF,U,UPrim &
#if PARABOLIC
                                  ,gradUx,gradUy,gradUz &
#endif
                                  ,F,G,H,Mf,Mg,Mh)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER                              ,INTENT(IN)  :: nDOF                   !< number of degrees of freedom in arrays
REAL,DIMENSION(PP_nVar       ,1:nDOF),INTENT(IN)  :: U                      !< Conservative solution
REAL,DIMENSION(PP_nVarPrim   ,1:nDOF),INTENT(IN)  :: UPrim                  !< Primitive solution
REAL,DIMENSION(3             ,1:nDOF),INTENT(IN)  :: Mf,Mg,Mh               !< Metrics in x,y,z
REAL,DIMENSION(PP_nVar       ,1:nDOF),INTENT(OUT) :: F,G,H                  !> Physical fluxes in x,y,z
#if PARABOLIC
REAL,DIMENSION(PP_nVarLifting,1:nDOF),INTENT(IN)  :: gradUx,gradUy,gradUz   !> gradients in x,y,z
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
!==================================================================================================================================

  DO i = 1,nDOF
#if PARABOLIC
    ! COME BACK LATER AND ADD THIS DURING WORK ON EvalTransformedDiffFlux
#else
    CALL EvalTransformedEulerFlux3D_fast_Host(U(:,i),UPrim(:,i),F(:,i),G(:,i),H(:,i),Mf(:,i),Mg(:,i),Mh(:,i))
#endif
  END DO

END SUBROUTINE EvalTransformedFlux3D_Host

!==================================================================================================================================
!> CPU backend core of the computation advection part of the Navier-Stokes fluxes in all space dimensions using the conservative  
!> and primitive variables.
!==================================================================================================================================
SUBROUTINE EvalTransformedEulerFlux3D_fast_Host(U,UPrim,F,G,H,Mf,Mg,Mh)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U        !< Conservative solution
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim    !< Primitive solution
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: F,G,H    !> Physical fluxes in x/y/z direction
REAL,DIMENSION(3          ),INTENT(IN)  :: Mf,Mg,Mh !> Metrics in x/y/z direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Ep
REAL                :: Mmom
!==================================================================================================================================

  ! auxiliary variables
  Ep = (U(ENER) + UPrim(PRES))/U(DENS)
  
  ! fluxes in x-direction
  Mmom    = DOT_PRODUCT(Mf(:),U(MOMV)) ! metrics times momentum vector
  F(DENS) =  Mmom
  F(MOM1) =  Mmom*UPrim(VEL1) &
          + Mf(1)*UPrim(PRES)
  F(MOM2) =  Mmom*UPrim(VEL2) &
          + Mf(2)*UPrim(PRES)
  F(MOM3) =  Mmom*UPrim(VEL3) &
          + Mf(3)*UPrim(PRES)
  F(ENER) = Mmom*Ep
  
  ! fluxes in y-direction
  Mmom    = DOT_PRODUCT(Mg(:),U(MOMV)) ! metrics times momentum vector
  G(DENS) =  Mmom
  G(MOM1) =  Mmom*UPrim(VEL1) &
          + Mg(1)*UPrim(PRES)
  G(MOM2) =  Mmom*UPrim(VEL2) &
          + Mg(2)*UPrim(PRES)
  G(MOM3) =  Mmom*UPrim(VEL3) &
          + Mg(3)*UPrim(PRES)
  G(ENER) =  Mmom*Ep
  
  ! fluxes in z-direction
  Mmom    = DOT_PRODUCT(Mh(:),U(MOMV)) ! metrics times momentum vector
  H(DENS) =  Mmom
  H(MOM1) =  Mmom*UPrim(VEL1) &
          + Mh(1)*UPrim(PRES)
  H(MOM2) =  Mmom*UPrim(VEL2) &
          + Mh(2)*UPrim(PRES)
  H(MOM3) =  Mmom*UPrim(VEL3) &
          + Mh(3)*UPrim(PRES)
  H(ENER) =  Mmom*Ep

END SUBROUTINE EvalTransformedEulerFlux3D_fast_Host
#endif /* USE_ACCEL == ACCEL_OFF */

!----------------------------------------------------------------------------------------------------------------------------------
! EvalDiffFlux3D METHODS
!----------------------------------------------------------------------------------------------------------------------------------
#if PARABOLIC
!==================================================================================================================================
!> Compute Navier-Stokes diffusive flux using the primitive variables and derivatives.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D_Point(UPrim,gradUx,gradUy,gradUz,f,g,h &
#if EDDYVISCOSITY
                                      ,muSGS &
#endif
)
! MODULES
USE MOD_Equation_Vars,ONLY: s23,s43
USE MOD_EOS_Vars,     ONLY: cp,Pr
USE MOD_Viscosity
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: PrSGS
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)   ,INTENT(IN)  :: UPrim                 !< Solution vector
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx,gradUy,gradUz  !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar)       ,INTENT(OUT) :: f,g,h                 !> Physical fluxes in x,y,z directions
#if EDDYVISCOSITY
REAL                          ,INTENT(IN)  :: muSGS                 !< SGS viscosity
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: muS,lambda
REAL                :: tau_xx,tau_yy,tau_xy
#if PP_dim==3
REAL                :: tau_zz,tau_xz,tau_yz
#endif
!==================================================================================================================================
! ideal gas law
muS    = VISCOSITY_PRIM(UPrim)
lambda = THERMAL_CONDUCTIVITY_H(muS)
!Add turbulent sub grid scale viscosity to mu
#if EDDYVISCOSITY
muS    = muS    + muSGS
lambda = lambda + muSGS*cp/PrSGS
#endif

ASSOCIATE( v1     => UPrim(VEL1),       v2     => UPrim(VEL2),       v3     => UPrim(VEL3), &
           gradT1 => GradUx(LIFT_TEMP), gradT2 => GradUy(LIFT_TEMP), gradT3 => GradUz(LIFT_TEMP) )
#if PP_dim==3
! gradients of primitive variables are directly available gradU = (/ drho, dv1, dv2, dv3, dT /)

! viscous fluxes in x-direction
tau_xx = muS * ( s43 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = muS * (-s23 * gradUx(LIFT_VEL1) + s43 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_zz = muS * (-s23 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) + s43 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x-2/3*mu*v_y +4/3*mu*w*z
tau_xy = muS * (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))               !mu*(u_y+v_x)
tau_xz = muS * (gradUz(LIFT_VEL1) + gradUx(LIFT_VEL3))               !mu*(u_z+w_x)
tau_yz = muS * (gradUz(LIFT_VEL2) + gradUy(LIFT_VEL3))               !mu*(y_z+w_y)

f(DENS) = 0.
f(MOM1) = -tau_xx                                       ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
f(MOM2) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
f(MOM3) = -tau_xz                                       ! F_euler-mu*(u_z+w_x)
f(ENER) = -tau_xx*v1-tau_xy*v2-tau_xz*v3-lambda*gradT1  ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) q_x=-lambda*T_x
! viscous fluxes in y-direction
g(DENS) = 0.
g(MOM1) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
g(MOM2) = -tau_yy                                       ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
g(MOM3) = -tau_yz                                       ! F_euler-mu*(y_z+w_y)
g(ENER) = -tau_xy*v1-tau_yy*v2-tau_yz*v3-lambda*gradT2  ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) q_y=-lambda*T_y
! viscous fluxes in z-direction
h(DENS) = 0.
h(MOM1) = -tau_xz                                       ! F_euler-mu*(u_z+w_x)
h(MOM2) = -tau_yz                                       ! F_euler-mu*(y_z+w_y)
h(MOM3) = -tau_zz                                       ! F_euler-4/3*mu*w_z+2/3*mu*(u_x+v_y)
h(ENER) = -tau_xz*v1-tau_yz*v2-tau_zz*v3-lambda*gradT3  ! F_euler-(tau_zx*u+tau_zy*v+tau_zz*w-q_z) q_z=-lambda*T_z
#else
! gradients of primitive variables are directly available gradU = (/ drho, dv1, dv2, dv3, dT /)

! viscous fluxes in x-direction
tau_xx = muS * ( s43 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2))  ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = muS * (-s23 * gradUx(LIFT_VEL1) + s43 * gradUy(LIFT_VEL2))  !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_xy = muS * (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))               !mu*(u_y+v_x)

f(DENS) = 0.
f(MOM1) = -tau_xx                                       ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
f(MOM2) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
f(MOM3) = 0.
f(ENER) = -tau_xx*v1-tau_xy*v2-lambda*gradT1            ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) q_x=-lambda*T_x
! viscous fluxes in y-direction
g(DENS) = 0.
g(MOM1) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
g(MOM2) = -tau_yy                                       ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
g(MOM3) = 0.
g(ENER) = -tau_xy*v1-tau_yy*v2-lambda*gradT2            ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) q_y=-lambda*T_y
! viscous fluxes in z-direction
h    = 0.
#endif
END ASSOCIATE
END SUBROUTINE EvalDiffFlux3D_Point


!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D_Volume(UPrim,gradUx,gradUy,gradUz,f,g,h,iElem)
! MODULES
USE MOD_PreProc
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: muSGS
#endif /*EDDYVISCOSITY*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim,   0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)  :: UPrim                !< Solution vector
REAL,DIMENSION(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)  :: gradUx,gradUy,gradUz !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar,       0:PP_N,0:PP_N,0:PP_NZ),INTENT(OUT) :: f,g,h                !> Physical fluxes in x,y,z directions
INTEGER                                             ,INTENT(IN)  :: iElem                !< element index in global array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k
!==================================================================================================================================
DO k=0,PP_NZ;  DO j=0,PP_N; DO i=0,PP_N
  CALL EvalDiffFlux3D_Point(Uprim(:,i,j,k),gradUx(:,i,j,k),gradUy(:,i,j,k),gradUz(:,i,j,k), &
                                                f(:,i,j,k),     g(:,i,j,k),     h(:,i,j,k)  &
#if EDDYVISCOSITY
                            ,muSGS(1,i,j,k,iElem)&
#endif /*EDDYVISCOSITY*/
                            )
END DO; END DO; END DO ! i,j,k

#if !EDDYVISCOSITY
NO_OP(iElem)
#endif /*!EDDYVISCOSITY*/
END SUBROUTINE EvalDiffFlux3D_Volume

#if FV_ENABLED
!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D_Volume_FV(UPrim,gradUx,gradUy,gradUz,f,g,h,iElem,PP_N_xi,PP_N_eta,PP_N_zeta)
! MODULES
USE MOD_PreProc
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: muSGS
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PRIM,0:PP_N_xi,0:PP_N_eta,0:PP_N_zeta),INTENT(IN)           :: UPrim                !< Solution vector
!> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVarLifting,0:PP_N_xi,0:PP_N_eta,0:PP_N_zeta),INTENT(IN) :: gradUx,gradUy,gradUz
!> Physical fluxes in x,y,z directions
REAL,DIMENSION(CONS,0:PP_N_xi,0:PP_N_eta,0:PP_N_zeta),INTENT(OUT)          :: f,g,h
INTEGER,INTENT(IN)                                                         :: iElem                !< element index in global array
INTEGER,INTENT(IN)                                                         :: PP_N_xi
INTEGER,INTENT(IN)                                                         :: PP_N_eta
INTEGER,INTENT(IN)                                                         :: PP_N_zeta
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k
!==================================================================================================================================
DO k=0,PP_N_zeta;  DO j=0,PP_N_eta; DO i=0,PP_N_xi
  CALL EvalDiffFlux3D_Point(Uprim(:,i,j,k),gradUx(:,i,j,k),gradUy(:,i,j,k),gradUz(:,i,j,k), &
                                                f(:,i,j,k),     g(:,i,j,k),     h(:,i,j,k)  &
#if EDDYVISCOSITY
                            ,muSGS(1,i,j,k,iElem)&
#endif
                            )
END DO; END DO; END DO ! i,j,k

#if !EDDYVISCOSITY
NO_OP(iElem)
#endif /*!EDDYVISCOSITY*/
END SUBROUTINE EvalDiffFlux3D_Volume_FV
#endif /*FV_ENABLED*/

!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single side
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D_Surface(Nloc,UPrim,gradUx,gradUy,gradUz,f,g,h &
#if EDDYVISCOSITY
                                 ,muSGS &
#endif /*EDDYVISCOSITY*/
                                 )
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER                                           ,INTENT(IN)  :: Nloc                 !< Polynomial degree of input solution
REAL,DIMENSION(PP_nVarPrim   ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim                !< Solution vector
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: gradUx,gradUy,gradUz !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar       ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: f,g,h                !> Physical fluxes in x,y,z directions
#if EDDYVISCOSITY
REAL,DIMENSION(1             ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: muSGS                !< SGS viscosity
#endif /*EDDYVISCOSITY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j
!==================================================================================================================================
DO j=0,ZDIM(Nloc); DO i=0,Nloc
  CALL EvalDiffFlux3D_Point(Uprim(:,i,j),gradUx(:,i,j),gradUy(:,i,j),gradUz(:,i,j), &
                                               f(:,i,j),     g(:,i,j),     h(:,i,j)  &
#if EDDYVISCOSITY
                            ,muSGS(1,i,j) &
#endif /*EDDYVISCOSITY*/
                            )
END DO; END DO ! i,j
END SUBROUTINE EvalDiffFlux3D_Surface
#endif /*PARABOLIC*/



!----------------------------------------------------------------------------------------------------------------------------------
! EvalEulerFlux1D_fast METHOD
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
!> Computes 1D Euler flux using the conservative and primitive variables (for better performance)
!==================================================================================================================================
PPURE SUBROUTINE EvalEulerFlux1D_Host(U,F)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: U(PP_2Var) !< vector of conservative and primitive variables
REAL,INTENT(OUT)    :: F(PP_nVar) !< Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

  ! Euler fluxes x-direction
  F(DENS)= U(EXT_MOM1)                         ! rho*u
  F(MOM1)= U(EXT_MOM1)*U(EXT_VEL1)+U(EXT_PRES) ! rho*u²+p
  F(MOM2)= U(EXT_MOM1)*U(EXT_VEL2)             ! rho*u*v
  F(MOM3)= U(EXT_MOM1)*U(EXT_VEL3)             ! rho*u*w
  F(ENER)=(U(EXT_ENER)+U(EXT_PRES))*U(EXT_VEL1)! (rho*e+p)*u

END SUBROUTINE EvalEulerFlux1D_Host

END MODULE MOD_Flux
