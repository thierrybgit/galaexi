!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! Copyright (c) 2016-2017 Gregor Gassner (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Florian Hindenlang (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Andrew Winters (github.com/project-fluxo/fluxo)
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

!==================================================================================================================================
!> Contains the routines for computing the fluxes for the Split-DG algorithm
!> Attention: The volume split fluxes have to incorporate a factor of 2 coming from the splitDG formulation. This can be used to
!> cancel out a single factor of 1/2 in the averages, thus reducing the required operations.
!> See the given references for details of the split fluxes or e.g. Gassner, Gregor J., Andrew R. Winters, and David A. Kopriva.
!> "Split form nodal discontinuous Galerkin schemes with summation-by-parts property for the compressible Euler equations."
!> Journal of Computational Physics 327 (2016): 39-66. for an overview.
!==================================================================================================================================
#include "flexi.h"
#include "eos.h"

MODULE MOD_SplitFlux
! MODULES
USE ISO_C_BINDING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
ABSTRACT INTERFACE
  PPURE SUBROUTINE VolumeFlux(URef,UPrimRef,U,UPrim,MRef,M,Flux)
    REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U
    REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim
    REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,M
    REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  PPURE SUBROUTINE SurfaceFlux(U_LL,U_RR,F)
    REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR
    REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F
  END SUBROUTINE
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
! DEVICE METHOD INTERFACES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE
  SUBROUTINE InitSplitDG_Device(splitDGOpt) BIND(C, NAME="InitSplitDG_Device")
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(C_INT),VALUE :: splitDGOpt
  END SUBROUTINE InitSplitDG_Device
END INTERFACE

PROCEDURE(VolumeFlux),POINTER    :: SplitDGVolume_pointer    !< pointer defining the SpliDG formulation beeing used
PROCEDURE(SurfaceFlux),POINTER   :: SplitDGSurface_pointer   !< pointer defining the SpliDG formulation beeing used

INTERFACE InitSplitDG
  MODULE PROCEDURE InitSplitDG
END INTERFACE

PUBLIC::InitSplitDG,DefineParametersSplitDG
PUBLIC::SplitDGSurface_pointer,SplitDGVolume_pointer
PUBLIC::GetLogMean
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersSplitDG()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("SplitDG")
CALL prms%CreateIntFromStringOption('SplitDG',"SplitDG formulation to be used: SD, MO, DU, KG, PI, CH","PI")
CALL addStrListEntry('SplitDG','sd',           PRM_SPLITDG_SD)
CALL addStrListEntry('SplitDG','mo',           PRM_SPLITDG_MO)
CALL addStrListEntry('SplitDG','du',           PRM_SPLITDG_DU)
CALL addStrListEntry('SplitDG','kg',           PRM_SPLITDG_KG)
CALL addStrListEntry('SplitDG','pi',           PRM_SPLITDG_PI)
CALL addStrListEntry('SplitDG','ch',           PRM_SPLITDG_CH)

END SUBROUTINE DefineParametersSplitDG

!==================================================================================================================================!
!> Initialize function pointers for the specific split version in use
!==================================================================================================================================!
SUBROUTINE InitSplitDG(SplitDG_in)
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR
USE MOD_DG_Vars     ,ONLY: SplitDG
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL :: SplitDG_in
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: dev_err
!==================================================================================================================================
! set pointers
IF (PRESENT(SplitDG_in)) THEN
  SplitDG = SplitDG_in
ELSE
  SplitDG = GETINTFROMSTR('SplitDG')
END IF
SELECT CASE(SplitDG)
CASE(PRM_SPLITDG_SD)
  SplitDGVolume_pointer  => SplitVolumeFluxSD_Host
  SplitDGSurface_pointer => SplitSurfaceFluxSD_Host
CASE(PRM_SPLITDG_MO)
  SplitDGVolume_pointer  => SplitVolumeFluxMO_Host
  SplitDGSurface_pointer => SplitSurfaceFluxMO_Host
CASE(PRM_SPLITDG_DU)
  SplitDGVolume_pointer  => SplitVolumeFluxDU_Host
  SplitDGSurface_pointer => SplitSurfaceFluxDU_Host
CASE(PRM_SPLITDG_KG)
  SplitDGVolume_pointer  => SplitVolumeFluxKG_Host
  SplitDGSurface_pointer => SplitSurfaceFluxKG_Host
CASE(PRM_SPLITDG_PI)
  SplitDGVolume_pointer  => SplitVolumeFluxPI_Host
  SplitDGSurface_pointer => SplitSurfaceFluxPI_Host
CASE(PRM_SPLITDG_CH)
  SplitDGVolume_pointer  => SplitVolumeFluxCH_Host
  SplitDGSurface_pointer => SplitSurfaceFluxCH_Host
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'SplitDG formulation not defined!')
END SELECT

! Initialize the Split DG func pointers on the device side as well
#if (USE_ACCEL != ACCEL_OFF)
CALL InitSplitDG_Device(SplitDG)
#endif

END SUBROUTINE InitSplitDG

!==================================================================================================================================
!> Computes the Split-Flux retaining the standard NS-Equations
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxSD_Host(URef,UPrimRef,U,UPrim,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(CONS),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(CONS),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(1:3 ),INTENT(IN)  :: MRef          !< metric terms
REAL,DIMENSION(1:3 ),INTENT(IN)  :: M             !< metric terms
REAL,DIMENSION(CONS),INTENT(OUT) :: Flux          !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: rhoEpRef,rhoEp ! auxiliary variable for (rho*E+p)
REAL,DIMENSION(PP_nVar)                 :: fTilde,gTilde  ! flux in physical space
REAL,DIMENSION(PP_nVar)                 :: hTilde         ! flux in physical space
!==================================================================================================================================

  ! compute auxiliary variables, total energy density plus pressure
  rhoEpRef = URef(ENER) + UPrimRef(PRES)
  rhoEp    = U(ENER) + UPrim(PRES)

  ! local Euler fluxes x-direction
  fTilde(DENS) = (URef(MOM1) + U(MOM1))                                                       ! {rho*u}
  fTilde(MOM1) = (URef(MOM1)*UPrimRef(VEL1)+UPrimRef(PRES) + U(MOM1)*UPrim(VEL1)+UPrim(PRES)) ! {rho*u²}+{p}
  fTilde(MOM2) = (URef(MOM1)*UPrimRef(VEL2) + U(MOM1)*UPrim(VEL2))                            ! {rho*u*v}
  fTilde(MOM3) = (URef(MOM1)*UPrimRef(VEL3) + U(MOM1)*UPrim(VEL3))                            ! {rho*u*w}
  fTilde(ENER) = (rhoEpRef*UPrimRef(VEL1) + rhoEp*UPrim(VEL1))                                ! {(rho*E+p)*u}
  ! local Euler fluxes y-direction
  gTilde(DENS) = (URef(MOM2) + U(MOM2))                                                       ! {rho*v}
  gTilde(MOM1) = (URef(MOM1)*UPrimRef(VEL2) + U(MOM1)*UPrim(VEL2))                            ! {rho*u*v}
  gTilde(MOM2) = (URef(MOM2)*UPrimRef(VEL2)+UPrimRef(PRES) + U(MOM2)*UPrim(VEL2)+UPrim(PRES)) ! {rho*v²}+{p}
  gTilde(MOM3) = (URef(MOM2)*UPrimRef(VEL3) + U(MOM2)*UPrim(VEL3))                            ! {rho*v*w}
  gTilde(ENER) = (rhoEpRef*UPrimRef(VEL2) + rhoEp*UPrim(VEL2))                                ! {(rho*E+p)*v}
  ! local Euler fluxes z-direction
  hTilde(DENS) = (URef(MOM3) + U(MOM3))                                                       ! {rho*w}
  hTilde(MOM1) = (URef(MOM1)*UPrimRef(VEL3) + U(MOM1)*UPrim(VEL3))                            ! {rho*u*w}
  hTilde(MOM2) = (URef(MOM2)*UPrimRef(VEL3) + U(MOM2)*UPrim(VEL3))                            ! {rho*v*w}
  hTilde(MOM3) = (URef(MOM3)*UPrimRef(VEL3)+UPrimRef(PRES) + U(MOM3)*UPrim(VEL3)+UPrim(PRES)) ! {rho*v²+p}
  hTilde(ENER) = (rhoEpRef*UPrimRef(VEL3) + rhoEp*UPrim(VEL3))                                ! {(rho*E+p)*w}

  ! transform into reference space
  Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
            0.5*(MRef(3)+M(3))*hTilde(:) + &
            0.5*(MRef(2)+M(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxSD_Host

!==================================================================================================================================
!> Computes the surface flux for the split formulation retaining the standard NS-Equations
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxSD_Host(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL      !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR      !< variables at the right surfaces
REAL,DIMENSION(CONS   ),INTENT(OUT) :: F         !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

  F(DENS)= 0.5*(U_LL(EXT_MOM1)+U_RR(EXT_MOM1))                                                                 ! {rho*u}
  F(MOM1)= 0.5*(U_LL(EXT_MOM1)*U_LL(EXT_VEL1)+U_LL(EXT_PRES)+U_RR(EXT_MOM1)*U_RR(EXT_VEL1)+U_RR(EXT_PRES))     ! {rho*u²}+{p}
  F(MOM2)= 0.5*(U_LL(EXT_MOM1)*U_LL(EXT_VEL2)+U_RR(EXT_MOM1)*U_RR(EXT_VEL2))                                   ! {rho*u*v}
  F(MOM3)= 0.5*(U_LL(EXT_MOM1)*U_LL(EXT_VEL3)+U_RR(EXT_MOM1)*U_RR(EXT_VEL3))                                   ! {rho*u*w}
  F(ENER)= 0.5*((U_LL(EXT_ENER)+U_LL(EXT_PRES))*U_LL(EXT_VEL1)+(U_RR(EXT_ENER)+U_RR(EXT_PRES))*U_RR(EXT_VEL1)) ! {(rho*E+p)*u}

END SUBROUTINE SplitSurfaceFluxSD_Host

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Ducros
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!> Uses quadratic split forms in all equations. Reference: Ducros, F., et al. "High-order fluxes for conservative
!> skew-symmetric-like schemes in structured meshes: application to compressible flows."
!> Journal of Computational Physics 161.1 (2000): 114-139.
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxDU_Host(URef,UPrimRef,U,UPrim,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(CONS),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(CONS),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(1:3 ),INTENT(IN)  :: MRef          !< metric terms
REAL,DIMENSION(1:3 ),INTENT(IN)  :: M             !< metric terms
REAL,DIMENSION(CONS),INTENT(OUT) :: Flux          !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar)             :: fTilde,gTilde     ! flux in physical space
REAL,DIMENSION(PP_nVar)             :: hTilde            ! flux in physical space
!==================================================================================================================================

  ! local Euler fluxes x-direction
  fTilde(DENS) = 0.5*(URef(DENS)+U(DENS))*(UPrimRef(VEL1)+UPrim(VEL1))                                ! {rho}*{u}
  fTilde(MOM1) = 0.5*(URef(MOM1)+U(MOM1))*(UPrimRef(VEL1)+UPrim(VEL1)) + (UPrimRef(PRES)+UPrim(PRES)) ! {rho*u}*{u}+{p}
  fTilde(MOM2) = 0.5*(URef(MOM2)+U(MOM2))*(UPrimRef(VEL1)+UPrim(VEL1))                                ! {rho*v}*{u}
  fTilde(MOM3) = 0.5*(URef(MOM3)+U(MOM3))*(UPrimRef(VEL1)+UPrim(VEL1))                                ! {rho*w}*{u}
  fTilde(ENER) = 0.5*(URef(ENER)+U(ENER)+UPrimRef(PRES)+UPrim(PRES))*(UPrimRef(VEL1)+UPrim(VEL1))     ! ({rho*E}+{p})*{u}
  ! local Euler fluxes y-direction
  gTilde(DENS) = 0.5*(URef(DENS)+U(DENS))*(UPrimRef(VEL2)+UPrim(VEL2))                                ! {rho}*{v}
  gTilde(MOM1) = 0.5*(URef(MOM1)+U(MOM1))*(UPrimRef(VEL2)+UPrim(VEL2))                                ! {rho*u}*{v}
  gTilde(MOM2) = 0.5*(URef(MOM2)+U(MOM2))*(UPrimRef(VEL2)+UPrim(VEL2)) + (UPrimRef(PRES)+UPrim(PRES)) ! {rho*v}*{v}+{p}
  gTilde(MOM3) = 0.5*(URef(MOM3)+U(MOM3))*(UPrimRef(VEL2)+UPrim(VEL2))                                ! {rho*w}*{v}
  gTilde(ENER) = 0.5*(URef(ENER)+U(ENER)+UPrimRef(PRES)+UPrim(PRES))*(UPrimRef(VEL2)+UPrim(VEL2))     ! ({rho*E}+{p})*{v}
  ! local Euler fluxes z-direction
  hTilde(DENS) = 0.5*(URef(DENS)+U(DENS))*(UPrimRef(VEL3)+UPrim(VEL3))                                ! {rho}*{w}
  hTilde(MOM1) = 0.5*(URef(MOM1)+U(MOM1))*(UPrimRef(VEL3)+UPrim(VEL3))                                ! {rho*u}*{w}
  hTilde(MOM2) = 0.5*(URef(MOM2)+U(MOM2))*(UPrimRef(VEL3)+UPrim(VEL3))                                ! {rho*v}*{w}
  hTilde(MOM3) = 0.5*(URef(MOM3)+U(MOM3))*(UPrimRef(VEL3)+UPrim(VEL3)) + (UPrimRef(PRES)+UPrim(PRES)) ! {rho*w}*{w}+{p}
  hTilde(ENER) = 0.5*(URef(ENER)+U(ENER)+UPrimRef(PRES)+UPrim(PRES))*(UPrimRef(VEL3)+UPrim(VEL3))     ! ({rho*E}+{p})*{w}

  ! transform into reference space
  Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
            0.5*(MRef(3)+M(3))*hTilde(:) + &
            0.5*(MRef(2)+M(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxDU_Host

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Ducros
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxDU_Host(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL      !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR      !< variables at the right surfaces
REAL,DIMENSION(CONS   ),INTENT(OUT) :: F         !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

  F(DENS)= 0.25*(U_LL(EXT_DENS)+U_RR(EXT_DENS))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))                                       ! {rho}*{u}
  F(MOM1)= 0.25*(U_LL(EXT_MOM1)+U_RR(EXT_MOM1))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1)) + 0.5*(U_LL(EXT_PRES)+U_RR(EXT_PRES)) ! {rho*u}*{u}+{p}
  F(MOM2)= 0.25*(U_LL(EXT_MOM2)+U_RR(EXT_MOM2))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))                                       ! {rho*v}*{u}
  F(MOM3)= 0.25*(U_LL(EXT_MOM3)+U_RR(EXT_MOM3))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))                                       ! {rho*w}*{u}
  F(ENER)= 0.25*(U_LL(EXT_ENER)+U_RR(EXT_ENER)+U_LL(EXT_PRES)+U_RR(EXT_PRES))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))         ! ({rho*E}+{p})*{u}

END SUBROUTINE SplitSurfaceFluxDU_Host

!==================================================================================================================================
!> Computes the Split-Flux retaining the KEP formulation of Kennedy and Gruber
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!> Reference: Kennedy, Christopher A., and Andrea Gruber. "Reduced aliasing formulations of the convective terms within the
!> Navier–Stokes equations for a compressible fluid." Journal of Computational Physics 227.3 (2008): 1676-1700.
!> Uses a quadratic splitting for u*p and a cubic splitting for rho*e*u in the energy equation.
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxKG_Host(URef,UPrimRef,U,UPrim,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(CONS),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(CONS),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(1:3 ),INTENT(IN)  :: MRef          !< metric terms
REAL,DIMENSION(1:3 ),INTENT(IN)  :: M             !< metric terms
REAL,DIMENSION(CONS),INTENT(OUT) :: Flux          !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: E,ERef        ! auxiliary variables for the specific total energy
REAL,DIMENSION(PP_nVar)                 :: fTilde,gTilde ! flux in physical space
REAL,DIMENSION(PP_nVar)                 :: hTilde        ! flux in physical space
!==================================================================================================================================

  ! specific total energy
  ERef = URef(ENER)/URef(DENS)
  E    = U(ENER)/U(DENS)

  ! local Euler fluxes x-direction
  fTilde(DENS) = 0.5* (URef(DENS)+U(DENS))*(UPrimRef(VEL1)+UPrim(VEL1))                                   ! {rho}*{u}
  fTilde(MOM1) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL1)+UPrim(VEL1))**2 + (UPrimRef(PRES)+UPrim(PRES)) ! {rho}*{u}²+{p}
  fTilde(MOM2) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL1)+UPrim(VEL1))*(UPrimRef(VEL2)+UPrim(VEL2))      ! {rho}*{u}*{v}
  fTilde(MOM3) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL1)+UPrim(VEL1))*(UPrimRef(VEL3)+UPrim(VEL3))      ! {rho}*{u}*{w}
  fTilde(ENER) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL1)+UPrim(VEL1))*(eRef+e) + &
                0.5* (UPrimRef(PRES)+UPrim(PRES))*(UPrimRef(VEL1)+UPrim(VEL1))                           ! {rho}*{E}*{u}+{p}*{u}
  ! local Euler fluxes y-direction
  gTilde(DENS) = 0.5 *(URef(DENS)+U(DENS))*(UPrimRef(VEL2)+UPrim(VEL2))                                   ! {rho}*{v}
  gTilde(MOM1) = fTilde(MOM2)                                                                             ! {rho}*{v}*{u}
  gTilde(MOM2) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL2)+UPrim(VEL2))**2 + (UPrimRef(PRES)+UPrim(PRES)) ! {rho}*{v}²+{p}
  gTilde(MOM3) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL2)+UPrim(VEL2))*(UPrimRef(VEL3)+UPrim(VEL3))      ! {rho}*{v}*{w}
  gTilde(ENER) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL2)+UPrim(VEL2))*(eRef+e) + &
                0.5* (UPrimRef(PRES)+UPrim(PRES))*(UPrimRef(VEL2)+UPrim(VEL2))                           ! {rho}*{E}*{v}+{p}*{v}
  ! local Euler fluxes z-direction
  hTilde(DENS) = 0.5 *(URef(DENS)+U(DENS))*(UPrimRef(VEL3)+UPrim(VEL3))                                   ! {rho}*{w}
  hTilde(MOM1) = fTilde(MOM3)                                                                             ! {rho}*{w}*{u}
  hTilde(MOM2) = gTilde(MOM3)                                                                             ! {rho}*{w}*{v}
  hTilde(MOM3) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL3)+UPrim(VEL3))**2 + (UPrimRef(PRES)+UPrim(PRES)) ! {rho}*{w}²+{p}
  hTilde(ENER) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL3)+UPrim(VEL3))*(eRef+e) + &
                0.5 *(UPrimRef(PRES)+UPrim(PRES))*(UPrimRef(VEL3)+UPrim(VEL3))                           ! {rho}*{E}*{w}+{p}*{w}

  ! transform into reference space
  Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
            0.5*(MRef(3)+M(3))*hTilde(:) + &
            0.5*(MRef(2)+M(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxKG_Host

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Kennedy and Gruber
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxKG_Host(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL      !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR      !< variables at the right surfaces
REAL,DIMENSION(CONS   ),INTENT(OUT) :: F         !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: E_LL,E_RR ! auxiliary variables for the specific total energy
!==================================================================================================================================

  ! specific total energy
  E_LL = U_LL(EXT_ENER)/U_LL(EXT_DENS)
  E_RR = U_RR(EXT_ENER)/U_RR(EXT_DENS)
  !compute flux
  F(DENS)= 0.25* (U_LL(EXT_DENS)+U_RR(EXT_DENS))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))                                          ! {rho}*{u}
  F(MOM1)= 0.125*(U_LL(EXT_DENS)+U_RR(EXT_DENS))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))**2 + 0.5*(U_LL(EXT_PRES)+U_RR(EXT_PRES)) ! {rho}*{u}²+{p}
  F(MOM2)= 0.125*(U_LL(EXT_DENS)+U_RR(EXT_DENS))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))*(U_LL(EXT_VEL2)+U_RR(EXT_VEL2))          ! {rho}*{u}*{v}
  F(MOM3)= 0.125*(U_LL(EXT_DENS)+U_RR(EXT_DENS))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))*(U_LL(EXT_VEL3)+U_RR(EXT_VEL3))          ! {rho}*{u}*{w}
  F(ENER)= 0.125*(U_LL(EXT_DENS)+U_RR(EXT_DENS))*(E_LL+E_RR)*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1)) + &
          0.25 *(U_LL(EXT_PRES)+U_RR(EXT_PRES))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))                                          ! {rho}*{E}*{u}+{p}*{u}

END SUBROUTINE SplitSurfaceFluxKG_Host

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Morinishi
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!> Uses quadratic splittings in the momentum equations, but pure advection form in the energy equation. This seems to lead to a
!> rather unstable formulation!
!> Reference: Morinishi, Yohei. "Skew-symmetric form of convective terms and fully conservative finite difference schemes for
!> variable density low-Mach number flows." Journal of Computational Physics 229.2 (2010): 276-300.
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxMO_Host(URef,UPrimRef,U,UPrim,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(CONS),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(CONS),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(1:3 ),INTENT(IN)  :: MRef          !< metric terms
REAL,DIMENSION(1:3 ),INTENT(IN)  :: M             !< metric terms
REAL,DIMENSION(CONS),INTENT(OUT) :: Flux          !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: rhoepRef,rhoep ! auxiliary variable for rho*inner energy + pressure
REAL,DIMENSION(PP_nVar)                 :: fTilde,gTilde  ! flux in physical space
REAL,DIMENSION(PP_nVar)                 :: hTilde         ! flux in physical space
!==================================================================================================================================

  ! rho*internal energy + pressure
  rhoepRef = URef(ENER)-0.5*URef(DENS)*(UPrimRef(VEL1)**2+UPrimRef(VEL2)**2+UPrimRef(VEL3)**2)+UPrimRef(PRES)
  rhoep    = U(ENER)-0.5*U(DENS)*(UPrim(VEL1)**2+UPrim(VEL2)**2+UPrim(VEL3)**2)+UPrim(PRES)

  ! local Euler fluxes x-direction
  fTilde(DENS) =     (URef(MOM1)+U(MOM1))                                                             ! {rho*u}
  fTilde(MOM1) = 0.5*(URef(MOM1)+U(MOM1))*(UPrimRef(VEL1)+UPrim(VEL1)) + (UPrimRef(PRES)+UPrim(PRES)) ! {rho*u}*{u}+{p}
  fTilde(MOM2) = 0.5*(URef(MOM1)+U(MOM1))*(UPrimRef(VEL2)+UPrim(VEL2))                                ! {rho*u}*{v}
  fTilde(MOM3) = 0.5*(URef(MOM1)+U(MOM1))*(UPrimRef(VEL3)+UPrim(VEL3))                                ! {rho*u}*{w}
  fTilde(ENER) = (rhoepRef*UPrimRef(VEL1)+rhoep*UPrim(VEL1)) + &                                      !{(rho*e+p)*u} +
              0.5*(URef(MOM1)*UPrimRef(VEL1)+U(MOM1)*UPrim(VEL1))*(UPrimRef(VEL1)+UPrim(VEL1)) + &    !{rho*u²}*{u} +
              0.5*(URef(MOM1)*UPrimRef(VEL2)+U(MOM1)*UPrim(VEL2))*(UPrimRef(VEL2)+UPrim(VEL2)) + &    !{rho*u*v}*{v} +
              0.5*(URef(MOM1)*UPrimRef(VEL3)+U(MOM1)*UPrim(VEL3))*(UPrimRef(VEL3)+UPrim(VEL3)) - &    !{rho*u*w}*{w} -
              0.5*(URef(MOM1)*UPrimRef(VEL1)*UPrimRef(VEL1)+U(MOM1)*UPrim(VEL1)*UPrim(VEL1))   - &    !1/2*({rho*u³} +
              0.5*(URef(MOM1)*UPrimRef(VEL2)*UPrimRef(VEL2)+U(MOM1)*UPrim(VEL2)*UPrim(VEL2))   - &    !{rho*u*v²} +
              0.5*(URef(MOM1)*UPrimRef(VEL3)*UPrimRef(VEL3)+U(MOM1)*UPrim(VEL3)*UPrim(VEL3))          !{rho*u*w²})
  ! local Euler fluxes y-direction
  gTilde(DENS) =     (URef(MOM2)+U(MOM2))                                                             ! {rho*v}
  gTilde(MOM1) = 0.5*(URef(MOM2)+U(MOM2))*(UPrimRef(VEL1)+UPrim(VEL1))                                ! {rho*v}*{u}
  gTilde(MOM2) = 0.5*(URef(MOM2)+U(MOM2))*(UPrimRef(VEL2)+UPrim(VEL2)) + (UPrimRef(PRES)+UPrim(PRES)) ! {rho*v}*{v}+{p}
  gTilde(MOM3) = 0.5*(URef(MOM2)+U(MOM2))*(UPrimRef(VEL3)+UPrim(VEL3))                                ! {rho*v}*{w}
  gTilde(ENER) = (rhoepRef*UPrimRef(VEL2)+rhoep*UPrim(VEL2)) + &                                      !{(rho*e+p)*v} +
              0.5*(URef(MOM2)*UPrimRef(VEL1)+U(MOM2)*UPrim(VEL1))*(UPrimRef(VEL1)+UPrim(VEL1)) + &    !{rho*v*u}*{u} +
              0.5*(URef(MOM2)*UPrimRef(VEL2)+U(MOM2)*UPrim(VEL2))*(UPrimRef(VEL2)+UPrim(VEL2)) + &    !{rho*v²}*{v} +
              0.5*(URef(MOM2)*UPrimRef(VEL3)+U(MOM2)*UPrim(VEL3))*(UPrimRef(VEL3)+UPrim(VEL3)) - &    !{rho*v*w}*{w} -
              0.5*(URef(MOM2)*UPrimRef(VEL1)*UPrimRef(VEL1)+U(MOM2)*UPrim(VEL1)*UPrim(VEL1))   - &    !1/2*({rho*v*u²} +
              0.5*(URef(MOM2)*UPrimRef(VEL2)*UPrimRef(VEL2)+U(MOM2)*UPrim(VEL2)*UPrim(VEL2))   - &    !{rho*v³} +
              0.5*(URef(MOM2)*UPrimRef(VEL3)*UPrimRef(VEL3)+U(MOM2)*UPrim(VEL3)*UPrim(VEL3))          !{rho*v*w²})
  ! local Euler fluxes z-direction
  hTilde(DENS) =     (URef(MOM3)+U(MOM3))                                                             ! {rho*w}
  hTilde(MOM1) = 0.5*(URef(MOM3)+U(MOM3))*(UPrimRef(VEL1)+UPrim(VEL1))                                ! {rho*w}*{u}
  hTilde(MOM2) = 0.5*(URef(MOM3)+U(MOM3))*(UPrimRef(VEL2)+UPrim(VEL2))                                ! {rho*w}*{v}
  hTilde(MOM3) = 0.5*(URef(MOM3)+U(MOM3))*(UPrimRef(VEL3)+UPrim(VEL3)) + (UPrimRef(PRES)+UPrim(PRES)) ! {rho*w}*{w}+{p}
  hTilde(ENER) = (rhoepRef*UPrimRef(VEL3)+rhoep*UPrim(VEL3)) + &                                      !{(rho*e+p)*w} +
              0.5*(URef(MOM3)*UPrimRef(VEL1)+U(MOM3)*UPrim(VEL1))*(UPrimRef(VEL1)+UPrim(VEL1)) + &    !{rho*w*u}*{u} +
              0.5*(URef(MOM3)*UPrimRef(VEL2)+U(MOM3)*UPrim(VEL2))*(UPrimRef(VEL2)+UPrim(VEL2)) + &    !{rho*w*v}*{v} +
              0.5*(URef(MOM3)*UPrimRef(VEL3)+U(MOM3)*UPrim(VEL3))*(UPrimRef(VEL3)+UPrim(VEL3)) - &    !{rho*w²}*{w} -
              0.5*(URef(MOM3)*UPrimRef(VEL1)*UPrimRef(VEL1)+U(MOM3)*UPrim(VEL1)*UPrim(VEL1))   - &    !1/2*({rho*w*u²} +
              0.5*(URef(MOM3)*UPrimRef(VEL2)*UPrimRef(VEL2)+U(MOM3)*UPrim(VEL2)*UPrim(VEL2))   - &    !{rho*w*v²} +
              0.5*(URef(MOM3)*UPrimRef(VEL3)*UPrimRef(VEL3)+U(MOM3)*UPrim(VEL3)*UPrim(VEL3))          !{rho*w³})

  ! transform into reference space
  Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
            0.5*(MRef(3)+M(3))*hTilde(:) + &
            0.5*(MRef(2)+M(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxMO_Host

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Morinishi
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxMO_Host(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR !< variables at the right surfaces
REAL,DIMENSION(CONS   ),INTENT(OUT) :: F    !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: rhoep_LL,rhoep_RR
!==================================================================================================================================

  ! rho*internal energy + pressure
  rhoep_LL = U_LL(EXT_ENER)-0.5*U_LL(EXT_DENS)*(U_LL(EXT_VEL1)**2+U_LL(EXT_VEL2)**2+U_LL(EXT_VEL3)**2)+U_LL(EXT_PRES)
  rhoep_RR = U_RR(EXT_ENER)-0.5*U_RR(EXT_DENS)*(U_RR(EXT_VEL1)**2+U_RR(EXT_VEL2)**2+U_RR(EXT_VEL3)**2)+U_RR(EXT_PRES)

  ! compute flux
  F(DENS)= 0.5 *(U_LL(EXT_MOM1)+U_RR(EXT_MOM1))                                                                       ! {rho*u}
  F(MOM1)= 0.25*(U_LL(EXT_MOM1)+U_RR(EXT_MOM1))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1)) + 0.5*(U_LL(EXT_PRES)+U_RR(EXT_PRES)) ! {rho*u}*{u}+{p}
  F(MOM2)= 0.25*(U_LL(EXT_MOM1)+U_RR(EXT_MOM1))*(U_LL(EXT_VEL2)+U_RR(EXT_VEL2))                                       ! {rho*u}*{v}
  F(MOM3)= 0.25*(U_LL(EXT_MOM1)+U_RR(EXT_MOM1))*(U_LL(EXT_VEL3)+U_RR(EXT_VEL3))                                       ! {rho*u}*{w}
  F(ENER)= 0.5 *(rhoep_LL*U_LL(EXT_VEL1)+rhoep_RR*U_RR(EXT_VEL1)) +  &                                                !{(rho*e+p)*u} +
          0.25*(U_LL(EXT_MOM1)*U_LL(EXT_VEL1)+U_RR(EXT_MOM1)*U_RR(EXT_VEL1))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1)) + &     !{rho*u²}*{u} +
          0.25*(U_LL(EXT_MOM1)*U_LL(EXT_VEL2)+U_RR(EXT_MOM1)*U_RR(EXT_VEL2))*(U_LL(EXT_VEL2)+U_RR(EXT_VEL2)) + &     !{rho*u*v}*{v} +
          0.25*(U_LL(EXT_MOM1)*U_LL(EXT_VEL3)+U_RR(EXT_MOM1)*U_RR(EXT_VEL3))*(U_LL(EXT_VEL3)+U_RR(EXT_VEL3)) - &     !{rho*u*w}*{w} -
          0.25*(U_LL(EXT_MOM1)*U_LL(EXT_VEL1)*U_LL(EXT_VEL1)+U_RR(EXT_MOM1)*U_RR(EXT_VEL1)*U_RR(EXT_VEL1)) - &       !1/2*({rho*u³} -
          0.25*(U_LL(EXT_MOM1)*U_LL(EXT_VEL2)*U_LL(EXT_VEL2)+U_RR(EXT_MOM1)*U_RR(EXT_VEL2)*U_RR(EXT_VEL2)) - &       !{rho*u*v²} -
          0.25*(U_LL(EXT_MOM1)*U_LL(EXT_VEL3)*U_LL(EXT_VEL3)+U_RR(EXT_MOM1)*U_RR(EXT_VEL3)*U_RR(EXT_VEL3))           !{rho*u*w²})

END SUBROUTINE SplitSurfaceFluxMO_Host

!==================================================================================================================================
!> Computes the Split-Flux retaining the KEP formulation of Pirozzoli
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!> This is a slight variation of the form proposed by Kennedy and Gruber, and is based on a cubic splitting of rho*H*u in the energy
!> equation. It was presented in the (worth reading) overview: Pirozzoli, Sergio. "Numerical methods for high-speed flows."
!> Annual review of fluid mechanics 43 (2011): 163-194.
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxPI_Host(URef,UPrimRef,U,UPrim,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(CONS),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(CONS),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(1:3 ),INTENT(IN)  :: MRef          !< metric terms
REAL,DIMENSION(1:3 ),INTENT(IN)  :: M             !< metric terms
REAL,DIMENSION(CONS),INTENT(OUT) :: Flux          !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: H,HRef        ! auxiliary variables for the specific enthalpy
REAL,DIMENSION(PP_nVar)                 :: fTilde,gTilde ! flux in physical space
REAL,DIMENSION(PP_nVar)                 :: hTilde        ! flux in physical space
!==================================================================================================================================

  ! specific enthalpy, H=E+p/rho=(rhoE+p)/rho
  HRef = (URef(ENER)+UPrimRef(PRES))/URef(DENS)
  H    = (U(ENER)+UPrim(PRES))/U(DENS)

  ! local Euler fluxes x-direction
  fTilde(DENS) = 0.5 *(URef(DENS)+U(DENS))*(UPrimRef(VEL1)+UPrim(VEL1))                                   ! {rho}*{u}
  fTilde(MOM1) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL1)+UPrim(VEL1))**2 + (UPrimRef(PRES)+UPrim(PRES)) ! {rho}*{u}²+{p}
  fTilde(MOM2) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL1)+UPrim(VEL1))*(UPrimRef(VEL2)+UPrim(VEL2))      ! {rho}*{u}*{v}
  fTilde(MOM3) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL1)+UPrim(VEL1))*(UPrimRef(VEL3)+UPrim(VEL3))      ! {rho}*{u}*{w}
  fTilde(ENER) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL1)+UPrim(VEL1))*(HRef+H)                          ! {rho}*{H}*{u}
  ! local Euler fluxes y-direction
  gTilde(DENS) = 0.5 *(URef(DENS)+U(DENS))*(UPrimRef(VEL2)+UPrim(VEL2))                                   ! {rho}*{v}
  gTilde(MOM1) = fTilde(MOM2)                                                                             ! {rho}*{v}*{u}
  gTilde(MOM2) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL2)+UPrim(VEL2))**2 + (UPrimRef(PRES)+UPrim(PRES)) ! {rho}*{v}²+{p}
  gTilde(MOM3) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL2)+UPrim(VEL2))*(UPrimRef(VEL3)+UPrim(VEL3))      ! {rho}*{v}*{w}
  gTilde(ENER) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL2)+UPrim(VEL2))*(HRef+H)                          ! {rho}*{H}*{v}
  ! local Euler fluxes z-direction
  hTilde(DENS) = 0.5 *(URef(DENS)+U(DENS))*(UPrimRef(VEL3)+UPrim(VEL3))                                   ! {rho}*{w}
  hTilde(MOM1) = fTilde(MOM3)                                                                             ! {rho}*{w}*{u}
  hTilde(MOM2) = gTilde(MOM3)                                                                             ! {rho}*{w}*{v}
  hTilde(MOM3) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL3)+UPrim(VEL3))**2 + (UPrimRef(PRES)+UPrim(PRES)) ! {rho}*{w}²+{p}
  hTilde(ENER) = 0.25*(URef(DENS)+U(DENS))*(UPrimRef(VEL3)+UPrim(VEL3))*(HRef+H)                          ! {rho}*{H}*{w}

  ! transform into reference space
  Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
            0.5*(MRef(3)+M(3))*hTilde(:) + &
            0.5*(MRef(2)+M(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxPI_Host

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Pirozzoli
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxPI_Host(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL      !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR      !< variables at the right surfaces
REAL,DIMENSION(CONS   ),INTENT(OUT) :: F         !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: H_LL,H_RR ! auxiliary variables for the specific energy
!==================================================================================================================================

  ! specific energy, H=E+p/rho=(rhoE+p)/rho
  H_LL = (U_LL(EXT_ENER)+U_LL(EXT_PRES))/U_LL(EXT_DENS)
  H_RR = (U_RR(EXT_ENER)+U_RR(EXT_PRES))/U_RR(EXT_DENS)
  !compute flux
  F(DENS)= 0.25* (U_LL(EXT_DENS)+U_RR(EXT_DENS))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))                                          ! {rho}*{u}
  F(MOM1)= 0.125*(U_LL(EXT_DENS)+U_RR(EXT_DENS))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))**2 + 0.5*(U_LL(EXT_PRES)+U_RR(EXT_PRES)) ! {rho}*{u}²+{p}
  F(MOM2)= 0.125*(U_LL(EXT_DENS)+U_RR(EXT_DENS))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))*(U_LL(EXT_VEL2)+U_RR(EXT_VEL2))          ! {rho}*{u}*{v}
  F(MOM3)= 0.125*(U_LL(EXT_DENS)+U_RR(EXT_DENS))*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))*(U_LL(EXT_VEL3)+U_RR(EXT_VEL3))          ! {rho}*{u}*{w}
  F(ENER)= 0.125*(U_LL(EXT_DENS)+U_RR(EXT_DENS))*(H_LL+H_RR)*(U_LL(EXT_VEL1)+U_RR(EXT_VEL1))                              ! {rho}*{H}*{u}

END SUBROUTINE SplitSurfaceFluxPI_Host

!==================================================================================================================================
!> Computes the Split-Flux retaining the entropy conserving (and formally KEP) formulation of Chandrashekar
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!> The flux after Chanrashekar uses a special computation of the pressure, based on the averages of density and inverse
!> temperature, which correspondonds to using the harmonic average of the temperature when applying the ideal gas law.
!> Reference: Chandrashekar, Praveen. "Kinetic energy preserving and entropy stable finite volume schemes for compressible Euler
!> and Navier-Stokes equations." Communications in Computational Physics 14.5 (2013): 1252-1286.
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxCH_Host(URef,UPrimRef,U,UPrim,MRef,M,Flux)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars, ONLY:sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(CONS),INTENT(IN)  :: URef     !< conserved variables
REAL,DIMENSION(CONS),INTENT(IN)  :: U        !< conserved variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrimRef !< primitive variables
REAL,DIMENSION(PRIM),INTENT(IN)  :: UPrim    !< primitive variables
REAL,DIMENSION(1:3 ),INTENT(IN)  :: MRef     !< metric terms
REAL,DIMENSION(1:3 ),INTENT(IN)  :: M        !< metric terms
REAL,DIMENSION(CONS),INTENT(OUT) :: Flux     !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: beta,betaRef            ! auxiliary variables for the inverse Temperature
REAL                                    :: pHatMean,HMean          ! auxiliary variable for the mean pressure and specific enthalpy
REAL                                    :: uMean,vMean,wMean       ! auxiliary variable for the average velocities
REAL                                    :: rhoLogMean,betaLogMean  ! auxiliary variable for the logarithmic means
REAL,DIMENSION(PP_nVar)                 :: fTilde,gTilde           ! flux in physical space
REAL,DIMENSION(PP_nVar)                 :: hTilde                  ! flux in physical space
!==================================================================================================================================

  ! average velocities
  uMean = 0.5*(UPrimRef(VEL1) + UPrim(VEL1))
  vMean = 0.5*(UPrimRef(VEL2) + UPrim(VEL2))
  wMean = 0.5*(UPrimRef(VEL3) + UPrim(VEL3))

  ! inverse temperature
  betaRef  = 0.5*URef(DENS)/UPrimRef(PRES)
  beta     = 0.5*U(DENS)/UPrim(PRES)

  ! Density and inverse temperature logarithmic average
  CALL GetLogMean(URef(DENS),U(DENS),rhoLogMean)
  CALL GetLogMean(betaRef,beta,betaLogMean)
  ! Average of pressure and specific enthalpy
  pHatMean = 0.5*(URef(DENS)+U(DENS))/(betaRef+beta)
  HMean    = 0.5*sKappaM1/betaLogMean + pHatMean/rhoLogMean + &
            0.5*DOT_PRODUCT(UPrimRef(VELV),UPrim(VELV))

  ! local Euler fluxes x-direction
  fTilde(DENS) = rhoLogMean*uMean                                      ! {rho}_log*{u}
  fTilde(MOM1) = rhoLogMean*uMean**2 + pHatMean                        ! {rho}_log*{u}²+{pHat}
  fTilde(MOM2) = rhoLogMean*uMean*vMean                                ! {rho}_log*{u}*{v}
  fTilde(MOM3) = rhoLogMean*uMean*wMean                                ! {rho}_log*{u}*{w}
  fTilde(ENER) = rhoLogMean*HMean*uMean                                ! {rho}_log*{H}*{u}
  ! local Euler fluxes y-direction
  gTilde(DENS) = rhoLogMean*vMean                                      ! {rho}_log*{v}
  gTilde(MOM1) = rhoLogMean*vMean*uMean                                ! {rho}_log*{v}*{u}
  gTilde(MOM2) = rhoLogMean*vMean**2 +pHatMean                         ! {rho}_log*{v}²+{pHat}
  gTilde(MOM3) = rhoLogMean*vMean*wMean                                ! {rho}_log*{v}*{w}
  gTilde(ENER) = rhoLogMean*HMean*vMean                                ! {rho}_log*{H}*{v}
  ! local Euler fluxes z-direction
  hTilde(DENS) = rhoLogMean*wMean                                      ! {rho}_log*{w}
  hTilde(MOM1) = rhoLogMean*wMean*uMean                                ! {rho}_log*{w}*{u}
  hTilde(MOM2) = rhoLogMean*wMean*vMean                                ! {rho}_log*{w}*{v}
  hTilde(MOM3) = rhoLogMean*wMean**2 + pHatMean                        ! {rho}_log*{w}²+{pHat}
  hTilde(ENER) = rhoLogMean*HMean*wMean                                ! {rho}_log*{H}*{w}

  ! transform into reference space
  Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
            0.5*(MRef(3)+M(3))*hTilde(:) + &
            0.5*(MRef(2)+M(2))*gTilde(:)

  ! Acount for factor of 2
  Flux = Flux*2.

END SUBROUTINE SplitVolumeFluxCH_Host

!==================================================================================================================================
!> Computes the surface flux for the entropy conserving formulation of Chandrashekar.
!> The flux after Chanrashekar uses a special computation of the pressure, based on the averages of density and inverse
!> temperature, which correspondonds to using the harmonic average of the temperature when applying the ideal gas law.
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxCH_Host(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars, ONLY:sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR !< variables at the right surfaces
REAL,DIMENSION(CONS   ),INTENT(OUT) :: F    !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: beta_LL,beta_RR        ! auxiliary variables for the inverse Temperature
REAL                                :: pHatMean,HMean         ! auxiliary variable for the mean pressure and specific enthalpy
REAL                                :: uMean,vMean,wMean      ! auxiliary variable for the average velocities
REAL                                :: rhoLogMean,betaLogMean ! auxiliary variable for the logarithmic mean
!==================================================================================================================================

  ! average velocities
  uMean = 0.5*(U_LL(EXT_VEL1) + U_RR(EXT_VEL1))
  vMean = 0.5*(U_LL(EXT_VEL2) + U_RR(EXT_VEL2))
  wMean = 0.5*(U_LL(EXT_VEL3) + U_RR(EXT_VEL3))

  ! inverse temperature
  beta_LL = 0.5*U_LL(EXT_DENS)/U_LL(EXT_PRES)
  beta_RR = 0.5*U_RR(EXT_DENS)/U_RR(EXT_PRES)

  ! average pressure, enthalpy, density and inverse temperature
  ! logarithmic mean
  CALL GetLogMean(U_LL(EXT_DENS),U_RR(EXT_DENS),rhoLogMean)
  CALL GetLogMean(beta_LL,beta_RR,betaLogMean)
  ! "standard" average
  pHatMean = 0.5*(U_LL(EXT_DENS)+U_RR(EXT_DENS))/(beta_LL+beta_RR)
  HMean    = 0.5*sKappaM1/betaLogMean + pHatMean/rhoLogMean + &
            0.5*(U_LL(EXT_VEL1)*U_RR(EXT_VEL1) + U_LL(EXT_VEL2)*U_RR(EXT_VEL2) + U_LL(EXT_VEL3)*U_RR(EXT_VEL3))

  !compute flux
  F(DENS) = rhoLogMean*uMean
  F(MOM1) = F(DENS)*uMean + pHatMean
  F(MOM2) = F(DENS)*vMean
  F(MOM3) = F(DENS)*wMean
  F(ENER) = F(DENS)*HMean

END SUBROUTINE SplitSurfaceFluxCH_Host

!==================================================================================================================================
!> auxilary function for calculating the logarithmic mean numerically stable according to Ismail and Roe
!==================================================================================================================================
ELEMENTAL SUBROUTINE GetLogMean(U_L,U_R,UMean)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: U_L   !< variables at the left surfaces
REAL,INTENT(IN)  :: U_R   !< variables at the right surfaces
REAL,INTENT(OUT) :: UMean !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER   :: epsilon = 0.01
REAL             :: chi,f,u,N ! auxiliary variables
!==================================================================================================================================
chi = U_L/U_R
f = (chi-1)/(chi+1)
u = f*f

IF (u .LT. epsilon) THEN
  N = 1.0+u/3.0+u*u/5.0+u*u*u/7.0
ELSE
  N = log(chi)/(2.*f)
ENDIF

UMean = (U_L+U_R)/(2.*N)
END SUBROUTINE getLogMean

END MODULE MOD_SplitFlux
