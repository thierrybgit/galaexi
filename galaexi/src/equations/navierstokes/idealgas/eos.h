#define COPYRIGHT "!==============================================================================================================  \
! Copyright (c) 2010-2022  Prof. Claus-Dieter Munz                                                                                  \
! Copyright (c) 2022-2024 Prof. Andrea Beck                                                                                         \
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.     \
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/                                      \
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License                \
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.             \
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty              \
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.                 \
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.   \
!=================================================================================================================================  \
! Define variables for normal and extended state vector                                                                             \
! Normal   U(1:5)  with conservative variables                                                                                      \
! Extended U(1:11) with conservative and primitive variables"   

#define CONS 1:PP_nVar          /* all cons variables */
#define PRIM 1:PP_nVarPrim      /* all prim variables */

#define PP_2Var PP_nVar+PP_nVarPrim

#define DENS  1             /* density */
#define MOM1  2             /* momentum x */
#define MOM2  3             /* momentum y */
#define MOM3  4             /* momentum z */
#define MOMV  MOM1:MOM3     /* momentum vector */
#define MMV2  MOM1:1+PP_dim /* momentum vector */
#define ENER  5             /* energy */

#define VEL1  2             /* velocity x */
#define VEL2  3             /* velocity y */
#define VEL3  4             /* velocity z */
#define VELV  VEL1:VEL3     /* velocity range */
#define VLV2  VEL1:6+PP_dim /* velocity range */
#define PRES  5             /* pressure */
#define TEMP  6             /* temperature */
#define VELVTEMP (/VEL1,VEL2,VEL3,TEMP/) /* velocity range and temperature */

#define EOS_KAPPA    1  /* isentropic coefficient kappa */
#define EOS_R        2  /* specific gas constant R      */
#define EOS_KAPPAM1  3
#define EOS_SKAPPAM1 4
#define EOS_KAPPAP1  5
#define EOS_SKAPPAP1 6
#define EOS_CP       7
#define EOS_CV       8
#if PARABOLIC /* EOS Parameters */
#define EOS_PR       9  /* Prandtl number Pr            */
#define EOS_KAPPASPR 10
#define EOS_MU0      11  /* reference viscosity mu0      */
#if PP_VISC == 0
#define PP_nVarEOS 11
#elif PP_VISC == 1
#define EOS_TREF     12 /* Sutherland constant */
#define EOS_EXPOSUTH 13 /* Sutherland constant */
#define EOS_TS       14 /* Sutherland constant */
#define EOS_CSUTH    15 /* Sutherland constant */
#define PP_nVarEOS   15
#elif PP_VISC == 2
#define EOS_TREF     12 /* Sutherland constant */
#define EOS_EXPOSUTH 13 /* Sutherland constant */
#define PP_nVarEOS   13
#endif /*PP_VISC == 1*/
#else /* not PARABOLIC */
#define PP_nVarEOS 8
#endif

#define KAPPASPR_MAX_TIMESTEP_H()      (MAX(4./3.,KappasPr)) /*routines to compute physical quantities*/
#define THERMAL_CONDUCTIVITY_H(mu)     (mu*cp/Pr)
#define TOTAL_TEMPERATURE_H(T,Mach)    (T*(1+0.5*(kappa-1)*Mach**2))
#define TOTAL_PRESSURE_H(p,Mach)       (p/((1+0.5*(kappa-1)*Mach**2)**(-kappa/(kappa-1.))))
#define BETA_RIEMANN_H()               (SQRT(0.5*kappaM1/kappa))
#define ROEC_RIEMANN_H(RoeH,RoeVel)    (SQRT(kappaM1*(RoeH-0.5*DOT_PRODUCT(RoeVel,RoeVel))))
#define ALPHA2_RIEMANN_H(RoeH,RoeVel,Roec,Delta_U)     (kappaM1/(Roec*Roec) * (Delta_U(1)*(RoeH-RoeVel(1)*RoeVel(1)) - Delta_U(6) + RoeVel(1)*Delta_U(2)))

#define VELOCITY_H(U,sRho)             (U(MOMV)*sRho) /*routines to compute physical quantities from conservative variables or extended variables*/
#define SPEEDOFSOUND_H(p,sRho)         (SQRT(Kappa*p*sRho))
#define TOTALENERGY_H(U,sRho,Vel)      (U(ENER)/U(DENS))
#define TOTALENTHALPY_H(U,p,sRho)      ((U(ENER)+p)*sRho)
#define ENTROPY_H(U,T)                 (R*(sKappaM1*LOG(T)-LOG(U(DENS))))
#define TEMPERATURE_H(U)               ((U(ENER)-0.5*DOT_PRODUCT(U(MOMV),U(MOMV))/U(DENS))/(U(DENS)*cv))
#define ENTROPY_HE(UE)                 R*(sKappaM1*LOG(UE(EXT_TEMP)) - LOG(UE(EXT_DENS)))

#define VELOCITY_HE(UE)                (UE(EXT_MOMV)*UE(EXT_SRHO))
#define PRESSURE_HE(UE)                (KappaM1*(UE(EXT_ENER)-0.5*DOT_PRODUCT(UE(EXT_VELV),UE(EXT_MOMV))))
#define SPEEDOFSOUND_HE(UE)            (SQRT(Kappa*UE(EXT_PRES)*UE(EXT_SRHO)))
#define TOTALENERGY_HE(UE)             (UE(EXT_ENER)*UE(EXT_SRHO))
#define TOTALENTHALPY_HE(UE)           ((UE(EXT_ENER)+UE(EXT_PRES))*UE(EXT_SRHO))
#define TEMPERATURE_HE(UE)             (UE(EXT_PRES)*UE(EXT_SRHO)/R)
#define ENERGY_HE(UE)                  (sKappaM1*UE(EXT_PRES)+0.5*DOT_PRODUCT(UE(EXT_MOMV),UE(EXT_VELV)))

#if PP_VISC == 0
#define VISCOSITY_PRIM(U)              mu0
#elif PP_VISC == 1
#define VISCOSITY_PRIM(U)              muSuth(U(TEMP))
#elif PP_VISC == 2
#define VISCOSITY_PRIM(U)              mu0*U(TEMP)**ExpoSuth
#endif

#if PP_VISC == 0
#define VISCOSITY_TEMPERATURE(T)       mu0
#elif PP_VISC == 1
#define VISCOSITY_TEMPERATURE(T)       muSuth(T)
#elif PP_VISC == 2
#define VISCOSITY_TEMPERATURE(T)       mu0*T**ExpoSuth
#endif

#define EXT_CONS    1:PP_nVar                  /* all ext cons variables */
#define EXT_PRIM    PP_nVarPrim:PP_2Var        /* all ext prim variables */

#define EXT_DENS    DENS                       /* density */
#define EXT_MOM1    MOM1                       /* momentum x */
#define EXT_MOM2    MOM2                       /* momentum y */
#define EXT_MOM3    MOM3                       /* momentum z */
#define EXT_MOMV    MOMV                       /* momentum vector */
#define EXT_ENER    ENER                       /* energy */

#define EXT_SRHO    PP_nVar+DENS               /* specific volume (1./density) */
#define EXT_VEL1    PP_nVar+VEL1               /* velocity x */
#define EXT_VEL2    PP_nVar+VEL2               /* velocity y */
#define EXT_VEL3    PP_nVar+VEL3               /* velocity z */
#define EXT_VELV    EXT_VEL1:EXT_VEL3          /* velocity range */
#define EXT_PRES    PP_nVar+PRES               /* pressure */
#define EXT_TEMP    PP_nVar+TEMP               /* temperature */

#if PP_OPTLIFT == 0 /*lifting variables*/
#define PP_nVarLifting               5
#define LIFT_DENS                    1
#define LIFT_VEL1                    2
#define LIFT_VEL2                    3
#define LIFT_VEL3                    4
#define LIFT_TEMP                    5
#define LIFT_VELV                    LIFT_VEL1:LIFT_VEL3
#define LIFT_VARS                    (/LIFT_DENS,LIFT_VEL1,LIFT_VEL2,LIFT_VEL3,LIFT_TEMP/)
#define PRIM_LIFT                    (/1,2,3,4,6/) /* density velocity range pressure and temperature */
#else
#define PP_nVarLifting               4
#define LIFT_VEL1                    1
#define LIFT_VEL2                    2
#define LIFT_VEL3                    3
#define LIFT_TEMP                    4
#define LIFT_VELV                    LIFT_VEL1:LIFT_VEL3
#define LIFT_VARS                    (/LIFT_VEL1,LIFT_VEL2,LIFT_VEL3,LIFT_TEMP/)
#define PRIM_LIFT                    (/2,3,4,6/) /* velocity range and temperature */
#endif

#define DELTA_U1                     1 /*Riemann Differences*/
#define DELTA_U2                     2
#define DELTA_U3                     3
#define DELTA_U4                     4
#define DELTA_UV                     DELTA_U2:DELTA_U4
#define DELTA_U5                     5
#define DELTA_U6                     6

#undef COPYRIGHT