/*================================================================================================================================
  Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
  Copyright (c) 2022-2024 Prof. Andrea Beck
  This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
  For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
 
  FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 
  FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
 
  You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
================================================================================================================================*/

//----------------------------------------------------------------------------------------------------------------------------------
// External C prototypes for procedures that bind to Fortran
//----------------------------------------------------------------------------------------------------------------------------------
extern "C"
{
    void InitRiemann_Device(int riemannOpt, int riemannOptBC);
}

//----------------------------------------------------------------------------------------------------------------------------------
// Riemann solver function pointers
//----------------------------------------------------------------------------------------------------------------------------------
/// Function pointer for Riemann solver
typedef void (*RiemannPointer_t)(double* F_L, double* F_R, double* U_LL, double* U_RR, double* riemannFlux);
__device__ RiemannPointer_t Riemann_Pointer_Device;

/// Function pointer to Riemann solver on boundaries
__device__ RiemannPointer_t RiemannBC_Pointer_Device;

//----------------------------------------------------------------------------------------------------------------------------------
// Advective flux methods
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Finds the Riemann flux for a single DOF
 **********************************************************************************************************************************/
__device__ void Riemann_Point_Device(double* Fout, double* U_L, double* U_R, double* UPrim_L, double* UPrim_R, 
                                                double* nv, double* t1, double* t2, bool doBC)
{
    double F_L[PP_nVar], F_R[PP_nVar], riemannFlux[PP_nVar];
    double U_LL[PP_2Var], U_RR[PP_2Var];

    // Momentum has to be rotated using the normal system individual for each
    // left state: U_L
    U_LL[EXT_DENS-1] = U_L[DENS-1];
    U_LL[EXT_SRHO-1] = 1.0/U_LL[EXT_DENS-1];
    U_LL[EXT_ENER-1] = U_L[ENER-1];
    U_LL[EXT_PRES-1] = UPrim_L[PRES-1];
    U_LL[EXT_TEMP-1] = UPrim_L[TEMP-1];

    // rotate velocity in normal and tangential direction
    U_LL[EXT_VEL1-1] = UPrim_L[VEL1-1]*nv[0] + UPrim_L[VEL2-1]*nv[1] + UPrim_L[VEL3-1]*nv[2];
    U_LL[EXT_VEL2-1] = UPrim_L[VEL1-1]*t1[0] + UPrim_L[VEL2-1]*t1[1] + UPrim_L[VEL3-1]*t1[2];
    U_LL[EXT_MOM1-1] = U_LL[EXT_DENS-1]*U_LL[EXT_VEL1-1];
    U_LL[EXT_MOM2-1] = U_LL[EXT_DENS-1]*U_LL[EXT_VEL2-1];
    U_LL[EXT_VEL3-1] = UPrim_L[VEL1-1]*t2[0] + UPrim_L[VEL2-1]*t2[1] + UPrim_L[VEL3-1]*t2[2];
    U_LL[EXT_MOM3-1] = U_LL[EXT_DENS-1]*U_LL[EXT_VEL3-1];

    // right state: U_R
    U_RR[EXT_DENS-1] = U_R[DENS-1];
    U_RR[EXT_SRHO-1] = 1.0/U_RR[EXT_DENS-1];
    U_RR[EXT_ENER-1] = U_R[ENER-1];
    U_RR[EXT_PRES-1] = UPrim_R[PRES-1];
    U_RR[EXT_TEMP-1] = UPrim_R[TEMP-1];

    // rotate momentum in normal and tangential direction
    U_RR[EXT_VEL1-1] = UPrim_R[VEL1-1]*nv[0] + UPrim_R[VEL2-1]*nv[1] + UPrim_R[VEL3-1]*nv[2];
    U_RR[EXT_VEL2-1] = UPrim_R[VEL1-1]*t1[0] + UPrim_R[VEL2-1]*t1[1] + UPrim_R[VEL3-1]*t1[2];
    U_RR[EXT_MOM1-1] = U_RR[EXT_DENS-1]*U_RR[EXT_VEL1-1];
    U_RR[EXT_MOM2-1] = U_RR[EXT_DENS-1]*U_RR[EXT_VEL2-1];
    U_RR[EXT_VEL3-1] = UPrim_R[VEL1-1]*t2[0] + UPrim_R[VEL2-1]*t2[1] + UPrim_R[VEL3-1]*t2[2];
    U_RR[EXT_MOM3-1] = U_RR[EXT_DENS-1]*U_RR[EXT_VEL3-1];

#ifndef SPLIT_DG
    EvalEulerFlux1D_Device(U_LL,F_L);
    EvalEulerFlux1D_Device(U_RR,F_R);
#endif /*SPLIT_DG*/

    if (doBC)
    {
        RiemannBC_Pointer_Device(F_L,F_R,U_LL,U_RR,riemannFlux);
    }
    else
    {
        Riemann_Pointer_Device(F_L,F_R,U_LL,U_RR,riemannFlux);
    }

    // Back rotate the normal flux into Cartesian direction
    Fout[DENS-1]=riemannFlux[DENS-1];
    for (int n = 0; n < 3; n++)
    {
        Fout[MOM1-1+n] = nv[n]*riemannFlux[MOM1-1] + t1[n]*riemannFlux[MOM2-1] + t2[n]*riemannFlux[MOM3-1];
    }
    Fout[ENER-1]=riemannFlux[ENER-1];
}


//----------------------------------------------------------------------------------------------------------------------------------
// Riemann solver support methods
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Device function rendering of the Fortran macro BETA_RIEMANN_H from eos.h
 **********************************************************************************************************************************/
__device__ double BetaRiemann_H(void)
{
    return sqrt(0.5*d_KappaM1 / d_Kappa);
}

/***********************************************************************************************************************************
 * @brief Device function rendering of the Fortran macro ROEC_RIEMANN_H from eos.h
 **********************************************************************************************************************************/
__device__ double RoecRiemann_H(double RoeH, double* RoeVel)
{
    return sqrt(d_KappaM1*(RoeH - 0.5*(RoeVel[0]*RoeVel[0] + RoeVel[1]*RoeVel[1] + RoeVel[2]*RoeVel[2])));
}

/***********************************************************************************************************************************
 * @brief Device function rendering of the Fortran macro ALPHA2_RIEMANN_H from eos.h
 **********************************************************************************************************************************/
__device__ double Alpha2Riemann_H(double RoeH, double* RoeVel, double Roec, double* Delta_U)
{
    double first_term = d_KappaM1/(Roec*Roec);
    double second_term = Delta_U[0]*(RoeH-RoeVel[0]*RoeVel[0]);
    double third_term = Delta_U[5] + RoeVel[0]*Delta_U[1];
    return first_term * (second_term - third_term);
}

/***********************************************************************************************************************************
 * @brief Type to hold all Roe fluxes, mean values and eigenvalues/vectors for Riemann solvers that need to compute them
 **********************************************************************************************************************************/
struct RoeData
{
    double H_L,H_R;
    double SqrtRho_L,SqrtRho_R,sSqrtRho;
    double RoeVel[3],RoeH,Roec,absVel;
    double a[PP_nVar],r1[PP_nVar],r2[PP_nVar],r3[PP_nVar],r4[PP_nVar],r5[PP_nVar];
};

/***********************************************************************************************************************************
 * @brief Method that hides computation of Roe fluxes, mean values and eigenvalues/vectors for all Riemann solvers that share it
 **********************************************************************************************************************************/
__device__ RoeData FindRoeData(double* U_LL, double* U_RR)
{
    struct RoeData roeData;

    roeData.H_L       = TotalEnthalpy_HE(U_LL);
    roeData.H_R       = TotalEnthalpy_HE(U_RR);
    roeData.SqrtRho_L = sqrt(U_LL[EXT_DENS-1]);
    roeData.SqrtRho_R = sqrt(U_RR[EXT_DENS-1]);

    roeData.sSqrtRho  = 1.0/(roeData.SqrtRho_L+roeData.SqrtRho_R);
    // Roe mean values
    for (int n = 0; n < 3; n++)
    {
        roeData.RoeVel[n] = roeData.sSqrtRho*(roeData.SqrtRho_R*U_RR[EXT_VEL1-1+n] + roeData.SqrtRho_L*U_LL[EXT_VEL1-1+n]);
    }
    roeData.absVel = roeData.RoeVel[0]*roeData.RoeVel[0] + roeData.RoeVel[1]*roeData.RoeVel[1] + roeData.RoeVel[2]*roeData.RoeVel[2];
    roeData.RoeH   = (roeData.SqrtRho_R*roeData.H_R+roeData.SqrtRho_L*roeData.H_L) * roeData.sSqrtRho;
    roeData.Roec   = RoecRiemann_H(roeData.RoeH,roeData.RoeVel);

    // mean eigenvalues and eigenvectors
    roeData.a[0] = roeData.RoeVel[0]-roeData.Roec;
    roeData.a[1] = roeData.RoeVel[0];
    roeData.a[2] = roeData.RoeVel[0];
    roeData.a[3] = roeData.RoeVel[0];
    roeData.a[4] = roeData.RoeVel[0]+roeData.Roec;

    roeData.r1[0] = 1.0;
    roeData.r1[1] = roeData.a[0];
    roeData.r1[2] = roeData.RoeVel[1];
    roeData.r1[3] = roeData.RoeVel[2];
    roeData.r1[4] = roeData.RoeH-roeData.RoeVel[0]*roeData.Roec;

    roeData.r2[0] = 1.0;
    roeData.r2[1] = roeData.RoeVel[0];
    roeData.r2[2] = roeData.RoeVel[1];
    roeData.r2[3] = roeData.RoeVel[2];
    roeData.r2[4] = 0.5*roeData.absVel;

    roeData.r3[0] = 0.0;
    roeData.r3[1] = 0.0;
    roeData.r3[2] = 1.0;
    roeData.r3[3] = 0.0;
    roeData.r3[4] = roeData.RoeVel[1];

    roeData.r4[0] = 0.0;
    roeData.r4[1] = 0.0;
    roeData.r4[2] = 0.0;
    roeData.r4[3] = 1.0;
    roeData.r4[4] = roeData.RoeVel[2];

    roeData.r5[0] = 1.0;
    roeData.r5[1] = roeData.a[4];
    roeData.r5[2] = roeData.RoeVel[1];
    roeData.r5[3] = roeData.RoeVel[2];
    roeData.r5[4] = roeData.RoeH+roeData.RoeVel[0]*roeData.Roec;

    return roeData;
}

//----------------------------------------------------------------------------------------------------------------------------------
// Riemann solver methods
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Local Lax-Friedrichs (Rusanov) Riemann solver
 **********************************************************************************************************************************/
__device__ void Riemann_LF_Device(double* F_L, double* F_R, double* U_LL, double* U_RR, double* riemannFlux)
{
    double lambdaMax = fmax(fabs(U_RR[EXT_VEL1-1]), fabs(U_LL[EXT_VEL1-1])) + fmax(SpeedOfSound_HE(U_LL), SpeedOfSound_HE(U_RR));

#ifndef SPLIT_DG
    for (int n = 0; n < PP_nVar; n++)
    {
        riemannFlux[n] = 0.5*((F_L[n] + F_R[n]) - lambdaMax*(U_RR[n] - U_LL[n]));
    }
#else
    for (int n = 0; n < PP_nVar; n++)
    {
        riemannFlux[n] = SplitDGSurface_Pointer_Device(n, U_LL, U_RR) - 0.5*lambdaMax*(U_RR[n] - U_LL[n]);
    }

#endif /* SPLIT_DG */
}

/***********************************************************************************************************************************
 * @brief Roe's approximate Riemann solver
 **********************************************************************************************************************************/
__device__ void Riemann_Roe_Device(double* F_L, double* F_R, double* U_LL, double* U_RR, double* riemannFlux)
{
    struct RoeData roeData;
    double Alpha[5];
    double Delta_U[PP_nVar+1];

    // Find Roe fluxes, mean values and eigenvalues/vectors
    roeData = FindRoeData(U_LL, U_RR);

    // calculate differences
    for (int n = 0; n < 5; n++)
    {
        Delta_U[n] = U_RR[n] - U_LL[n];
    }
    Delta_U[DELTA_U6-1] = Delta_U[DELTA_U5-1]-(Delta_U[DELTA_U3-1]-roeData.RoeVel[DELTA_U2-1]*Delta_U[DELTA_U1-1])*roeData.RoeVel[1] -
                        (Delta_U[DELTA_U4-1]-roeData.RoeVel[DELTA_U3-1]*Delta_U[DELTA_U1-1])*roeData.RoeVel[DELTA_U3-1];
    // calculate factors
    Alpha[2] = Delta_U[DELTA_U3-1] - roeData.RoeVel[DELTA_U2-1]*Delta_U[DELTA_U1-1];
    Alpha[3] = Delta_U[DELTA_U4-1] - roeData.RoeVel[DELTA_U3-1]*Delta_U[DELTA_U1-1];
    Alpha[1] = Alpha2Riemann_H(roeData.RoeH,roeData.RoeVel,roeData.Roec,Delta_U);
    Alpha[0] = 0.5/roeData.Roec * (Delta_U[DELTA_U1-1]*(roeData.RoeVel[0]+roeData.Roec) - Delta_U[DELTA_U2-1] - roeData.Roec*Alpha[1]);
    Alpha[4] = Delta_U[DELTA_U1-1] - Alpha[0] - Alpha[1];
#ifndef SPLIT_DG
    // assemble Roe flux
    for (int n = 0; n < PP_nVar; n++)
    {
        riemannFlux[n] = 0.5*((F_L[n] + F_R[n]) - 
                    Alpha[0]*fabs(roeData.a[0])*roeData.r1[n] -
                    Alpha[1]*fabs(roeData.a[1])*roeData.r2[n] -
                    Alpha[2]*fabs(roeData.a[2])*roeData.r3[n] -
                    Alpha[3]*fabs(roeData.a[3])*roeData.r4[n] -
                    Alpha[4]*fabs(roeData.a[4])*roeData.r5[n]);
    }
#else
    for (int n = 0; n < PP_nVar; n++)
    {
        riemannFlux[n] = SplitDGSurface_Pointer_Device(n, U_LL,U_RR); // get split flux
        // assemble Roe flux
        riemannFlux[n] -= 0.5*(Alpha[0]*fabs(roeData.a[0])*roeData.r1[n] +
                    Alpha[1]*fabs(roeData.a[1])*roeData.r2[n] +
                    Alpha[2]*fabs(roeData.a[2])*roeData.r3[n] +
                    Alpha[3]*fabs(roeData.a[3])*roeData.r4[n] +
                    Alpha[4]*fabs(roeData.a[4])*roeData.r5[n]);
    }
#endif /*SPLIT_DG*/ 
}

/***********************************************************************************************************************************
 * @brief Roe's approximate Riemann solver using the Harten and Hymen II entropy fix, see
 *        Pelanti, Marica & Quartapelle, Luigi & Vigevano, L & Vigevano, Luigi. (2018):
 *        A review of entropy fixes as applied to Roe's linearization.
 **********************************************************************************************************************************/
__device__ void Riemann_RoeEntropyFix_Device(double* F_L, double* F_R, double* U_LL, double* U_RR, double* riemannFlux)
{
    double c_L,c_R;
    double al[PP_nVar],ar[PP_nVar];
    double Alpha[PP_nVar];
    double Delta_U[PP_nVar];
    double tmp, da;
    double RoeDens;
    struct RoeData roeData;

    roeData = FindRoeData(U_LL, U_RR);
    c_L = SpeedOfSound_HE(U_LL);
    c_R = SpeedOfSound_HE(U_RR);
    RoeDens= sqrt(U_LL[EXT_DENS-1]*U_RR[EXT_DENS-1]);

    // Calculate jump
    Delta_U[DELTA_U1-1] = U_RR[EXT_DENS-1] - U_LL[EXT_DENS-1];
    Delta_U[DELTA_U2-1] = U_RR[EXT_VEL1-1] - U_LL[EXT_VEL1-1];
    Delta_U[DELTA_U3-1] = U_RR[EXT_VEL2-1] - U_LL[EXT_VEL2-1];
    Delta_U[DELTA_U4-1] = U_RR[EXT_VEL3-1] - U_LL[EXT_VEL3-1];
    Delta_U[DELTA_U5-1] = U_RR[EXT_PRES-1] - U_LL[EXT_PRES-1];

    // calculate wave strenghts
    tmp      = 0.5/(roeData.Roec*roeData.Roec);
    Alpha[0] = tmp*(Delta_U[DELTA_U5-1]-RoeDens*roeData.Roec*Delta_U[DELTA_U2-1]);
    Alpha[1] = Delta_U[DELTA_U1-1] - Delta_U[DELTA_U5-1]*2.0*tmp;
    Alpha[2] = RoeDens*Delta_U[DELTA_U3-1];
    Alpha[3] = RoeDens*Delta_U[DELTA_U4-1];
    Alpha[4] = tmp*(Delta_U[DELTA_U5-1]+RoeDens*roeData.Roec*Delta_U[DELTA_U2-1]);

    // Harten+Hyman entropy fix (apply only for acoustic waves, don't fix r)
    al[0] = U_LL[EXT_VEL1-1] - c_L;
    al[1] = U_LL[EXT_VEL1-1];
    al[2] = U_LL[EXT_VEL1-1];
    al[3] = U_LL[EXT_VEL1-1];
    al[4] = U_LL[EXT_VEL1-1] + c_L;
    ar[0] = U_RR[EXT_VEL1-1] - c_R;
    ar[1] = U_RR[EXT_VEL1-1];
    ar[2] = U_RR[EXT_VEL1-1];
    ar[3] = U_RR[EXT_VEL1-1];
    ar[4] = U_RR[EXT_VEL1-1] + c_R;

    for (int n = 0; n < 5; n++)
    {
        da = fmax(roeData.a[n]-al[n], ar[n]-roeData.a[n]);
        da = fmax(da, 0.0);
        if (fabs(roeData.a[n]) < da)
        {
            roeData.a[n] = 0.5*(roeData.a[n]*roeData.a[n]/da+da);
        }
        else
        {
            roeData.a[n] = fabs(roeData.a[n]);
        }
    }

#ifndef SPLIT_DG
    // assemble Roe flux
    for (int n = 0; n < PP_nVar; n++)
    {
        riemannFlux[n] = 0.5*((F_L[n] + F_R[n]) - 
                    Alpha[0]*roeData.a[0]*roeData.r1[n] -
                    Alpha[1]*roeData.a[1]*roeData.r2[n] -
                    Alpha[2]*roeData.a[2]*roeData.r3[n] -
                    Alpha[3]*roeData.a[3]*roeData.r4[n] -
                    Alpha[4]*roeData.a[4]*roeData.r5[n]);
    }
#else
    for (int n = 0; n < PP_nVar; n++)
    {
        riemannFlux[n] = SplitDGSurface_Pointer_Device(n,U_LL,U_RR); // get split flux
        // assemble Roe flux
        riemannFlux[n] -= 0.5*(Alpha[0]*roeData.a[0]*roeData.r1[n] +
                    Alpha[1]*roeData.a[1]*roeData.r2[n] +
                    Alpha[2]*roeData.a[2]*roeData.r3[n] +
                    Alpha[3]*roeData.a[3]*roeData.r4[n] +
                    Alpha[4]*roeData.a[4]*roeData.r5[n]);
    }
#endif /*SPLIT_DG*/ 
}

/***********************************************************************************************************************************
 * @brief Low mach number Roe's approximate Riemann solver according to Oßwald(2015)
 **********************************************************************************************************************************/
__device__ void Riemann_RoeL2_Device(double* F_L, double* F_R, double* U_LL, double* U_RR, double* riemannFlux)
{
    double Ma_loc;
    double Alpha[5];
    double Delta_U[PP_nVar+1];
    struct RoeData roeData;

    roeData = FindRoeData(U_LL, U_RR);

    // calculate differences
    Ma_loc = sqrt(roeData.absVel) / (roeData.Roec*sqrt(d_Kappa));
    for (int n = 0; n < 5; n++)
    {
        Delta_U[n] = U_RR[n] - U_LL[n];
        // Low Mach number fix
        if (n > 0 && n < 4)
        {
            Delta_U[n] *= Ma_loc;
        }
    }
    Delta_U[DELTA_U6-1] = Delta_U[DELTA_U5-1]-(Delta_U[DELTA_U3-1]-roeData.RoeVel[1]*Delta_U[DELTA_U1-1])*roeData.RoeVel[1] -
                        (Delta_U[DELTA_U4-1]-roeData.RoeVel[2]*Delta_U[DELTA_U1-1])*roeData.RoeVel[2];
    
    // calculate factors
    Alpha[2] = Delta_U[DELTA_U3-1] - roeData.RoeVel[DELTA_U2-1]*Delta_U[DELTA_U1-1];
    Alpha[3] = Delta_U[DELTA_U4-1] - roeData.RoeVel[DELTA_U3-1]*Delta_U[DELTA_U1-1];
    Alpha[1] = Alpha2Riemann_H(roeData.RoeH,roeData.RoeVel,roeData.Roec,Delta_U);
    Alpha[0] = 0.5/roeData.Roec * (Delta_U[DELTA_U1-1]*(roeData.RoeVel[0]+roeData.Roec) - Delta_U[DELTA_U2-1] - roeData.Roec*Alpha[1]);
    Alpha[4] = Delta_U[DELTA_U1-1] - Alpha[0] - Alpha[1];

#ifndef SPLIT_DG
    // assemble Roe flux
    for (int n = 0; n < PP_nVar; n++)
    {
        riemannFlux[n] = 0.5*((F_L[n] + F_R[n]) - 
                    Alpha[0]*fabs(roeData.a[0])*roeData.r1[n] -
                    Alpha[1]*fabs(roeData.a[1])*roeData.r2[n] -
                    Alpha[2]*fabs(roeData.a[2])*roeData.r3[n] -
                    Alpha[3]*fabs(roeData.a[3])*roeData.r4[n] -
                    Alpha[4]*fabs(roeData.a[4])*roeData.r5[n]);
    }
#else
    for (int n = 0; n < PP_nVar; n++)
    {
        riemannFlux[n] = SplitDGSurface_Pointer_Device(n,U_LL,U_RR); // get split flux
        // assemble Roe flux
        riemannFlux[n] -= 0.5*(Alpha[0]*fabs(roeData.a[0])*roeData.r1[n] +
                    Alpha[1]*fabs(roeData.a[1])*roeData.r2[n] +
                    Alpha[2]*fabs(roeData.a[2])*roeData.r3[n] +
                    Alpha[3]*fabs(roeData.a[3])*roeData.r4[n] +
                    Alpha[4]*fabs(roeData.a[4])*roeData.r5[n]);
    }
#endif /*SPLIT_DG*/ 
}

#ifndef SPLIT_DG
/***********************************************************************************************************************************
 * @brief Standard Harten-Lax-Van-Leer Riemann solver without contact discontinuity
 **********************************************************************************************************************************/
__device__ void Riemann_HLL_Device(double* F_L, double* F_R, double* U_LL, double* U_RR, double* riemannFlux)
{
    struct RoeData roeData;
    double ssl, ssr;

    roeData = FindRoeData(U_LL, U_RR);

    // HLL flux
    ssl = roeData.RoeVel[0] - roeData.Roec;
    ssr = roeData.RoeVel[1] + roeData.Roec;

    for (int n = 0; n < PP_nVar; n++)
    {
        if (ssl >= 0.0)
        {
            riemannFlux[n] = F_L[n];
        }
        else if (ssr <= 0.0)
        {
            riemannFlux[n] = F_R[n];
        }
        else
        {
            riemannFlux[n] = (ssr*F_L[n]-ssl*F_R[n]+ssl*ssr*(U_RR[n]-U_LL[n]))/(ssr-ssl);
        }
    }
}

/***********************************************************************************************************************************
 * @brief Harten-Lax-Van-Leer Riemann solver resolving contact discontinuity
 **********************************************************************************************************************************/
__device__ void Riemann_HLLC_Device(double* F_L, double* F_R, double* U_LL, double* U_RR, double* riemannFlux)
{
    struct RoeData roeData;
    double ssl, ssr;
    double sMu_L, sMu_R;
    double SStar; // Coolest variable name I've ever seen, personally. (Spencer Starr, 14.08.24)
    double U_Star, U_Star_factor[PP_nVar];

    roeData = FindRoeData(U_LL, U_RR);

    // HLL flux
    ssl = roeData.RoeVel[0] - roeData.Roec;
    ssr = roeData.RoeVel[1] + roeData.Roec;

    sMu_L = ssl - U_LL[EXT_VEL1-1];
    sMu_R = ssr - U_RR[EXT_VEL1-1];
    SStar = (U_RR[EXT_PRES-1] - U_LL[EXT_PRES-1] + U_LL[EXT_MOM1-1]*sMu_L - U_RR[EXT_MOM1-1]*sMu_R) / 
                        (U_LL[EXT_DENS-1]*sMu_L - U_RR[EXT_DENS-1]*sMu_R);
    if (ssl <= 0.0 && SStar >= 0.0)
    {
        U_Star_factor[0] = 1.0;
        U_Star_factor[1] = SStar;
        U_Star_factor[2] = U_LL[EXT_VEL2-1];
        U_Star_factor[3] = U_LL[EXT_VEL3-1];
        U_Star_factor[4] = TotalEnergy_HE(U_LL) + 
                        (SStar-U_LL[EXT_VEL1-1])*(SStar + U_LL[EXT_PRES-1]*U_LL[EXT_SRHO-1]/sMu_L);
    }
    else
    {
        U_Star_factor[0] = 1.0;
        U_Star_factor[1] = SStar;
        U_Star_factor[2] = U_RR[EXT_VEL2-1];
        U_Star_factor[3] = U_RR[EXT_VEL3-1];
        U_Star_factor[4] = TotalEnergy_HE(U_RR) + 
                        (SStar-U_LL[EXT_VEL1-1])*(SStar + U_RR[EXT_PRES-1]*U_RR[EXT_SRHO-1]/sMu_R);
    }

    for (int n = 0; n < PP_nVar; n++)
    {
        if (ssl >= 0.0)
        {
            riemannFlux[n] = F_L[n];
        }
        else if (ssr <= 0.0)
        {
            riemannFlux[n] = F_R[n];
        }
        else
        {
            if (ssl <= 0.0 && SStar >= 0.0)
            {
                U_Star = U_LL[EXT_DENS-1] * sMu_L/(ssl-SStar) * U_Star_factor[n];
                riemannFlux[n] = F_L[n] + ssl*(U_Star-U_LL[n]);
            }
            else
            {
                U_Star = U_RR[EXT_DENS-1] * sMu_R/(ssr-SStar) * U_Star_factor[n];
                riemannFlux[n] = F_R[n] + ssr*(U_Star-U_RR[n]);
            }
        }
    }
}

/***********************************************************************************************************************************
 * @brief Harten-Lax-Van-Leer-Einfeldt Riemann solver
 **********************************************************************************************************************************/
__device__ void Riemann_HLLE_Device(double* F_L, double* F_R, double* U_LL, double* U_RR, double* riemannFlux)
{
    struct RoeData roeData;
    double ssl, ssr;
    double beta;

    roeData = FindRoeData(U_LL, U_RR);

    // HLL flux
    beta = BetaRiemann_H();
    ssl = fmin(roeData.RoeVel[0]-roeData.Roec, U_LL[EXT_VEL1-1]-beta*SpeedOfSound_HE(U_LL));
    ssl = fmin(ssl, 0.0);
    ssr = fmax(roeData.RoeVel[0]-roeData.Roec, U_RR[EXT_VEL1-1]+beta*SpeedOfSound_HE(U_RR));
    ssr = fmax(ssl, 0.0);

    for (int n = 0; n < PP_nVar; n++)
    {
        if (ssl >= 0.0)
        {
            riemannFlux[n] = F_L[n];
        }
        else if (ssr <= 0.0)
        {
            riemannFlux[n] = F_R[n];
        }
        else
        {
            riemannFlux[n] = (ssr*F_L[n]-ssl*F_R[n]+ssl*ssr*(U_RR[n]-U_LL[n]))/(ssr-ssl);
        }
    }
}

/***********************************************************************************************************************************
 * @brief Harten-Lax-Van-Leer-Einfeldt-Munz Riemann solver
 **********************************************************************************************************************************/
__device__ void Riemann_HLLEM_Device(double* F_L, double* F_R, double* U_LL, double* U_RR, double* riemannFlux)
{
    struct RoeData roeData;
    double RoeDens;
    double ssl, ssr;
    double beta, delta, corr_term;
    double Alpha[3];

    roeData = FindRoeData(U_LL, U_RR);
    RoeDens= sqrt(U_LL[EXT_DENS-1]*U_RR[EXT_DENS-1]);

    // HLL flux
    beta = BetaRiemann_H();
    ssl = fmin(roeData.RoeVel[0]-roeData.Roec, U_LL[EXT_VEL1-1]-beta*SpeedOfSound_HE(U_LL));
    ssl = fmin(ssl, 0.0);
    ssr = fmax(roeData.RoeVel[0]-roeData.Roec, U_RR[EXT_VEL1-1]+beta*SpeedOfSound_HE(U_RR));
    ssr = fmax(ssl, 0.0);

    // Mean eigenvectors
    Alpha[0]   = (U_RR[EXT_DENS-1]-U_LL[EXT_DENS-1])  - (U_RR[EXT_PRES-1]-U_LL[EXT_PRES-1])/(roeData.Roec*roeData.Roec);
    Alpha[1] = RoeDens*(U_RR[EXT_VEL2-1] - U_LL[EXT_VEL2-1]);
    Alpha[1] = RoeDens*(U_RR[EXT_VEL3-1] - U_LL[EXT_VEL3-1]);

    for (int n = 0; n < PP_nVar; n++)
    {
        if (ssl >= 0.0)
        {
            riemannFlux[n] = F_L[n];
        }
        else if (ssr <= 0.0)
        {
            riemannFlux[n] = F_R[n];
        }
        else
        {
            corr_term = delta*(roeData.r2[n]*Alpha[0]+roeData.r3[n]*Alpha[1]+roeData.r4[n]*Alpha[2]);
            riemannFlux[n] = (ssr*F_L[n]-ssl*F_R[n]+ssl*ssr*(U_RR[n]-U_LL[n] - corr_term))/(ssr-ssl);
        }
    }
}

#else /* SPLIT_DG */
/***********************************************************************************************************************************
 * @brief Riemann solver using purely the average fluxes
 **********************************************************************************************************************************/
__device__ void Riemann_FluxAverage_Device(double* F_L, double* F_R, double* U_LL, double* U_RR, double* riemannFlux)
{
    for(int n = 0; n < PP_nVar; n++)
    {
        riemannFlux[n] = SplitDGSurface_Pointer_Device(n,U_LL, U_RR);
    }
}

/***********************************************************************************************************************************
 * @brief Kinetic energy preserving and entropy consistent flux according to Chandrashekar (2012)
 **********************************************************************************************************************************/
__device__ void Riemann_CH_Device(double* F_L, double* F_R, double* U_LL, double* U_RR, double* riemannFlux)
{
    double lambdaMax = fmax(fabs(U_RR[EXT_VEL1-1]), fabs(U_LL[EXT_VEL1-1])) + fmax(SpeedOfSound_HE(U_LL), SpeedOfSound_HE(U_RR));

    double betaLogMean = 0.0;

    // average velocities
    double rhoMean = 0.5*(U_LL[EXT_DENS-1] + U_RR[EXT_DENS-1]);
    double uMean = 0.5*(U_LL[EXT_VEL1-1] + U_RR[EXT_VEL1-1]);
    double vMean = 0.5*(U_LL[EXT_VEL2-1] + U_RR[EXT_VEL2-1]);
    double wMean = 0.5*(U_LL[EXT_VEL3-1] + U_RR[EXT_VEL3-1]);

    // inverse temperature
    double beta_LL = 0.5*U_LL[EXT_DENS-1]/U_LL[EXT_PRES-1];
    double beta_RR = 0.5*U_RR[EXT_DENS-1]/U_RR[EXT_PRES-1];

    // average pressure, enthalpy, density and inverse temperature
    // logarithmic mean
    betaLogMean = GetLogMean_Device(&beta_LL,&beta_RR);

    // Compute the flux
    for(int n = 0; n < PP_nVar; n++)
    {
        riemannFlux[n] = SplitDGSurface_Pointer_Device(n,U_LL, U_RR);
        if (n < MOM3)
        {
            riemannFlux[n] -= 0.5*lambdaMax*(U_RR[n] - U_LL[n]);
        }
    }
    riemannFlux[ENER-1] -= 0.5*lambdaMax*(
                                          (U_RR[EXT_DENS-1]-U_LL[EXT_DENS-1])*(0.5*d_sKappaM1/betaLogMean + 
                                           0.5*(U_RR[EXT_VEL1-1]*U_LL[EXT_VEL1-1]+U_RR[EXT_VEL2-1]*U_LL[EXT_VEL2-1]+U_RR[EXT_VEL3-1]*U_LL[EXT_VEL3-1])) +
                                           rhoMean*uMean*(U_RR[EXT_VEL1-1]-U_LL[EXT_VEL1-1]) + 
                                           rhoMean*vMean*(U_RR[EXT_VEL2-1]-U_LL[EXT_VEL2-1]) + 
                                           rhoMean*wMean*(U_RR[EXT_VEL3-1]-U_LL[EXT_VEL3-1]) +
                                           0.5*rhoMean*d_sKappaM1*(1./beta_RR - 1./beta_LL)
                                         );
}
#endif /* SPLIT_DG */

//----------------------------------------------------------------------------------------------------------------------------------
// Viscous flux methods
//----------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------
// Initialization functions
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Initialize function pointers for Riemann solver methods on the device side. This method runs on the HOST only
 * @param riemannOpt User chosen option for the Riemann solver from the input file
 * @param riemannOptBC User chosen option for the BC Riemann solver from the input file
 * @remark We trust that illegal value checking is done by the CPU before this kernel is ever called, so we don't check here.
 **********************************************************************************************************************************/
__global__ void InitRiemann_Kernel(int riemannOpt, int riemannOptBC)
{
    // Set pointer for general Riemann solver
    switch (riemannOpt)
    {
        case PRM_RIEMANN_LF:
            Riemann_Pointer_Device = Riemann_LF_Device;
            break;
        case PRM_RIEMANN_ROE:
            Riemann_Pointer_Device = Riemann_Roe_Device;
            break;
        case PRM_RIEMANN_ROEENTROPYFIX:
            Riemann_Pointer_Device = Riemann_RoeEntropyFix_Device;
            break;
        case PRM_RIEMANN_ROEL2:
            Riemann_Pointer_Device = Riemann_RoeL2_Device;
            break;
#ifndef SPLIT_DG
        case PRM_RIEMANN_HLL:
            Riemann_Pointer_Device = Riemann_HLL_Device;
            break;
        case PRM_RIEMANN_HLLC:
            Riemann_Pointer_Device = Riemann_HLLC_Device;
            break;
        case PRM_RIEMANN_HLLE:
            Riemann_Pointer_Device = Riemann_HLLE_Device;
            break;
        case PRM_RIEMANN_HLLEM:
            Riemann_Pointer_Device = Riemann_HLLEM_Device;
            break;
#else /* SPLIT_DG */
        case PRM_RIEMANN_CH:
            Riemann_Pointer_Device = Riemann_CH_Device;
            break;
        case PRM_RIEMANN_Average:
            Riemann_Pointer_Device = Riemann_FluxAverage_Device;
            break;
#endif /* SPLIT_DG */
    }

    // Set pointer for BC Riemann solver
    switch (riemannOptBC)
    {
        case PRM_RIEMANN_SAME:
            RiemannBC_Pointer_Device = Riemann_Pointer_Device;
            break;
        case PRM_RIEMANN_LF:
            RiemannBC_Pointer_Device = Riemann_LF_Device;
            break;
        case PRM_RIEMANN_ROE:
            RiemannBC_Pointer_Device = Riemann_Roe_Device;
            break;
        case PRM_RIEMANN_ROEENTROPYFIX:
            RiemannBC_Pointer_Device = Riemann_RoeEntropyFix_Device;
            break;
        case PRM_RIEMANN_ROEL2:
            RiemannBC_Pointer_Device = Riemann_RoeL2_Device;
            break;
#ifndef SPLIT_DG
        case PRM_RIEMANN_HLL:
            RiemannBC_Pointer_Device = Riemann_HLL_Device;
            break;
        case PRM_RIEMANN_HLLC:
            RiemannBC_Pointer_Device = Riemann_HLLC_Device;
            break;
        case PRM_RIEMANN_HLLE:
            RiemannBC_Pointer_Device = Riemann_HLLE_Device;
            break;
        case PRM_RIEMANN_HLLEM:
            RiemannBC_Pointer_Device = Riemann_HLLEM_Device;
            break;
#else /* SPLIT_DG */
        case PRM_RIEMANN_CH:
            RiemannBC_Pointer_Device = Riemann_CH_Device;
            break;
        case PRM_RIEMANN_Average:
            RiemannBC_Pointer_Device = Riemann_FluxAverage_Device;
            break;
#endif  /* SPLIT_DG */
    }

}

/***********************************************************************************************************************************
 * @brief Entry point for the device side initialization of the Riemann solver function pointers
 **********************************************************************************************************************************/
void InitRiemann_Device(int riemannOpt, int riemannOptBC)
{
  INVOKE_KERNEL( InitRiemann_Kernel, 1, 1, 0, streams[0], riemannOpt, riemannOptBC );
}