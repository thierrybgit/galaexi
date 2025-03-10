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
    void EvalTransformedFlux3D_Device(int Nloc,int nElems,int d_U,int d_UPrim
#if PARABOLIC
                                          ,int d_gradUx,int d_gradUy,int d_gradUz
#endif
                                          ,int d_F,int d_G,int d_H
                                          ,int d_Metrics_fTilde,int d_Metrics_gTilde,int d_Metrics_hTilde, int myStream);
}

//----------------------------------------------------------------------------------------------------------------------------------
// EvalTransformedFlux3D METHODS
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Device backend core of the computation advection part of the Navier-Stokes fluxes in all space dimensions using the 
 *        conservative and primitive variables
 **********************************************************************************************************************************/
__device__ void EvalTransformedEulerFlux3D_fast_Device(double* U, double* UPrim, double* F, double* G, double* H, double* Mf, double* Mg, double* Mh)
{
    // LOCAL VARIABLES
    double Ep;
    double Mmom;

    // auxiliary variables
    Ep = (U[ENER-1] + UPrim[PRES-1])/U[DENS-1];

    // fluxes in x-direction
    Mmom    = Mf[0]*U[MOM1-1] + Mf[1]*U[MOM2-1] + Mf[2]*U[MOM3-1];
    F[DENS-1] =  Mmom;
    F[MOM1-1] =  Mmom*UPrim[VEL1-1] + Mf[0]*UPrim[PRES-1];
    F[MOM2-1] =  Mmom*UPrim[VEL2-1] + Mf[1]*UPrim[PRES-1];
    F[MOM3-1] =  Mmom*UPrim[VEL3-1] + Mf[2]*UPrim[PRES-1];
    F[ENER-1] =  Mmom*Ep;

    // fluxes in y-direction
    Mmom    = Mg[0]*U[MOM1-1] + Mg[1]*U[MOM2-1] + Mg[2]*U[MOM3-1];
    G[DENS-1] =  Mmom;
    G[MOM1-1] =  Mmom*UPrim[VEL1-1] + Mg[0]*UPrim[PRES-1];
    G[MOM2-1] =  Mmom*UPrim[VEL2-1] + Mg[1]*UPrim[PRES-1];
    G[MOM3-1] =  Mmom*UPrim[VEL3-1] + Mg[2]*UPrim[PRES-1];
    G[ENER-1] =  Mmom*Ep;

    // fluxes in z-direction
    Mmom    = Mh[0]*U[MOM1-1] + Mh[1]*U[MOM2-1] + Mh[2]*U[MOM3-1];
    H[DENS-1] =  Mmom;
    H[MOM1-1] =  Mmom*UPrim[VEL1-1] + Mh[0]*UPrim[PRES-1];
    H[MOM2-1] =  Mmom*UPrim[VEL2-1] + Mh[1]*UPrim[PRES-1];
    H[MOM3-1] =  Mmom*UPrim[VEL3-1] + Mh[2]*UPrim[PRES-1];
    H[ENER-1] =  Mmom*Ep;
}

/***********************************************************************************************************************************
 * @brief Computes the advection part of the Navier-Stokes fluxes on the device
 **********************************************************************************************************************************/
__global__ void EvalTransFormedFlux3D_Kernel(int Nloc, int nElems, double* U, double* UPrim
#if PARABOLIC
                                                                             , double* gradUx, double* gradUy, double* gradUz
#endif
                                                                             , double* F , double* G , double* H
                                                                             , double* Mf, double* Mg, double* Mh)
{
    int U_Offset = 0;
    int UPrim_Offset = 0;
    int Metrics_Offset = 0;
    int threadID = blockIdx.x*blockDim.x + threadIdx.x;
    int nDOF = (Nloc+1)*(Nloc+1)*(Nloc+1)*nElems;
    if (threadID < nDOF)
    {
        // Unlike on the CPU side, we have the full arrays here. We need to index to the start of the desired element block
        // The offset to the DOF the current thread is working on is then the threadID*nVars in the array being index.
        U_Offset = threadID*PP_nVar;
        UPrim_Offset = threadID*PP_nVarPrim;
        Metrics_Offset = threadID*3;
#if PARABOLIC
        // WILL ADD IN THE KERNEL CALL HERE FOR EvalTransformedEulerDiffFlux3D AFTER IT IS WRITTEN
#else
        EvalTransformedEulerFlux3D_fast_Device(&U[U_Offset], &UPrim[UPrim_Offset],
                                          &F[U_Offset], &G[U_Offset], &H[U_Offset],
                                          &Mf[Metrics_Offset], &Mg[Metrics_Offset], &Mh[Metrics_Offset]
                                         );
#endif
    }
}

/***********************************************************************************************************************************
 * @brief Wrapper method for the device-side kernel for the transformed flux computation
 **********************************************************************************************************************************/
void EvalTransformedFlux3D_Device(int Nloc,int nElems,int d_U,int d_UPrim
#if PARABOLIC
                                          ,int d_gradUx,int d_gradUy,int d_gradUz
#endif
                                          ,int d_F,int d_G,int d_H
                                          ,int d_Metrics_fTilde,int d_Metrics_gTilde,int d_Metrics_hTilde, int myStream)
{

#ifdef USE_NVTX /* If profiling, add range */
    nvtxRangePushA("EvalTransFormedFlux3D");
#endif

    int nDOF = (Nloc+1)*(Nloc+1)*(Nloc+1)*nElems;
    INVOKE_KERNEL( EvalTransFormedFlux3D_Kernel, nDOF/256+1, 256, 0, streams[myStream], Nloc, nElems, (double*)DeviceVars[d_U], (double*)DeviceVars[d_UPrim]
#if PARABOLIC
                                                                             , (double*)DeviceVars[d_gradUx], (double*)DeviceVars[d_gradUy], (double*)DeviceVars[d_gradUz]
#endif
                                                                             , (double*)DeviceVars[d_F] , (double*)DeviceVars[d_G] , (double*)DeviceVars[d_H]
                                                                             , (double*)DeviceVars[d_Metrics_fTilde], (double*)DeviceVars[d_Metrics_gTilde], (double*)DeviceVars[d_Metrics_hTilde] );

#ifdef USE_NVTX /* If profiling, end range */
    nvtxRangePop(); 
#endif

}

//----------------------------------------------------------------------------------------------------------------------------------
// EvalDiffFlux3D METHODS
//----------------------------------------------------------------------------------------------------------------------------------
#if PARABOLIC
#endif

//----------------------------------------------------------------------------------------------------------------------------------
// EvalEulerFlux1D_fast METHOD
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Device backend for computing the 1D Euler flux. CPU backend found in flux.f90
 * @remark Method only accessible from device code as it is only ever called from within another kernel.
 * @param U Vector of cons and prim variables at a single DOF
 * @param F Cartesian flux in x-direction at that same single DOF
 **********************************************************************************************************************************/
__device__ void EvalEulerFlux1D_Device(double* U, double* F)
{
    F[DENS-1] = U[EXT_MOM1-1];
    F[MOM1-1] = U[EXT_MOM1-1]*U[EXT_VEL1-1]+U[EXT_PRES-1];  // rho*u²+p
    F[MOM2-1] = U[EXT_MOM1-1]*U[EXT_VEL2-1];                // rho*u*v
    F[MOM3-1] = U[EXT_MOM1-1]*U[EXT_VEL3-1];                // rho*u*w
    F[ENER-1] =(U[EXT_ENER-1]+U[EXT_PRES-1])*U[EXT_VEL1-1]; // (rho*e+p)*u
}