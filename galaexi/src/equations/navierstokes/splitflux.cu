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
    void InitSplitDG_Device(int splitDGOpt);

}


//----------------------------------------------------------------------------------------------------------------------------------
// Split DG formulation function pointers
//----------------------------------------------------------------------------------------------------------------------------------
/// Function pointer to Split DG Volume method
typedef double (*SplitDGVolumePointer_t)(int n, double* Uref, double* UPrimRef, double* U, double* UPrim, double* Mref, double* M);
__device__ SplitDGVolumePointer_t SplitDGVolume_Pointer_Device;

/// Function pointer to Split DG Surface method
typedef double (*SplitDGSurfacePointer_t)(int n, double* U_LL, double* U_RR);
__device__ SplitDGSurfacePointer_t SplitDGSurface_Pointer_Device;


//----------------------------------------------------------------------------------------------------------------------------------
// Other variables local to this file
//----------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------
// Utility methods for Split Flux methods
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Util function for calculating the logarithmic mean numerically stable according to Ismail and Roe
 **********************************************************************************************************************************/
__device__ double GetLogMean_Device(double* U_L, double* U_R)
{
  double UMean = 0.0;
  double N = 0.0;
  double epsilon = 0.01;

  double chi = *U_L / *U_R;
  double f = (chi-1)/(chi+1);
  double u = f*f;

  if (u < epsilon)
  {
    N = 1.0+u/3.0+u*u/5.0+u*u*u/7.0;
  }
  else
  {
    N = log(chi)/(2.0*f);
  }

  UMean = (*U_L+*U_R)/(2.0*N);
  return UMean;
}


//----------------------------------------------------------------------------------------------------------------------------------
// Split DG Volume Formulation Device Methods
//
// NOTE: All of these methods operate on a single variable at once. This is because it is easier to do than
//       return an array from a C function. It also works in a more clean way with the logic of the methods that
//       call these methods.
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Computes the Split-Flux retaining the standard NS-Equations
 **********************************************************************************************************************************/
__device__ double SplitVolumeFluxSD_Device(int n, double* Uref, double* UPrimRef, double* U, double* UPrim, double* Mref, double* M)
{
  double fTilde, gTilde, hTilde;
  double metricCoeff1 = 0.5*(Mref[0]+M[0]);
  double metricCoeff2 = 0.5*(Mref[1]+M[1]);
  double metricCoeff3 = 0.5*(Mref[2]+M[2]);
  double rhoEpRef = Uref[ENER-1] + UPrimRef[PRES-1];
  double rhoEp = U[ENER-1] + UPrim[PRES-1];

  if (n == DENS-1)
  {
    fTilde = Uref[MOM1-1] + U[MOM1-1];
    gTilde = Uref[MOM2-1] + U[MOM2-1];
    hTilde = Uref[MOM3-1] + U[MOM3-1];
  }
  else if (n == MOM1-1)
  {
    fTilde = Uref[MOM1-1]*UPrimRef[VEL1-1]+UPrimRef[PRES-1] + U[MOM1-1]*UPrim[VEL1-1]+UPrim[PRES-1];
    gTilde = Uref[MOM1-1]*UPrimRef[VEL2-1] + U[MOM1-1]*UPrim[VEL2-1];
    hTilde = Uref[MOM1-1]*UPrimRef[VEL3-1] + U[MOM1-1]*UPrim[VEL3-1];
  }
  else if (n == MOM2-1)
  {
    fTilde = Uref[MOM1-1]*UPrimRef[VEL2-1] + U[MOM1-1]*UPrim[VEL2-1];
    gTilde = Uref[MOM2-1]*UPrimRef[VEL2-1]+UPrimRef[PRES-1] + U[MOM2-1]*UPrim[VEL2-1]+UPrim[PRES-1];
    hTilde = Uref[MOM2-1]*UPrimRef[VEL3-1] + U[MOM2-1]*UPrim[VEL3-1];
  }
  else if (n == MOM3-1)
  {
    fTilde = Uref[MOM1-1]*UPrimRef[VEL3-1] + U[MOM1-1]*UPrim[VEL3-1];
    gTilde = Uref[MOM2-1]*UPrimRef[VEL3-1] + U[MOM2-1]*UPrim[VEL3-1];
    hTilde = Uref[MOM3-1]*UPrimRef[VEL3-1]+UPrimRef[PRES-1] + U[MOM3-1]*UPrim[VEL3-1]+UPrim[PRES-1];
  }
  else if (n == ENER-1)
  {
    fTilde = rhoEpRef*UPrimRef[VEL1-1] + rhoEp*UPrim[VEL1-1];
    gTilde = rhoEpRef*UPrimRef[VEL2-1] + rhoEp*UPrim[VEL2-1];
    hTilde = rhoEpRef*UPrimRef[VEL3-1] + rhoEp*UPrim[VEL3-1];
  }

  return metricCoeff1*fTilde + metricCoeff2*gTilde + metricCoeff3*hTilde;
}

/***********************************************************************************************************************************
 * @brief Computes the Split-Flux retaining the formulation of Morinishi
 **********************************************************************************************************************************/
__device__ double SplitVolumeFluxMO_Device(int n, double* Uref, double* UPrimRef, double* U, double* UPrim, double* Mref, double* M)
{
  double fTilde, gTilde, hTilde;
  double metricCoeff1 = 0.5*(Mref[0]+M[0]);
  double metricCoeff2 = 0.5*(Mref[1]+M[1]);
  double metricCoeff3 = 0.5*(Mref[2]+M[2]);

  double rhoEpRef = Uref[ENER-1]-0.5*Uref[DENS-1]*(UPrimRef[VEL1-1]*UPrimRef[VEL1-1] + 
                    UPrimRef[VEL2-1]*UPrimRef[VEL2-1]+UPrimRef[VEL3-1]*UPrimRef[VEL3-1])+UPrimRef[PRES-1];
  double rhoEp = U[ENER-1]-0.5*U[DENS-1]*(UPrim[VEL1-1]*UPrim[VEL1-1]+UPrim[VEL2-1]*UPrim[VEL2-1] +
                  UPrim[VEL3-1]*UPrim[VEL3-1])+UPrim[PRES-1];

  if (n == DENS-1)
  {
    fTilde = Uref[MOM1-1]+U[MOM1-1];
    gTilde = Uref[MOM2-1]+U[MOM2-1];
    hTilde = Uref[MOM3-1]+U[MOM3-1];
  }
  else if (n == MOM1-1)
  {
    fTilde = 0.5*(Uref[MOM1-1]+U[MOM1-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]) + (UPrimRef[PRES-1]+UPrim[PRES-1]);
    gTilde = 0.5*(Uref[MOM2-1]+U[MOM2-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]);
    hTilde = 0.5*(Uref[MOM3-1]+U[MOM3-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]);
  }
  else if (n == MOM2-1)
  {
    fTilde = 0.5*(Uref[MOM1-1]+U[MOM1-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    gTilde = 0.5*(Uref[MOM2-1]+U[MOM2-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]) + (UPrimRef[PRES-1]+UPrim[PRES-1]);
    hTilde = 0.5*(Uref[MOM3-1]+U[MOM3-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);    
  }
  else if (n == MOM3-1)
  {
    fTilde = 0.5*(Uref[MOM1-1]+U[MOM1-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
    gTilde = 0.5*(Uref[MOM2-1]+U[MOM2-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
    hTilde = 0.5*(Uref[MOM3-1]+U[MOM3-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]) + (UPrimRef[PRES-1]+UPrim[PRES-1]);
  }
  else if (n == ENER-1)
  {
    fTilde = (rhoEpRef*UPrimRef[VEL1-1]+rhoEp*UPrim[VEL1-1]) +
              0.5*(Uref[MOM1-1]*UPrimRef[VEL1-1]+U[MOM1-1]*UPrim[VEL1-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]) +
              0.5*(Uref[MOM1-1]*UPrimRef[VEL2-1]+U[MOM1-1]*UPrim[VEL2-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]) +
              0.5*(Uref[MOM1-1]*UPrimRef[VEL3-1]+U[MOM1-1]*UPrim[VEL3-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]) -
              0.5*(Uref[MOM1-1]*UPrimRef[VEL1-1]*UPrimRef[VEL1-1]+U[MOM1-1]*UPrim[VEL1-1]*UPrim[VEL1-1])   -
              0.5*(Uref[MOM1-1]*UPrimRef[VEL2-1]*UPrimRef[VEL2-1]+U[MOM1-1]*UPrim[VEL2-1]*UPrim[VEL2-1])   -
              0.5*(Uref[MOM1-1]*UPrimRef[VEL3-1]*UPrimRef[VEL3-1]+U[MOM1-1]*UPrim[VEL3-1]*UPrim[VEL3-1]);
    gTilde = (rhoEpRef*UPrimRef[VEL2-1]+rhoEp*UPrim[VEL2-1]) +
              0.5*(Uref[MOM2-1]*UPrimRef[VEL1-1]+U[MOM2-1]*UPrim[VEL1-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]) +
              0.5*(Uref[MOM2-1]*UPrimRef[VEL2-1]+U[MOM2-1]*UPrim[VEL2-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]) +
              0.5*(Uref[MOM2-1]*UPrimRef[VEL3-1]+U[MOM2-1]*UPrim[VEL3-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]) -
              0.5*(Uref[MOM2-1]*UPrimRef[VEL1-1]*UPrimRef[VEL1-1]+U[MOM2-1]*UPrim[VEL1-1]*UPrim[VEL1-1])   -
              0.5*(Uref[MOM2-1]*UPrimRef[VEL2-1]*UPrimRef[VEL2-1]+U[MOM2-1]*UPrim[VEL2-1]*UPrim[VEL2-1])   -
              0.5*(Uref[MOM2-1]*UPrimRef[VEL3-1]*UPrimRef[VEL3-1]+U[MOM2-1]*UPrim[VEL3-1]*UPrim[VEL3-1]);
    hTilde = (rhoEpRef*UPrimRef[VEL3-1]+rhoEp*UPrim[VEL3-1]) +
              0.5*(Uref[MOM3-1]*UPrimRef[VEL1-1]+U[MOM3-1]*UPrim[VEL1-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]) +
              0.5*(Uref[MOM3-1]*UPrimRef[VEL2-1]+U[MOM3-1]*UPrim[VEL2-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]) +
              0.5*(Uref[MOM3-1]*UPrimRef[VEL3-1]+U[MOM3-1]*UPrim[VEL3-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]) -
              0.5*(Uref[MOM3-1]*UPrimRef[VEL1-1]*UPrimRef[VEL1-1]+U[MOM3-1]*UPrim[VEL1-1]*UPrim[VEL1-1])   -
              0.5*(Uref[MOM3-1]*UPrimRef[VEL2-1]*UPrimRef[VEL2-1]+U[MOM3-1]*UPrim[VEL2-1]*UPrim[VEL2-1])   -
              0.5*(Uref[MOM3-1]*UPrimRef[VEL3-1]*UPrimRef[VEL3-1]+U[MOM3-1]*UPrim[VEL3-1]*UPrim[VEL3-1]);
  }

  return metricCoeff1*fTilde + metricCoeff2*gTilde + metricCoeff3*hTilde;
}

/***********************************************************************************************************************************
 * @brief Computes the Split-Flux retaining the formulation of Ducros
 **********************************************************************************************************************************/
__device__ double SplitVolumeFluxDU_Device(int n, double* Uref, double* UPrimRef, double* U, double* UPrim, double* Mref, double* M)
{
  double fTilde, gTilde, hTilde;
  double metricCoeff1 = 0.5*(Mref[0]+M[0]);
  double metricCoeff2 = 0.5*(Mref[1]+M[1]);
  double metricCoeff3 = 0.5*(Mref[2]+M[2]);

  if (n == DENS-1)
  {
    fTilde = 0.5*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]);
    gTilde = 0.5*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    hTilde = 0.5*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
  }
  else if (n == MOM1-1)
  {
    fTilde = 0.5*(Uref[MOM1-1]+U[MOM1-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]) + (UPrimRef[PRES-1]+UPrim[PRES-1]);
    gTilde = 0.5*(Uref[MOM1-1]+U[MOM1-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    hTilde = 0.5*(Uref[MOM1-1]+U[MOM1-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
  }
  else if (n == MOM2-1)
  {
    fTilde = 0.5*(Uref[MOM2-1]+U[MOM2-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]);
    gTilde = 0.5*(Uref[MOM2-1]+U[MOM2-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]) + (UPrimRef[PRES-1]+UPrim[PRES-1]);
    hTilde = 0.5*(Uref[MOM2-1]+U[MOM2-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
  }
  else if (n == MOM3-1)
  {
    fTilde = 0.5*(Uref[MOM3-1]+U[MOM3-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]);
    gTilde = 0.5*(Uref[MOM3-1]+U[MOM3-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    hTilde = 0.5*(Uref[MOM3-1]+U[MOM3-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]) + (UPrimRef[PRES-1]+UPrim[PRES-1]);
  }
  else if (n == ENER-1)
  {
    fTilde = 0.5*(Uref[ENER-1]+U[ENER-1]+UPrimRef[PRES-1]+UPrim[PRES-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]);
    gTilde = 0.5*(Uref[ENER-1]+U[ENER-1]+UPrimRef[PRES-1]+UPrim[PRES-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    hTilde = 0.5*(Uref[ENER-1]+U[ENER-1]+UPrimRef[PRES-1]+UPrim[PRES-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
  }

  return metricCoeff1*fTilde + metricCoeff2*gTilde + metricCoeff3*hTilde;
}

/***********************************************************************************************************************************
 * @brief Computes the Split-Flux retaining the KEP formulation of Kennedy and Gruber
 **********************************************************************************************************************************/
__device__ double SplitVolumeFluxKG_Device(int n, double* Uref, double* UPrimRef, double* U, double* UPrim, double* Mref, double* M)
{
  double fTilde, gTilde, hTilde;
  double metricCoeff1 = 0.5*(Mref[0]+M[0]);
  double metricCoeff2 = 0.5*(Mref[1]+M[1]);
  double metricCoeff3 = 0.5*(Mref[2]+M[2]);

  // specific total energy
  double ERef = Uref[ENER-1]/Uref[DENS-1];
  double E    = U[ENER-1]/U[DENS-1];

  double tmpVal = 0.0;

  if (n == DENS-1)
  {
    fTilde = 0.5* (Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]);
    gTilde = 0.5 *(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    hTilde = 0.5 *(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
  }
  else if (n == MOM1-1)
  {
    tmpVal = (UPrimRef[VEL1-1]+UPrim[VEL1-1]);
    fTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*tmpVal*tmpVal + (UPrimRef[PRES-1]+UPrim[PRES-1]);
    gTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    hTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
  }
  else if (n == MOM2-1)
  {
    fTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    tmpVal = (UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    gTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*tmpVal*tmpVal + (UPrimRef[PRES-1]+UPrim[PRES-1]);
    hTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
  }
  else if (n == MOM3-1)
  {
    fTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
    gTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
    tmpVal = (UPrimRef[VEL3-1]+UPrim[VEL3-1]);
    hTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*tmpVal*tmpVal + (UPrimRef[PRES-1]+UPrim[PRES-1]);
  }
  else if (n == ENER-1)
  {
    fTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1])*(ERef+E) +
                0.5* (UPrimRef[PRES-1]+UPrim[PRES-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]);
    gTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1])*(ERef+E) +
                0.5* (UPrimRef[PRES-1]+UPrim[PRES-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    hTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1])*(ERef+E) +
                0.5 *(UPrimRef[PRES-1]+UPrim[PRES-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
  }

  return metricCoeff1*fTilde + metricCoeff2*gTilde + metricCoeff3*hTilde;
}

/***********************************************************************************************************************************
 * @brief Computes the Split-Flux retaining the KEP formulation of Pirozzoli
 **********************************************************************************************************************************/
__device__ double SplitVolumeFluxPI_Device(int n, double* Uref, double* UPrimRef, double* U, double* UPrim, double* Mref, double* M)
{
  double fTilde, gTilde, hTilde;
  double metricCoeff1 = 0.5*(Mref[0]+M[0]);
  double metricCoeff2 = 0.5*(Mref[1]+M[1]);
  double metricCoeff3 = 0.5*(Mref[2]+M[2]);

  // specific enthalpy, H=E+p/rho=(rhoE+p)/rho
  double HRef = (Uref[ENER-1]+UPrimRef[PRES-1])/Uref[DENS-1];
  double H    = (U[ENER-1]+UPrim[PRES-1])/U[DENS-1];

  double tmpVal = 0.0;

  if (n == DENS-1)
  {
    fTilde = 0.5 *(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1]);
    gTilde = 0.5 *(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    hTilde = 0.5 *(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
  }
  else if (n == MOM1-1)
  {
    tmpVal = (UPrimRef[VEL1-1]+UPrim[VEL1-1]);
    fTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*tmpVal*tmpVal + (UPrimRef[PRES-1]+UPrim[PRES-1]);
    gTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    hTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
  }
  else if (n == MOM2-1)
  {
    fTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    tmpVal = (UPrimRef[VEL2-1]+UPrim[VEL2-1]);
    gTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*tmpVal*tmpVal + (UPrimRef[PRES-1]+UPrim[PRES-1]);
    hTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
  }
  else if (n == MOM3-1)
  {
    fTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
    gTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1]);
    tmpVal = (UPrimRef[VEL3-1]+UPrim[VEL3-1]);
    hTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*tmpVal*tmpVal + (UPrimRef[PRES-1]+UPrim[PRES-1]);
  }
  else if (n == ENER-1)
  {
    fTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL1-1]+UPrim[VEL1-1])*(HRef+H);
    gTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL2-1]+UPrim[VEL2-1])*(HRef+H);
    hTilde = 0.25*(Uref[DENS-1]+U[DENS-1])*(UPrimRef[VEL3-1]+UPrim[VEL3-1])*(HRef+H);
  }

  return metricCoeff1*fTilde + metricCoeff2*gTilde + metricCoeff3*hTilde;
}

/***********************************************************************************************************************************
 * @brief Computes the Split-Flux retaining the entropy conserving (and formally KEP) formulation of Chandrashekar
 **********************************************************************************************************************************/
__device__ double SplitVolumeFluxCH_Device(int n, double* Uref, double* UPrimRef, double* U, double* UPrim, double* Mref, double* M)
{
  double fTilde, gTilde, hTilde;
  double metricCoeff1 = 0.5*(Mref[0]+M[0]);
  double metricCoeff2 = 0.5*(Mref[1]+M[1]);
  double metricCoeff3 = 0.5*(Mref[2]+M[2]);

  // average velocities
  double uMean = 0.5*(UPrimRef[VEL1-1] + UPrim[VEL1-1]);
  double vMean = 0.5*(UPrimRef[VEL2-1] + UPrim[VEL2-1]);
  double wMean = 0.5*(UPrimRef[VEL3-1] + UPrim[VEL3-1]);

  // inverse temperature
  double betaRef  = 0.5*Uref[DENS-1]/UPrimRef[PRES-1];
  double beta     = 0.5*U[DENS-1]/UPrim[PRES-1];

  // Density and inverse temperature logarithmic average
  double rhoLogMean = GetLogMean_Device(&Uref[DENS-1],&U[DENS-1]);
  double betaLogMean = GetLogMean_Device(&betaRef,&beta);

  // Average of pressure and specific enthalpy
  double pHatMean = 0.5*(Uref[DENS-1]+U[DENS-1])/(betaRef+beta);
  double HMean    = 0.5*d_sKappaM1/betaLogMean + pHatMean/rhoLogMean +
            0.5*(UPrimRef[VEL1-1]*UPrim[VEL1-1] + UPrimRef[VEL2-1]*UPrim[VEL2-1] + UPrimRef[VEL3-1]*UPrim[VEL3-1]);

  if (n == DENS-1)
  {
    fTilde = rhoLogMean*uMean;
    gTilde = rhoLogMean*vMean;
    hTilde = rhoLogMean*wMean;
  }
  else if (n == MOM1-1)
  {
    fTilde = rhoLogMean*uMean*uMean + pHatMean;
    gTilde = rhoLogMean*vMean*uMean;
    hTilde = rhoLogMean*wMean*uMean;
  }
  else if (n == MOM2-1)
  {
    fTilde = rhoLogMean*uMean*vMean;
    gTilde = rhoLogMean*vMean*vMean +pHatMean;
    hTilde = rhoLogMean*wMean*vMean;
  }
  else if (n == MOM3-1)
  {
    fTilde = rhoLogMean*uMean*wMean;
    gTilde = rhoLogMean*vMean*wMean;
    hTilde = rhoLogMean*wMean*wMean + pHatMean;
  }
  else if (n == ENER-1)
  {
    fTilde = rhoLogMean*HMean*uMean;
    gTilde = rhoLogMean*HMean*vMean;
    hTilde = rhoLogMean*HMean*wMean;
  }

  return 2.0*(metricCoeff1*fTilde + metricCoeff2*gTilde + metricCoeff3*hTilde);
}


//----------------------------------------------------------------------------------------------------------------------------------
// Split DG Surface Formulation Device Methods
//
// NOTE: All of these methods operate on a single variable at once. This is because it is easier to do than
//       return an array from a C function. It also works in a more clean way with the logic of the methods that
//       call these methods.
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Computes the surface flux for the split formulation retaining the standard NS-Equations
 **********************************************************************************************************************************/
__device__ double SplitSurfaceFluxSD_Device(int n, double* U_LL, double* U_RR)
{
  double Flux;

  if (n == DENS-1)
  {
    Flux = 0.5*(U_LL[EXT_MOM1-1]+U_RR[EXT_MOM1-1]);
  }
  else if (n == MOM1-1)
  {
    Flux = 0.5*(U_LL[EXT_MOM1-1]*U_LL[EXT_VEL1-1]+U_LL[EXT_PRES-1]+U_RR[EXT_MOM1-1]*U_RR[EXT_VEL1-1]+U_RR[EXT_PRES-1]);
  }
  else if (n == MOM2-1)
  {
    Flux = 0.5*(U_LL[EXT_MOM1-1]*U_LL[EXT_VEL2-1]+U_RR[EXT_MOM1-1]*U_RR[EXT_VEL2-1]);
  }
  else if (n == MOM3-1)
  {
    Flux = 0.5*(U_LL[EXT_MOM1-1]*U_LL[EXT_VEL3-1]+U_RR[EXT_MOM1-1]*U_RR[EXT_VEL3-1]);
  }
  else if (n == ENER-1)
  {
    Flux = 0.5*((U_LL[EXT_ENER-1]+U_LL[EXT_PRES-1])*U_LL[EXT_VEL1-1]+(U_RR[EXT_ENER-1]+U_RR[EXT_PRES-1])*U_RR[EXT_VEL1-1]);
  }

  return Flux;
}

/***********************************************************************************************************************************
 * @brief Computes the surface flux for the split formulation of Morinishi
 **********************************************************************************************************************************/
__device__ double SplitSurfaceFluxMO_Device(int n, double* U_LL, double* U_RR)
{
  double Flux;

  double rhoep_LL = U_LL[EXT_ENER-1]-0.5*U_LL[EXT_DENS-1]*(U_LL[EXT_VEL1-1]*U_LL[EXT_VEL1-1] +
                      U_LL[EXT_VEL2-1]*U_LL[EXT_VEL2-1]+U_LL[EXT_VEL3-1]*U_LL[EXT_VEL3-1])+U_LL[EXT_PRES-1];
  double rhoep_RR = U_RR[EXT_ENER-1]-0.5*U_RR[EXT_DENS-1]*(U_RR[EXT_VEL1-1]*U_RR[EXT_VEL1-1] +
                      U_RR[EXT_VEL2-1]*U_RR[EXT_VEL2-1]+U_RR[EXT_VEL3-1]*U_RR[EXT_VEL3-1])+U_RR[EXT_PRES-1];

  // compute flux
  if (n == DENS-1)
  {
    Flux = 0.5 *(U_LL[EXT_MOM1-1]+U_RR[EXT_MOM1-1]);
  }
  else if (n == MOM1-1)
  {
    Flux = 0.25*(U_LL[EXT_MOM1-1]+U_RR[EXT_MOM1-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]) + 0.5*(U_LL[EXT_PRES-1]+U_RR[EXT_PRES-1]);
  }
  else if (n == MOM2-1)
  {
    Flux = 0.25*(U_LL[EXT_MOM1-1]+U_RR[EXT_MOM1-1])*(U_LL[EXT_VEL2-1]+U_RR[EXT_VEL2-1]);
  }
  else if (n == MOM3-1)
  {
    Flux = 0.25*(U_LL[EXT_MOM1-1]+U_RR[EXT_MOM1-1])*(U_LL[EXT_VEL3-1]+U_RR[EXT_VEL3-1]);
  }
  else if (n == ENER-1)
  {
    Flux = 0.5 *(rhoep_LL*U_LL[EXT_VEL1-1]+rhoep_RR*U_RR[EXT_VEL1-1]) +
          0.25*(U_LL[EXT_MOM1-1]*U_LL[EXT_VEL1-1]+U_RR[EXT_MOM1-1]*U_RR[EXT_VEL1-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]) +
          0.25*(U_LL[EXT_MOM1-1]*U_LL[EXT_VEL2-1]+U_RR[EXT_MOM1-1]*U_RR[EXT_VEL2-1])*(U_LL[EXT_VEL2-1]+U_RR[EXT_VEL2-1]) +
          0.25*(U_LL[EXT_MOM1-1]*U_LL[EXT_VEL3-1]+U_RR[EXT_MOM1-1]*U_RR[EXT_VEL3-1])*(U_LL[EXT_VEL3-1]+U_RR[EXT_VEL3-1]) -
          0.25*(U_LL[EXT_MOM1-1]*U_LL[EXT_VEL1-1]*U_LL[EXT_VEL1-1]+U_RR[EXT_MOM1-1]*U_RR[EXT_VEL1-1]*U_RR[EXT_VEL1-1]) -
          0.25*(U_LL[EXT_MOM1-1]*U_LL[EXT_VEL2-1]*U_LL[EXT_VEL2-1]+U_RR[EXT_MOM1-1]*U_RR[EXT_VEL2-1]*U_RR[EXT_VEL2-1]) -
          0.25*(U_LL[EXT_MOM1-1]*U_LL[EXT_VEL3-1]*U_LL[EXT_VEL3-1]+U_RR[EXT_MOM1-1]*U_RR[EXT_VEL3-1]*U_RR[EXT_VEL3-1]);
  }

  return Flux;
}

/***********************************************************************************************************************************
 * @brief Computes the surface flux for the split formulation of Ducros
 **********************************************************************************************************************************/
__device__ double SplitSurfaceFluxDU_Device(int n, double* U_LL, double* U_RR)
{
  double Flux;

  if (n == DENS-1)
  {
    Flux = 0.25*(U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]);
  }
  else if (n == MOM1-1)
  {
    Flux = 0.25*(U_LL[EXT_MOM1-1]+U_RR[EXT_MOM1-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]) + 0.5*(U_LL[EXT_PRES-1]+U_RR[EXT_PRES-1]);
  }
  else if (n == MOM2-1)
  {
    Flux = 0.25*(U_LL[EXT_MOM2-1]+U_RR[EXT_MOM2-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]);
  }
  else if (n == MOM3-1)
  {
    Flux = 0.25*(U_LL[EXT_MOM3-1]+U_RR[EXT_MOM3-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]);
  }
  else if (n == ENER-1)
  {
    Flux = 0.25*(U_LL[EXT_ENER-1]+U_RR[EXT_ENER-1]+U_LL[EXT_PRES-1]+U_RR[EXT_PRES-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]);
  }

  return Flux;
}

/***********************************************************************************************************************************
 * @brief Computes the surface flux for the split formulation of Kennedy and Gruber
 **********************************************************************************************************************************/
__device__ double SplitSurfaceFluxKG_Device(int n, double* U_LL, double* U_RR)
{
  double Flux;
  double tmpVal = 0.0;

  // specific total energy
  double E_LL = U_LL[EXT_ENER-1]/U_LL[EXT_DENS-1];
  double E_RR = U_RR[EXT_ENER-1]/U_RR[EXT_DENS-1];

  if (n == DENS-1)
  {
    Flux = 0.25* (U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]);
  }
  else if (n == MOM1-1)
  {
    tmpVal = (U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]);
    Flux = 0.125*(U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])*tmpVal*tmpVal + 0.5*(U_LL[EXT_PRES-1]+U_RR[EXT_PRES-1]);
  }
  else if (n == MOM2-1)
  {
    Flux = 0.125*(U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1])*(U_LL[EXT_VEL2-1]+U_RR[EXT_VEL2-1]);
  }
  else if (n == MOM3-1)
  {
    Flux = 0.125*(U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1])*(U_LL[EXT_VEL3-1]+U_RR[EXT_VEL3-1]);
  }
  else if (n == ENER-1)
  {
    Flux = 0.125*(U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])*(E_LL+E_RR)*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]) +
          0.25 *(U_LL[EXT_PRES-1]+U_RR[EXT_PRES-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]);
  }

  return Flux;
}

/***********************************************************************************************************************************
 * @brief Computes the surface flux for the split formulation of Pirozzoli
 **********************************************************************************************************************************/
__device__ double SplitSurfaceFluxPI_Device(int n, double* U_LL, double* U_RR)
{
  double Flux;
  double tmpVal = 0.0;

  // specific energy, H=E+p/rho=(rhoE+p)/rho
  double H_LL = (U_LL[EXT_ENER-1]+U_LL[EXT_PRES-1])/U_LL[EXT_DENS-1];
  double H_RR = (U_RR[EXT_ENER-1]+U_RR[EXT_PRES-1])/U_RR[EXT_DENS-1];

  // compute flux
  if (n == DENS-1)
  {
    Flux = 0.25* (U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]);
  }
  else if (n == MOM1-1)
  {
    tmpVal = (U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]);
    Flux = 0.125*(U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])*tmpVal*tmpVal + 0.5*(U_LL[EXT_PRES-1]+U_RR[EXT_PRES-1]);
  }
  else if (n == MOM2-1)
  {
    Flux = 0.125*(U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1])*(U_LL[EXT_VEL2-1]+U_RR[EXT_VEL2-1]);
  }
  else if (n == MOM3-1)
  {
    Flux = 0.125*(U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1])*(U_LL[EXT_VEL3-1]+U_RR[EXT_VEL3-1]);
  }
  else if (n == ENER-1)
  {
    Flux = 0.125*(U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])*(H_LL+H_RR)*(U_LL[EXT_VEL1-1]+U_RR[EXT_VEL1-1]);
  }

  return Flux;
}

/***********************************************************************************************************************************
 * @brief Computes the surface flux for the entropy conserving formulation of Chandrashekar.
 **********************************************************************************************************************************/
__device__ double SplitSurfaceFluxCH_Device(int n, double* U_LL, double* U_RR)
{
  double Flux;
  double betaLogMean = 0.0;
  double rhoLogMean = 0.0;

  // average velocities
  double uMean = 0.5*(U_LL[EXT_VEL1-1] + U_RR[EXT_VEL1-1]);
  double vMean = 0.5*(U_LL[EXT_VEL2-1] + U_RR[EXT_VEL2-1]);
  double wMean = 0.5*(U_LL[EXT_VEL3-1] + U_RR[EXT_VEL3-1]);

  // inverse temperature
  double beta_LL = 0.5*U_LL[EXT_DENS-1]/U_LL[EXT_PRES-1];
  double beta_RR = 0.5*U_RR[EXT_DENS-1]/U_RR[EXT_PRES-1];

  // average pressure, enthalpy, density and inverse temperature
  // logarithmic mean
  rhoLogMean  = GetLogMean_Device(&U_LL[EXT_DENS-1],&U_RR[EXT_DENS-1]);
  betaLogMean = GetLogMean_Device(&beta_LL,&beta_RR);
  // "standard" average
  double pHatMean = 0.5*(U_LL[EXT_DENS-1]+U_RR[EXT_DENS-1])/(beta_LL+beta_RR);
  double HMean    = 0.5*d_sKappaM1/betaLogMean + pHatMean/rhoLogMean +
            0.5*(U_LL[EXT_VEL1-1]*U_RR[EXT_VEL1-1] + U_LL[EXT_VEL2-1]*U_RR[EXT_VEL2-1] + U_LL[EXT_VEL3-1]*U_RR[EXT_VEL3-1]);

  //compute flux
  if (n == DENS-1)
  {
    Flux = rhoLogMean*uMean;
  }
  else if (n == MOM1-1)
  {
    Flux = rhoLogMean*uMean*uMean + pHatMean;
  }
  else if (n == MOM2-1)
  {
    Flux = rhoLogMean*uMean*vMean;
  }
  else if (n == MOM3-1)
  {
    Flux = rhoLogMean*uMean*wMean;
  }
  else if (n == ENER-1)
  {
    Flux = rhoLogMean*uMean*HMean;
  }

  return Flux;
}

//----------------------------------------------------------------------------------------------------------------------------------
// Initialization functions
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Initialize function pointers for Split DG methods on the device side. This method runs on the DEVICE only
 * @param splitDGOpt User chosen option for the Split DG formulation from the input file
 * @remark We trust that illegal value checking is done by the CPU before this kernel is ever called, so we don't check here.
 **********************************************************************************************************************************/
__global__ void InitSplitDG_Kernel(int splitDGOpt)
{
  switch (splitDGOpt)
  {
    case PRM_SPLITDG_SD:
      SplitDGVolume_Pointer_Device  = SplitVolumeFluxSD_Device;
      SplitDGSurface_Pointer_Device = SplitSurfaceFluxSD_Device;
      break;
    case PRM_SPLITDG_MO:
      SplitDGVolume_Pointer_Device  = SplitVolumeFluxMO_Device;
      SplitDGSurface_Pointer_Device = SplitSurfaceFluxMO_Device;
      break;
    case PRM_SPLITDG_DU:
      SplitDGVolume_Pointer_Device  = SplitVolumeFluxDU_Device;
      SplitDGSurface_Pointer_Device = SplitSurfaceFluxDU_Device;
      break;
    case PRM_SPLITDG_KG:
      SplitDGVolume_Pointer_Device  = SplitVolumeFluxKG_Device;
      SplitDGSurface_Pointer_Device = SplitSurfaceFluxKG_Device;
      break;
    case PRM_SPLITDG_PI:
      SplitDGVolume_Pointer_Device  = SplitVolumeFluxPI_Device;
      SplitDGSurface_Pointer_Device = SplitSurfaceFluxPI_Device;
      break;
    case PRM_SPLITDG_CH:
      SplitDGVolume_Pointer_Device  = SplitVolumeFluxCH_Device;
      SplitDGSurface_Pointer_Device = SplitSurfaceFluxCH_Device;
      break;
  }
}

/***********************************************************************************************************************************
 * @brief Entry point for the device side initialization of the SplitFlux function pointers
 **********************************************************************************************************************************/
void InitSplitDG_Device(int splitDGOpt)
{
  INVOKE_KERNEL( InitSplitDG_Kernel, 1, 1, 0, streams[0], splitDGOpt );
}




