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
  void ConsToPrim_Volume_Device(int primKey, int consKey, int Nloc, double kappaM1, double R, int nElems, int streamID);
  void ConsToPrim_Elem_Device(int primKey, int consKey, int Nloc, double kappaM1, double R, int elemIdx);
  void ConsToPrim_Side_Device(int primKey, int consKey, int Nloc, double kappaM1, double R, int elemIdx, int sideIdx);
  void ConsToPrim_Sides_Device(int primKey, int consKey, int Nloc, double kappaM1, double R, int elemIdx, int sideIdx, int nSides, int streamID);
  void ConsToPrim_Point_Device(int primKey, int consKey, int Nloc, double kappaM1, double R, int elemIdx, int i, int j, int k);

  void ConsToEntropy_Volume_Device(int Nloc, int nElems, int d_entropy, int d_cons, double kappa, double kappaM1, double sKappaM1, int streamID);
  void EntropyToCons_Side_Device(int Nloc, int sideID, int elemID, int d_entropy, int d_cons, double kappa, double kappaM1, double sKappaM1);

  void InitEOS_Device(double* EOS_Vars);
}

//----------------------------------------------------------------------------------------------------------------------------------
// Any device copies of variables from MOD_EOS_Vars -- Functionally replaces EOS_vars array in GLXI v1.0
//----------------------------------------------------------------------------------------------------------------------------------
__device__ double d_Kappa;
__device__ double d_KappaM1;
__device__ double d_sKappaM1;
__device__ double d_KappaP1;
__device__ double d_sKappaP1;
__device__ double d_R;

__device__ double d_cp;
__device__ double d_cv;

#if PARABOLIC
__device__ double d_mu0;
__device__ double d_Pr;
__device__ double d_KappasPr;
#if PP_VISC==1
__device__ double d_Ts;
__device__ double d_cSuth;
#endif
#if (PP_VISC==1) || (PP_VISC==2)
__device__ double d_Tref;
__device__ double d_ExpoSuth;
#endif
#endif /*PARABOLIC*/


//##################################################################################################################################
// ConsToPrim
//##################################################################################################################################
/***********************************************************************************************************************************
 * @brief Device backend for transformation for core variable transformation for ConsToPrim
 * @param prim Pointer to DOF in primitive variables array assigned to this device thread
 * @param cons Pointer to DOF in conservative variables array assigned to this device thread
 * @param kappaM1 Specific heat ratio - 1.0
 * @param R Gas constant
 **********************************************************************************************************************************/
__device__ void ConsToPrimCore_Device(double* prim, double* cons, double kappaM1, double R)
{
    double velSum;
    
    double sRho = 1.0/cons[DENS-1];

    // density	
    prim[DENS-1] = cons[DENS-1];

    // velocities
    prim[VEL1-1] = cons[MOM1-1]*sRho;
    prim[VEL2-1] = cons[MOM2-1]*sRho;
    prim[VEL3-1] = cons[MOM3-1]*sRho;

    // pressure
    velSum = cons[MOM1-1]*prim[VEL1-1] + cons[MOM2-1]*prim[VEL2-1] + cons[MOM3-1]*prim[VEL3-1];
    prim[PRES-1] = kappaM1*(cons[ENER-1] - 0.5*velSum);

    // temperature
    prim[TEMP-1] = prim[PRES-1]*sRho / R;
}

//----------------------------------------------------------------------------------------------------------------------------------
// ConsToPrim Kernels
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Device kernel for Volume overload of ConsToPrim
 * @param nDOF Total number of degrees of freedom the kernel will operate on during this call
 * @param prim Pointer to device copty of prim var data
 * @param cons Pointer to device copy of cons var data
 * @param kappaM1 Specific head ration - 1.0
 * @param R Gas constant
 **********************************************************************************************************************************/
__global__ void __launch_bounds__(256, 1) ConsToPrim_Volume_Kernel(int nDOF, double* prim, double* cons, double kappaM1, double R)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  int i = threadID*PP_nVar;
  int j = threadID*PP_nVarPrim;
  
  if (threadID < nDOF)
  {
    ConsToPrimCore_Device(&prim[j], &cons[i], kappaM1, R);
  }
}

/***********************************************************************************************************************************
 * @brief Device kernel for Element overload of ConsToPrim
 * @param nDOF Total number of degrees of freedom the kernel will operate on during this call
 * @param prim Pointer to device copty of prim var data
 * @param cons Pointer to device copy of cons var data
 * @param kappaM1 Specific head ration - 1.0
 * @param R Gas constant
 * @param con_offset Offset in bytes to the starting index of the desired element in the cons array
 * @param prim_offset Offset in bytes to the starting index of the desired element in the prim array
 **********************************************************************************************************************************/
__global__ void ConsToPrim_Elem_Kernel(int nDOF, double* prim, double* cons, double kappaM1, double R, int con_offset, int prim_offset)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  // Need to offset the start of the calculation to the element we are working on
  int i = con_offset + threadID*PP_nVar;
  int j = prim_offset + threadID*PP_nVarPrim;
  
  if (threadID < nDOF)
  {
    ConsToPrimCore_Device(&prim[j], &cons[i], kappaM1, R);
  }
}

/***********************************************************************************************************************************
 * @brief Device kernel for Side overload of ConsToPrim
 * @param nDOF Total number of degrees of freedom the kernel will operate on during this call
 * @param prim Pointer to device copty of prim var data
 * @param cons Pointer to device copy of cons var data
 * @param kappaM1 Specific head ration - 1.0
 * @param R Gas constant
 * @param con_offset Offset in bytes to the starting index of the desired side in the cons array
 * @param prim_offset Offset in bytes to the starting index of the desired side in the prim array
 **********************************************************************************************************************************/
__global__ void ConsToPrim_Side_Kernel(int nDOF, double* prim, double* cons, double kappaM1, double R, int con_offset, int prim_offset)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  int i = con_offset + threadID*PP_nVar;
  int j = prim_offset + threadID*PP_nVarPrim;
  
  if (threadID < nDOF)
  {
    ConsToPrimCore_Device(&prim[j], &cons[i], kappaM1, R);
  }
}

/***********************************************************************************************************************************
 * @brief Device kernel for Point overload of ConsToPrim
 * @param nDOF Total number of degrees of freedom the kernel will operate on during this call
 * @param prim Pointer to device copty of prim var data
 * @param cons Pointer to device copy of cons var data
 * @param kappaM1 Specific head ration - 1.0
 * @param R Gas constant
 * @param con_offset Offset in bytes to the starting index of the desired DOF in the cons array
 * @param prim_offset Offset in bytes to the starting index of the desired DOF in the prim array
 **********************************************************************************************************************************/
__global__ void ConsToPrim_Point_Kernel(double* prim, double* cons, double kappaM1, double R, int con_offset, int prim_offset)
{
  // Kernel only ever called with a single thread, so threadID = 0 always
  ConsToPrimCore_Device(&prim[prim_offset], &cons[con_offset], kappaM1, R);
}


//----------------------------------------------------------------------------------------------------------------------------------
// ConsToPrim HOST Interfaces - Call kernels from HOST
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Interface entry point method for Volume overload of ConsToPrim
 * @param primKey Hash key string to prim var device array
 * @param consKey Hash key string to cons var device array
 * @param Nloc Local ploynomial order
 * @param kappaM1 Specific head ration - 1.0
 * @param R Gas constant
 * @param nElems Number of elems with in the volume to convert
 * @param streamID Index of the stream to use for this call
 **********************************************************************************************************************************/
void ConsToPrim_Volume_Device(int primKey, int consKey, int Nloc, double kappaM1, double R, int nElems, int streamID)
{
#ifdef USE_NVTX
  nvtxRangePushA("ConsToPrim_Volume");
#endif

  int nDOF = (Nloc+1)*(Nloc+1)*(Nloc+1)*nElems;

  INVOKE_KERNEL( ConsToPrim_Volume_Kernel, nDOF/256+1, 256, 0, streams[streamID],nDOF,(double*)DeviceVars[primKey],(double*)DeviceVars[consKey],kappaM1,R );

#ifdef USE_NVTX
  nvtxRangePop();
#endif
}

/***********************************************************************************************************************************
 * @brief Interface entry point method for Element overload of ConsToPrim
 * @param primKey Hash key string to prim var device array
 * @param consKey Hash key string to cons var device array
 * @param Nloc Local ploynomial order
 * @param kappaM1 Specific head ration - 1.0
 * @param R Gas constant
 * @param elemIdx Index of the element that will be converted
 **********************************************************************************************************************************/
void ConsToPrim_Elem_Device(int primKey, int consKey, int Nloc, double kappaM1, double R, int elemIdx)
{
  int nDOF = (Nloc+1)*(Nloc+1)*(Nloc+1);
  int con_offset = FindElemOffset(PP_nVar, Nloc, elemIdx);
  int prim_offset = FindElemOffset(PP_nVarPrim, Nloc, elemIdx);

  INVOKE_KERNEL( ConsToPrim_Elem_Kernel, nDOF/256+1, 256, 0, streams[0], nDOF,(double*)DeviceVars[primKey],(double*)DeviceVars[consKey],kappaM1,R,con_offset,prim_offset );
}

/***********************************************************************************************************************************
 * @brief Interface entry point method for Side overload of ConsToPrim
 * @param primKey Hash key string to prim var device array
 * @param consKey Hash key string to cons var device array
 * @param Nloc Local ploynomial order
 * @param kappaM1 Specific head ration - 1.0
 * @param R Gas constant
 * @param elemIdx Index of the element that will be converted (set to 0 for sides only arrays)
 * @param sideIdx Index of the side that will be converted
 **********************************************************************************************************************************/
void ConsToPrim_Side_Device(int primKey, int consKey, int Nloc, double kappaM1, double R, int elemIdx, int sideIdx)
{
  int nDOF = (Nloc+1)*(Nloc+1);
  int con_offset = FindSideOffset(PP_nVar, Nloc, elemIdx, sideIdx);
  int prim_offset = FindSideOffset(PP_nVarPrim, Nloc, elemIdx, sideIdx);

  INVOKE_KERNEL( ConsToPrim_Side_Kernel, nDOF/64+1, 64, 0, streams[0], nDOF,(double*)DeviceVars[primKey],(double*)DeviceVars[consKey],kappaM1,R,con_offset,prim_offset );
}

/***********************************************************************************************************************************
 * @brief Interface entry point method for Sides overload of ConsToPrim
 * @param primKey Hash key string to prim var device array
 * @param consKey Hash key string to cons var device array
 * @param Nloc Local ploynomial order
 * @param kappaM1 Specific head ration - 1.0
 * @param R Gas constant
 * @param elemIdx Index of the element that will be converted (set to 0 for sides only arrays)
 * @param sideIdx Index of the side that will be converted
 * @param nSides Number of sides to convert within the element
 **********************************************************************************************************************************/
void ConsToPrim_Sides_Device(int primKey, int consKey, int Nloc, double kappaM1, double R, int elemIdx, int sideIdx, int nSides, int streamID)
{
  int nDOF = (Nloc+1)*(Nloc+1)*nSides;
  int con_offset = FindSideOffset(PP_nVar, Nloc, elemIdx, sideIdx);
  int prim_offset = FindSideOffset(PP_nVarPrim, Nloc, elemIdx, sideIdx);

  INVOKE_KERNEL( ConsToPrim_Side_Kernel, nDOF/256+1, 256, 0, streams[streamID], nDOF,(double*)DeviceVars[primKey],(double*)DeviceVars[consKey],kappaM1,R,con_offset,prim_offset );
}

/***********************************************************************************************************************************
 * @brief Interface entry point method for Point overload of ConsToPrim
 * @param primKey Hash key string to prim var device array
 * @param consKey Hash key string to cons var device array
 * @param Nloc Local ploynomial order
 * @param kappaM1 Specific head ration - 1.0
 * @param R Gas constant
 * @param elemIdx Index of the element that will be converted
 * @param i 1st dim index of DOF within element
 * @param j 2nd dim index of DOF within element
 * @param k 3rd dim index of DOF within element
 **********************************************************************************************************************************/
void ConsToPrim_Point_Device(int primKey, int consKey, int Nloc, double kappaM1, double R, int elemIdx, int i, int j, int k)
{
  int con_offset = FindPointOffset(PP_nVar, Nloc, elemIdx, i, j, k);
  int prim_offset = FindPointOffset(PP_nVarPrim, Nloc, elemIdx, i, j, k);
  INVOKE_KERNEL( ConsToPrim_Point_Kernel, 1, 1, 0, streams[0], (double*)DeviceVars[primKey],(double*)DeviceVars[consKey],kappaM1,R,con_offset,prim_offset );
}

//##################################################################################################################################
// ConsToEntropy
//##################################################################################################################################
/***********************************************************************************************************************************
 * @brief Conversion of conservative variables to entropy variables for a single DOF
 * @param entropy Pointer to the DOF assigned to this device thread in entropy variables array
 * @param cons Pointer to the DOF assigned to this device thread in conservative variables array
 * @param kappa Specific heat ratio
 * @param kappaM1 kappa - 1.0
 * @param sKappaM1 1.0/kappaM1
 **********************************************************************************************************************************/
__device__ void ConsToEntropyCore_Device(double* entropy, double* cons, double kappa, double kappaM1, double sKappaM1)
{
  double vel[3];
  double p,s,rho_p;

  for (int n = 0; n < 3; n++)
  {
    vel[n] = cons[MOM1-1+n] / cons[DENS-1];
  }
  p = kappaM1*(cons[ENER-1]-0.5*(cons[MOM1-1]*vel[0] + cons[MOM2-1]*vel[1] + cons[MOM3-1]*vel[2]));
  // entropy: log(p) - eq.γ * log(ρ)
  s = log(p) - kappa*log(cons[DENS-1]);
  rho_p  = cons[DENS-1]/p;

  // Convert to entropy variables
  // (γ - s) / (γ - 1) - (ρu^2 + ρv^2 + ρw^2) / ρ / 2p,
  entropy[DENS-1] = (kappa-s)*sKappaM1 - 0.5*rho_p*(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
  // ρu / p
  entropy[MOM1-1] = rho_p*vel[0];
  entropy[MOM2-1] = rho_p*vel[1];
  entropy[MOM3-1] = rho_p*vel[2];
  // -ρ / p
  entropy[ENER-1] = -rho_p;
}

//----------------------------------------------------------------------------------------------------------------------------------
// ConsToEntropy Kernels
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Conversion of conservative variables to entropy variables for the entire volume
 * @param nDOF Total number of DOFs in the volume being converted
 * @param entropy Pointer to the entire device copy of the entropy variables array
 * @param cons Pointer to the entire device copy of the conservative variables array
 * @param kappa Specific heat ratio
 * @param kappaM1 kappa - 1.0
 * @param sKappaM1 1.0/kappaM1
 **********************************************************************************************************************************/
__global__ void __launch_bounds__(256, 1) ConsToEntropy_Volume_Kernel(int nDOF, double* entropy, double* cons, double kappa, double kappaM1, double sKappaM1)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  int i = threadID*PP_nVar;
  
  if (threadID < nDOF)
    ConsToEntropyCore_Device(&entropy[i], &cons[i], kappa, kappaM1, sKappaM1);
}


//----------------------------------------------------------------------------------------------------------------------------------
// ConsToEntropy HOST Interfaces - Call kernels from HOST
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Entry point for ConsToEntropy volume kernel
 * @param nDOF Total number of DOFs in the volume being converted
 * @param nElems Total number of elements in the volume being converted
 * @param d_entropy Device variable key for the entropy variables array
 * @param d_cons Device variable key for the conservative variables array
 * @param kappa Specific heat ratio
 * @param kappaM1 kappa - 1.0
 * @param sKappaM1 1.0/kappaM1
 * @param streamID ID of device stream to run the kernel with
 **********************************************************************************************************************************/
void ConsToEntropy_Volume_Device(int Nloc, int nElems, int d_entropy, int d_cons, double kappa, double kappaM1, double sKappaM1, int streamID)
{
#ifdef USE_NVTX
  nvtxRangePushA("ConsToEntropy_Volume");
#endif

  int nDOF = (Nloc+1)*(Nloc+1)*(Nloc+1)*nElems;

  INVOKE_KERNEL( ConsToEntropy_Volume_Kernel, nDOF/256+1, 256, 0, streams[streamID], 
                            nDOF,(double*)DeviceVars[d_entropy],(double*)DeviceVars[d_cons],kappa,kappaM1,sKappaM1);

#ifdef USE_NVTX
  nvtxRangePop();
#endif
}

//##################################################################################################################################
// EntropyToCons
//##################################################################################################################################
/***********************************************************************************************************************************
 * @brief Conversion of entropy variables to conservative variables for a single DOF
 * @param entropy Pointer to the DOF assigned to this device thread in entropy variables array
 * @param cons Pointer to the DOF assigned to this device thread in conservative variables array
 * @param kappa Specific heat ratio
 * @param kappaM1 kappa - 1.0
 * @param sKappaM1 1.0/kappaM1
 **********************************************************************************************************************************/
__device__ void EntropyToConsCore_Device(double* entropy, double* cons, double kappa, double kappaM1, double sKappaM1)
{
  double ent2[PP_nVar];
  double s, rhoe;
  double sumMom2 = 0.0;

  for (int n = 0; n < PP_nVar; n++)
  {
    ent2[n] = entropy[n]*kappaM1;
  }
  sumMom2 = ent2[MOM1-1]*ent2[MOM1-1] + ent2[MOM2-1]*ent2[MOM2-1] + ent2[MOM3-1]*ent2[MOM3-1];
  s = kappa - ent2[DENS-1] + 0.5*sumMom2 / ent2[ENER-1];
  rhoe = pow(kappaM1 / pow(-ent2[ENER-1],kappa), sKappaM1) * exp(-s*sKappaM1);

  cons[DENS-1] = -rhoe * ent2[ENER-1]; // ρ = -p * W[5]
  cons[MOM1-1] =  rhoe * ent2[MOM1-1];
  cons[MOM2-1] =  rhoe * ent2[MOM2-1];
  cons[MOM3-1] =  rhoe * ent2[MOM3-1];
  cons[ENER-1] =  rhoe * (1.0 - sumMom2 * 0.5/ent2[ENER-1]);
}


//----------------------------------------------------------------------------------------------------------------------------------
// EntropyToCons Kernels
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Conversion of conservative variables to entropy variables for a given side
 * @param nDOF Total number of DOFs in the side being converted
 * @param offset Memory offset in bytes to the working side
 * @param entropy Pointer to the entire device copy of the entropy variables array
 * @param cons Pointer to the entire device copy of the conservative variables array
 * @param kappa Specific heat ratio
 * @param kappaM1 kappa - 1.0
 * @param sKappaM1 1.0/kappaM1
 **********************************************************************************************************************************/
__global__ void EntropyToCons_Side_Kernel(int nDOF, int offset, double* entropy, double* cons, double kappa, double kappaM1, double sKappaM1)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  // Get offset to this thread's DOF -- offset to side + offset to DOF in side
  int i = offset + threadID*PP_nVar;
  
  if (threadID < nDOF)
    EntropyToConsCore_Device(&entropy[i], &cons[i], kappa, kappaM1, sKappaM1);
}

//----------------------------------------------------------------------------------------------------------------------------------
// EntropyToCons HOST Interfaces - Call kernels from HOST
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Entry point for EntropyToCons side kernel
 * @param nDOF Total number of DOFs in the volume being converted
 * @param sideID Index of the side we would like to convert properties for
 * @param elemID Index of the element associated with the sideID. 0 for sides only arrays.
 * @param d_entropy Device variable key for the entropy variables array
 * @param d_cons Device variable key for the conservative variables array
 * @param kappa Specific heat ratio
 * @param kappaM1 kappa - 1.0
 * @param sKappaM1 1.0/kappaM1
 **********************************************************************************************************************************/
void EntropyToCons_Side_Device(int Nloc, int sideID, int elemID, int d_entropy, int d_cons, double kappa, double kappaM1, double sKappaM1)
{
  int nDOF = (Nloc+1)*(Nloc+1);
  int offset = FindSideOffset(PP_nVar, Nloc, elemID, sideID-1);

  INVOKE_KERNEL( EntropyToCons_Side_Kernel, nDOF/64+1, 64, 0, streams[0], 
                    nDOF,offset,(double*)DeviceVars[d_entropy],(double*)DeviceVars[d_cons],kappa,kappaM1,sKappaM1 );
}



//----------------------------------------------------------------------------------------------------------------------------------
// Device init method
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Init method for any device side stuff for EOS vars so we don't have to pass into kernels.
 * @param EOS_Vars Pointer to Fortran host copy of EOS_Vars array
 **********************************************************************************************************************************/
void InitEOS_Device(double* EOS_Vars)
{
  // Now copy the calculated values to the device.
#if (USE_ACCEL == ACCEL_CUDA)
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_Kappa, &EOS_Vars[EOS_KAPPA-1], sizeof(double), 0))
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_KappaM1, &EOS_Vars[EOS_KAPPAM1-1], sizeof(double), 0))
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_sKappaM1, &EOS_Vars[EOS_SKAPPAM1-1], sizeof(double), 0))
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_KappaP1, &EOS_Vars[EOS_KAPPAP1-1], sizeof(double), 0))
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_sKappaP1, &EOS_Vars[EOS_SKAPPAP1-1], sizeof(double), 0))
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_R, &EOS_Vars[EOS_R-1], sizeof(double), 0))

  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_cp, &EOS_Vars[EOS_CP-1], sizeof(double), 0))
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_cv, &EOS_Vars[EOS_CV-1], sizeof(double), 0))

#if PARABOLIC
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_mu0, &EOS_Vars[EOS_MU0-1], sizeof(double), 0))
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_Pr, &EOS_Vars[EOS_PR-1], sizeof(double), 0))
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_KappasPr, &EOS_Vars[EOS_KAPPASPR-1], sizeof(double), 0))
#if PP_VISC==1
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_Ts, &EOS_Vars[EOS_TS-1], sizeof(double), 0))
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_cSuth, &EOS_Vars[EOS_CSUTH-1], sizeof(double), 0))
#endif
#if (PP_VISC==1) || (PP_VISC==2)
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_Tref, &EOS_Vars[EOS_TREF-1], sizeof(double), 0))
  DEVICE_ERR_CHECK( cudaMemcpyToSymbol(d_ExpoSuth, &EOS_Vars[EOS_EXPOSUTH-1], sizeof(double), 0))
#endif
#endif /*PARABOLIC*/
#elif (USE_ACCEL == ACCEL_HIP)
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_Kappa), &EOS_Vars[EOS_KAPPA-1], sizeof(double), 0, hipMemcpyHostToDevice))
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_KappaM1), &EOS_Vars[EOS_KAPPAM1-1], sizeof(double), 0, hipMemcpyHostToDevice))
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_sKappaM1), &EOS_Vars[EOS_SKAPPAM1-1], sizeof(double), 0, hipMemcpyHostToDevice))
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_KappaP1), &EOS_Vars[EOS_KAPPAP1-1], sizeof(double), 0, hipMemcpyHostToDevice))
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_sKappaP1), &EOS_Vars[EOS_SKAPPAP1-1], sizeof(double), 0, hipMemcpyHostToDevice))
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_R), &EOS_Vars[EOS_R-1], sizeof(double), 0, hipMemcpyHostToDevice))

  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_cp), &EOS_Vars[EOS_CP-1], sizeof(double), 0, hipMemcpyHostToDevice))
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_cv), &EOS_Vars[EOS_CV-1], sizeof(double), 0, hipMemcpyHostToDevice))

#if PARABOLIC
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_mu0), &EOS_Vars[EOS_MU0-1], sizeof(double), 0, hipMemcpyHostToDevice))
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_Pr), &EOS_Vars[EOS_PR-1], sizeof(double), 0, hipMemcpyHostToDevice))
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_KappasPr), &EOS_Vars[EOS_KAPPASPR-1], sizeof(double), 0, hipMemcpyHostToDevice))
#if PP_VISC==1
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_Ts), &EOS_Vars[EOS_TS-1], sizeof(double), 0, hipMemcpyHostToDevice))
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_cSuth), &EOS_Vars[EOS_CSUTH-1], sizeof(double), 0, hipMemcpyHostToDevice))
#endif
#if (PP_VISC==1) || (PP_VISC==2)
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_Tref), &EOS_Vars[EOS_TREF-1], sizeof(double), 0, hipMemcpyHostToDevice))
  DEVICE_ERR_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(d_ExpoSuth), &EOS_Vars[EOS_EXPOSUTH-1], sizeof(double), 0, hipMemcpyHostToDevice))
#endif
#endif /*PARABOLIC*/
#endif /* USE_ACCEL */
}

//----------------------------------------------------------------------------------------------------------------------------------
// Device function renderings of _HE Macros from eos.h that compute physical quantities from conservative variables or extended variables
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Device function rendering of the VELOCITY_HE macro from eos.h. Coverts momentum to velocity
 * @param vel Client side allocated array in which the resulting converted velocities are stored
 **********************************************************************************************************************************/
__device__ void Velocity_HE(double* vel, double* U)
{
  for (int n = 0; n < 3; n++)
  {
    vel[n] = U[EXT_MOM1-1+n]*U[EXT_SRHO-1];
  }
}

/***********************************************************************************************************************************
 * @brief Device function rendering of the PRESSURE_HE macro from eos.h
 **********************************************************************************************************************************/
__device__ double Pressure_HE(double* U)
{
  return d_KappaM1*U[EXT_ENER-1] - 0.5*(U[EXT_VEL1-1]*U[EXT_MOM1-1] + U[EXT_VEL2-1]*U[EXT_MOM2-1] + U[EXT_VEL3-1]*U[EXT_MOM3-1]);
}

/***********************************************************************************************************************************
 * @brief Device function rendering of the SPEEDOFSOUND_HE macro from eos.h
 **********************************************************************************************************************************/
__device__ double SpeedOfSound_HE(double* U)
{
  return sqrt(d_Kappa*U[EXT_PRES-1]*U[EXT_SRHO-1]);
}

/***********************************************************************************************************************************
 * @brief Device function rendering of the TOTALENERGY_HE macro from eos.h
 **********************************************************************************************************************************/
__device__ double TotalEnergy_HE(double* U)
{
  return U[EXT_ENER-1]*U[EXT_SRHO-1];
}

/***********************************************************************************************************************************
 * @brief Device function rendering of the TOTALENTHALPY_HE macro from eos.h
 **********************************************************************************************************************************/
__device__ double TotalEnthalpy_HE(double* U)
{
  return (U[EXT_ENER-1]+U[EXT_PRES-1])*U[EXT_SRHO-1];
}

/***********************************************************************************************************************************
 * @brief Device function rendering of the TEMPERATURE_HE macro from eos.h
 **********************************************************************************************************************************/
__device__ double Temperature_HE(double* U)
{
  return U[EXT_PRES-1]*U[EXT_SRHO-1]/d_R;
}

/***********************************************************************************************************************************
 * @brief Device function rendering of the ENERGY_HE macro from eos.h
 **********************************************************************************************************************************/
__device__ double Energy_HE(double* U)
{
  return d_sKappaM1*U[EXT_PRES-1] + 0.5*(U[EXT_MOM1-1]*U[EXT_VEL1-1] + U[EXT_MOM2-1]*U[EXT_VEL2-1] + U[EXT_MOM3-1]*U[EXT_VEL3-1]);
}

#if PARABOLIC
#if PP_VISC == 1
/***********************************************************************************************************************************
 * @brief Sutherland's formula can be used to derive the dynamic viscosity of an ideal gas as a function of the temperature
 * @remark Initialization of d_mu0, d_Ts, d_Tref, d_ExpoSuth and cSuth takes place in InitEOS_Device
 * @remark The host copy of this method appears in viscosity.f90. It is put here to prevent creation of separate source file
 *         just for it. If it becomes necessary or logical to move it to a "viscosity.cu" file, do so.
 **********************************************************************************************************************************/
__device__ double MuSuth(double T)
{
  double TnoDim = T*d_Tref;
  // Attention: only valid for T < 550K. But we don't know what to do for higher temperatures...
  if (TnoDim >= d_Ts)
  {
    return pow(d_mu0*TnoDim, d_ExpoSuth*(1.0+d_Ts)/(TnoDim+d_Ts));
  }
  else
  {
    return d_mu0*TnoDim*d_cSuth;
  }
}
#endif /* PP_VISC == 1 */

/***********************************************************************************************************************************
 * @brief Device function rendering of the VISCOSITY_PRIM macro from eos.h
 **********************************************************************************************************************************/
__device__ double ViscosityPrim(double* U)
{
#if PP_VISC == 0
  return d_mu0;
#elif PP_VISC == 1
  return muSuth(U[TEMP-1]);
#elif PP_VISC == 2
  return pow(d_mu0*U[TEMP-1], d_ExpoSuth);
#endif
}
#endif /* PARABOLIC */
