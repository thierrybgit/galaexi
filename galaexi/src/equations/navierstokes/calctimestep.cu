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
void CalcTimeStep_Device(int Nloc,int nElems,int d_errType,int d_TSconv,int d_TSvisc,int d_U,
                          int d_sJ,int d_Mf,int d_Mg,int d_Mh, 
                          int d_CFLScale,int d_MetAdv,double Kappa,double KappaM1,double R
#if PARABOLIC
                          ,int d_MetVisc,int d_DFLScale
#endif /*PARABOLIC*/
#if FV_ENABLED
                          ,int d_FV_Elems
#if FV_ENABLED == 2
                          ,int d_FV_alpha,double FV_alpha_min
#endif /* FV_ENABLED == 2 */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
                          ,int d_muSGS
#endif /* EDDYVISCOSITY */
                        );
}


//----------------------------------------------------------------------------------------------------------------------------------
// Helper Methods
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Simple method to find the L_2 norm of an array.
 * @param arr Array of data
 * @param arrSize Number of elements in the data array
 * @returns L_2 norm of the passed array
 **********************************************************************************************************************************/
__device__ double norm2(double* arr, int arrSize)
{
    double sum = 0.0;

    for (int i = 0; i < arrSize; i++)
    {
        sum += arr[i]*arr[i];
    }

    return sqrt(sum);
}

#if PARABOLIC
/***********************************************************************************************************************************
 * @brief Perform the array summation needed in the viscous terms of InitCalcTimeStep.
 * Fortran equivalent -- SUM((Mf(:,i,j,k,iElem,FVE)*sJ(i,j,k,iElem,FVE))**2)
 * @param M Metrics for the DOF assigned to calling device thread
 * @param J Value of the Jacobian for the DOF assigned to calling device thread
 * @returns Computed sum ()
 **********************************************************************************************************************************/
__device__ double MetricsJacobianSum(double* M, double J)
{
    double sum = 0.0; 

    for (int i = 0; i < 3; i++)
    {
        sum += (M[i]*J)*(M[i]*J);
    }

    return sum;
}
#endif /* PARABOLIC */

/***********************************************************************************************************************************
 * @brief CUDA does not support atomicMin with doubles. This is a simple workaround using atomicCAS.
 **********************************************************************************************************************************/
__device__ double atomicMin_double(double* address, double val)
{
#if (USE_ACCEL == ACCEL_CUDA)
    unsigned long long int* address_as_ull = (unsigned long long int*) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
#elif (USE_ACCEL == ACCEL_HIP)
    // HIP contains an actual double* overload for atomicMin.
    // So instead of the hacky workaround used for CUDA above, just use the HIP lib function
    double old = atomicMin(address, val);
    return old;
#endif
}

//----------------------------------------------------------------------------------------------------------------------------------
// CalcTimeStep Methods
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Compute the max eigenvalue in a given coordinate direction for a single DOF
 * @param M Metrics at the current DOF
 * @param curr_val Current eigenvalue value at this DOF for this direction
 * @param MetAdv Advective metrics at this DOF for this direction
 * @param VsJ Component of physical velocity for this direction at this DOF
 * @param c Speed of Sound at this DOF
 * @returns Maximum eigenvalue at this DOF for a single coordinate direction
 **********************************************************************************************************************************/
__device__ double FindMaxEigenValue(double* M, double* curr_val, double* MetAdv, double* VsJ, double* c)
{
    double metSum = 0.0;
    for (int n = 0; n < 3; n++)
    {
        metSum += M[n];
    }
    metSum = fabs(metSum*VsJ[0]);
    return fmax(curr_val[0], metSum + c[0]*MetAdv[0]);
}

/***********************************************************************************************************************************
 * @brief Find maximum value for SGS viscosity for this element.
 * @param muSGS Pointer to start of element assigned to this device thread in muSGS array
 * muSGS_Max = MAXVAL(muSGS(1,:,:,:,iElem))
 **********************************************************************************************************************************/
__device__ double FindMaxMuSGS(double* muSGS, int nDOFsThisElem)
{
    double muSGS_Max = muSGS[0];

    for (int i = 1; i < nDOFsThisElem; i++)
    {
        if ( muSGS[i] > muSGS_Max ) {
            muSGS_Max = muSGS[i];
        }
    }

    return muSGS_Max;
}

#if FV_ENABLED == 3
/***********************************************************************************************************************************
 * @brief Determine if all FV_allpha values for an element are less than FV_alpha_min
 * @param FV_Alpha Pointer to start of the element in FV_alpha array
 * @param FV_alpha_min Threshold value to check FV_alpha alues against
 * @param nDOFsThisElem Number of DOFs in the element
 * @returns True if all FV_alpha values in this element are less than FV_alpha_min. Otherwise False.
 **********************************************************************************************************************************/
__device__ bool CheckFV_alpha(double* FV_alpha, double FV_alpha_min, int nDOFsThisElem)
{
    bool result = true;
    // Search DOFs in the element until a value greater than FV_alpha_min is found
    for (int i = 0; i < nDOFsThisElem; n++)
    {
        if (FV_alpha[i] > FV_alpha_min)
        {
            result = false;
            break;
        }
    }
    return result;
}
#endif

/***********************************************************************************************************************************
 * @brief Kernel for CalcTimeStep. Kernel is element-parallel.
 **********************************************************************************************************************************/
__global__ void __launch_bounds__(256, 1) CalcTimeStep_Kernel(int Nloc,int nElems,int* errType, double* TSconv,double* TSvisc,double* U,double* sJ,
                                        double* Mf,double* Mg,double* Mh,double* CFLScale,double* MetAdv, 
                                        double Kappa,double KappaM1,double R
#if PARABOLIC
                                        ,double* MetVisc,double* DFLScale
#endif /*PARABOLIC*/
#if FV_ENABLED
                                        ,int* FV_Elems
#if FV_ENABLED == 2
                                        ,double* FV_alpha,double FV_alpha_min
#endif /* FV_ENABLED == 2 */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
                                        ,double* muSGS
#endif /* EDDYVISCOSITY */
                                    )
{
    int U_offset = 0;
    int J_offset = 0;
    int M_offset = 0;
    int FVE;
    double Ue[PP_2Var];
    double c;
    double VsJ[3];
    int nDOFsThisElem = (Nloc+1)*(Nloc+1)*(Nloc+1);
    double MaxLambda[3] = {0.0, 0.0, 0.0};
    double SumMaxLambda = 0.0;
    int thisThreadErrType = 0;
#if PARABOLIC
    double mu;
    double MaxLambdaVisc[3] = {0.0, 0.0, 0.0};
    double SumMaxLambdaVisc = 0.0;
#endif
#if EDDYVISCOSITY
    int SGS_offset = 0;
    double muSGS_Max;
#endif
#if FV_ENABLED == 2
    int FV_alpha_offset = 0;
    double minCFL_Scale = 0.0;
#endif

    int threadID = blockIdx.x*blockDim.x + threadIdx.x;
    // For this kernel, every thread gets an element. This assignment is just for readability
    // +1 because all offset methods assumed that ElemID is 1-indexed (because it is on Fortran side)
    int ElemID = threadID + 1;
    *TSconv = DBL_MAX;
    *TSvisc = DBL_MAX;
    if (threadID == 0) *errType = 0;

    if (ElemID <= nElems)
    {
#if FV_ENABLED
        // Have to -1 to adjust for 1-indexing imposed above.
        FVE = FV_Elems[ElemID-1];
#else
        FVE = 0;
#endif
#if EDDYVISCOSITY
        SGS_offset = FindElemOffset(1, Nloc, ElemID);
        muSGS_Max = FindMaxMuSGS(muSGS[SGS_offset], nDOFsThisElem);
#endif
        // Galaexi v1.0 had this process as a separate, DOF-parallel kernel
        // I (Spencer Starr - 08.2024) am betting that, since each element's DOFS are contiguous
        // in memory, keeping this and the below min process as a single element-parallel kernel 
        // will be overall faster than launching a separate kernel for each process.
        // If that turns out to NOT be the case, then Galaexi v1.0 holds the pattern for how to
        // rewrite the CalcTimeStep kernel in the other way.
        for (int i = 0; i < nDOFsThisElem; i++)
        {
            U_offset = FindElemOffset(PP_nVar, Nloc, ElemID) + i*PP_nVar;
            for (int n = 0; n < PP_nVar; n++)
            {
                Ue[n] = U[U_offset+n];
            }
            Ue[EXT_SRHO-1] = 1.0/Ue[EXT_DENS-1];
            Velocity_HE(&Ue[EXT_VEL1-1],&Ue[0]);
            Ue[EXT_PRES-1] = Pressure_HE(&Ue[0]);
            Ue[EXT_TEMP-1] = Temperature_HE(&Ue[0]);

            if (isnan(U[U_offset]))
            {
                thisThreadErrType = 1;
            }
            atomicAdd(errType, thisThreadErrType);

            J_offset = FindElemJacobianOffset(Nloc, ElemID, FVE) + i;
            M_offset = FindElemMetricsOffset(Nloc, ElemID, FVE) + i*3;
            for (int n = 0; n < 3; n++)
            {
                VsJ[n] = Ue[EXT_VEL1-1+n]*sJ[J_offset];
            }
            c = SpeedOfSound_HE(&Ue[0]);

            MaxLambda[0] = FindMaxEigenValue(&Mf[M_offset],&MaxLambda[0], &MetAdv[M_offset+0],&VsJ[0],&c);
            MaxLambda[1] = FindMaxEigenValue(&Mg[M_offset],&MaxLambda[1], &MetAdv[M_offset+1],&VsJ[1],&c);
            MaxLambda[2] = FindMaxEigenValue(&Mh[M_offset],&MaxLambda[2], &MetAdv[M_offset+2],&VsJ[2],&c);
#if PARABOLIC
            // Viscous Eigenvalues
            mu = ViscosityPrim(Ue[EXT_PRIM-1]);
#if EDDYVISCOSITY
            mu += muSGS_Max;
#endif
            for (int n = 0; n < 3; n++)
            {
                MaxLambdaVisc[n] = fmax(MaxLambdaVisc[n], Ue[EXT_SRHO-1]*MetVisc[M_offset+n]);
            }
#endif /* PARABOLIC*/
        }

        SumMaxLambda = MaxLambda[0] + MaxLambda[1] + MaxLambda[2];
#if FV_ENABLED == 2
        // In the Fortran here, FV_alpha (a 5D array) is indexed as a 1D array
        // I (Spencer Starr 08.2024) have no idea what the hell that is about or how its actually working
        // or how who knows how many compilers have let that slide.
        // For now, just use the first value for this thread's element and move on
        FV_alpha_offset = FindElemOffset(1,Nloc,ElemID)
        if (FV_alpha[FV_alpha_offset] <= FV_alpha_min)
        {
            atomicMin_double(TSconv,CFLScale[0]*2.0/SumMaxLambda);
        }
        else
        {
            for (int n = 0; n < FV_SIZE; n++)
            {
                if (minCFL_Scale > CFLScale[n])
                {
                    minCFL_Scale = CFLScale[n];
                }
            }
            atomicMin_double(TSconv, minCFL_Scale*2.0/SumMaxLambda);
        }
#elif FV_ENABLED == 3
        // Unlike above, this is much more explicit on what is expected. The if condition must be true
        // for ALL the DOFs in the element.
        FV_alpha_offset = FindElemOffset(1,Nloc,ElemID)
        if (CheckFV_alphaMin(&FV_alpha[FV_alpha_offset],FV_alpha_min,nDOFsThisElem))
        {
            atomicMin_double(TSconv,CFLScale[0]*2.0/SumMaxLambda);
        }
        else
        {
            for (int n = 0; n < FV_SIZE; n++)
            {
                if (minCFL_Scale > CFLScale[n])
                {
                    minCFL_Scale = CFLScale[n];
                }
            }
            atomicMin_double(&TSconv, minCFL_Scale*2.0/SumMaxLambda);
        }
#else /* FV_ENABLED */
        atomicMin_double(TSconv, CFLScale[FVE]*2.0/SumMaxLambda);
#endif /* FV_ENABLED */

#if PARABOLIC
        SumMaxLambdaVisc = MaxLambdaVisc[0] + MaxLambdaVisc[1] + MaxLambdaVisc[2];
        if (SumMaxLambdaVisc > 0.0)
        {
#if FV_ENABLED == 2 || FV_ENABLED == 3
            for (int n = 0; n < FV_SIZE; n++)
            {
                if (minCFL_Scale > DFLScale[n])
                {
                    minCFL_Scale = DFLScale[n];
                }
            }
            atomicMin_double(TSvisc, minCFL_Scale*4.0/SumMaxLambdaVisc);
#else
            atomicMin_double(TSvisc, DFLScale[FVE]*4.0/SumMaxLambdaVisc);
#endif
        }
#endif /* PARABOLIC*/

    }
}

/***********************************************************************************************************************************
 * @brief Entry point for the CalcTimeStep device backend. Runs on host and calls device kernel.
 **********************************************************************************************************************************/
void CalcTimeStep_Device(int Nloc,int nElems,int d_errType,int d_TSconv,int d_TSvisc,int d_U,
                          int d_sJ,int d_Mf,int d_Mg,int d_Mh, 
                          int d_CFLScale,int d_MetAdv,double Kappa,double KappaM1,double R
#if PARABOLIC
                          ,int d_MetVisc,int d_DFLScale
#endif /*PARABOLIC*/
#if FV_ENABLED
                          ,int d_FV_Elems
#if FV_ENABLED == 2
                          ,int d_FV_alpha,double FV_alpha_min
#endif /* FV_ENABLED == 2 */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
                          ,int d_muSGS
#endif /* EDDYVISCOSITY */
                        )
{
#ifdef USE_NVTX
  nvtxRangePushA("CalcTimeStep");
#endif

    // Call the kernel
    INVOKE_KERNEL( CalcTimeStep_Kernel, nElems/256+1, 256, 0, streams[0], Nloc,nElems, (int*)DeviceVars[d_errType],
                          (double*)DeviceVars[d_TSconv],(double*)DeviceVars[d_TSvisc],
                          (double*)DeviceVars[d_U],(double*)DeviceVars[d_sJ],
                          (double*)DeviceVars[d_Mf],(double*)DeviceVars[d_Mg],(double*)DeviceVars[d_Mh],
                          (double*)DeviceVars[d_CFLScale],(double*)DeviceVars[d_MetAdv],
                          Kappa,KappaM1,R
#if PARABOLIC
                          ,(double*)DeviceVars[d_MetVisc],(double*)DeviceVars[d_DFLScale]
#endif /*PARABOLIC*/
#if FV_ENABLED
                          ,(int*)DeviceVars[d_FV_Elems]
#if FV_ENABLED == 2
                          ,(double*)DeviceVars[d_FV_alpha],double FV_alpha_min
#endif /* FV_ENABLED == 2 */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
                          ,(double*)DeviceVars[d_muSGS]
#endif /* EDDYVISCOSITY */
                 );

#ifdef USE_NVTX
    nvtxRangePop();
#endif
}
