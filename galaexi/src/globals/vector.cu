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
    void VAXPB_2STEP_Device(int nTotal,int d_VecOut,int d_VecOut2,int d_VecIn,double Const1,double Const2);
    void VAXPB_IN_OUT_Device(int nTotal,int d_VecOut,int d_VecIn,double ConstOut,double ConstIn);
    void VAXPB_OUT_Device(int nTotal,int d_VecOut,int d_VecIn,double Const);
    void VAXPB_IN_Device(int nTotal,int d_VecOut,int d_VecIn,double Const);
    void VAXPB_ADD_Device(int nTotal,int d_VecOut,int d_VecIn);
    void VAXPB_CONST_Device(int nTotal,int d_VecOut,double Const);
}


//----------------------------------------------------------------------------------------------------------------------------------
// Kernels
//----------------------------------------------------------------------------------------------------------------------------------
__global__ void VAXPB_2STEP_Kernel(int nTotal,double* VecOut,double* VecOut2,double* VecIn,double Const1,double Const2)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nTotal) return;
    
    VecOut[i]  = VecOut[i]*Const1 + VecIn[i];
    VecOut2[i] = VecOut2[i] + VecOut[i]*Const2;
}

__global__ void VAXPB_IN_OUT_Kernel(int nTotal,double* VecOut,double* VecIn,double ConstOut,double ConstIn)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nTotal) return;

    VecOut[i] = VecOut[i]*ConstOut + VecIn[i]*ConstIn;
}

__global__ void VAXPB_OUT_Kernel(int nTotal,double* VecOut,double* VecIn,double Const)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nTotal) return;

    VecOut[i] = VecOut[i]*Const + VecIn[i];
}

__global__ void VAXPB_IN_Kernel(int nTotal,double* VecOut,double* VecIn,double Const)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nTotal) return;

    VecOut[i] = VecOut[i] + VecIn[i]*Const;
}

__global__ void VAXPB_ADD_Kernel(int nTotal,double* VecOut,double* VecIn)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nTotal) return;

    VecOut[i]=VecOut[i]+VecIn[i];
}

__global__ void VAXPB_CONST_Kernel(int nTotal,double* VecOut,double Const)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nTotal) return;

    VecOut[i]=VecOut[i]*Const;
}


//----------------------------------------------------------------------------------------------------------------------------------
// Interfaces - Call kernels from HOST
// See Fortran methods in vector.f90 for methods documentation
//----------------------------------------------------------------------------------------------------------------------------------
void VAXPB_2STEP_Device(int nTotal,int d_VecOut,int d_VecOut2,int d_VecIn,double Const1,double Const2)
{
#ifdef USE_NVTX /* If profiling, add range */
    nvtxRangePushA("VAXPB_2STEP");
#endif

    INVOKE_KERNEL(VAXPB_2STEP_Kernel, nTotal/256+1, 256, 0, streams[0], nTotal, (double*)DeviceVars[d_VecOut], 
                                        (double*)DeviceVars[d_VecOut2], (double*)DeviceVars[d_VecIn], Const1, Const2);

#ifdef USE_NVTX
    nvtxRangePop();
#endif
}

void VAXPB_IN_OUT_Device(int nTotal,int d_VecOut,int d_VecIn,double ConstOut,double ConstIn)
{
#ifdef USE_NVTX /* If profiling, add range */
    nvtxRangePushA("VAXPB_IN_OUT");
#endif

    INVOKE_KERNEL(VAXPB_IN_OUT_Kernel, nTotal/256+1, 256, 0, streams[0], nTotal, (double*)DeviceVars[d_VecOut], 
                                        (double*)DeviceVars[d_VecIn], ConstOut, ConstIn);

#ifdef USE_NVTX
    nvtxRangePop();
#endif
}

void VAXPB_OUT_Device(int nTotal,int d_VecOut,int d_VecIn,double Const)
{
#ifdef USE_NVTX /* If profiling, add range */
    nvtxRangePushA("VAXPB_OUT");
#endif

    INVOKE_KERNEL(VAXPB_OUT_Kernel, nTotal/256+1, 256, 0, streams[0], nTotal, (double*)DeviceVars[d_VecOut], 
                                        (double*)DeviceVars[d_VecIn], Const);

#ifdef USE_NVTX
    nvtxRangePop();
#endif
}

void VAXPB_IN_Device(int nTotal,int d_VecOut,int d_VecIn,double Const)
{
#ifdef USE_NVTX /* If profiling, add range */
    nvtxRangePushA("VAXPB_IN");
#endif

    INVOKE_KERNEL(VAXPB_IN_Kernel, nTotal/256+1, 256, 0, streams[0], nTotal, (double*)DeviceVars[d_VecOut], 
                                        (double*)DeviceVars[d_VecIn], Const);

#ifdef USE_NVTX
    nvtxRangePop();
#endif
}

void VAXPB_ADD_Device(int nTotal,int d_VecOut,int d_VecIn)
{
#ifdef USE_NVTX /* If profiling, add range */
    nvtxRangePushA("VAXPB_ADD");
#endif

    INVOKE_KERNEL(VAXPB_ADD_Kernel, nTotal/256+1, 256, 0, streams[0], nTotal, (double*)DeviceVars[d_VecOut], 
                                        (double*)DeviceVars[d_VecIn]);

#ifdef USE_NVTX
    nvtxRangePop();
#endif
}

void VAXPB_CONST_Device(int nTotal,int d_VecOut,double Const)
{
#ifdef USE_NVTX /* If profiling, add range */
    nvtxRangePushA("VAXPB_CONST");
#endif

    INVOKE_KERNEL(VAXPB_CONST_Kernel, nTotal/256+1, 256, 0, streams[0], nTotal, (double*)DeviceVars[d_VecOut], Const);

#ifdef USE_NVTX
    nvtxRangePop();
#endif
}