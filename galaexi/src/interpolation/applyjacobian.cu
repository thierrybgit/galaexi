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
    void ApplyJacobian_Device(int TP_nVar, int Nloc, int nElems, int d_U, bool toPhysical, int d_sJ, int FVE
#if FV_ENABLED
                                   ,int d_FV_Elems, bool chooseFV
#endif                                  
                                   ,int streamID
                                  );

}

//----------------------------------------------------------------------------------------------------------------------------------
// Kernel
//----------------------------------------------------------------------------------------------------------------------------------
__global__ void ApplyJacobian_Kernel(int TP_nVar,int Nloc,int nDOFs,double* U,bool toPhysical,double* sJ, int FVE
#if FV_ENABLED
                                     ,int* FV_Elems, bool chooseFV
#endif
                                    )
{
    int FVidx = 0; // If FV is off, then we will always use FVidx = 0, so set that as the default
    int U_offset = 0;
    int J_offset = 0;
#if FV_ENABLED
    int elemIdx = 0;
    int nDOFsElem = (Nloc+1)*(Nloc+1)*(Nloc+1);
#endif

    // We're working on the entire volume, so just give a single DOF to a thread
    int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    if (threadID < nDOFs)
    {
        // Compute the offsets to this thread's DOF in the flattened arrays
        U_offset = threadID*TP_nVar;

#if FV_ENABLED
        // If running with FV, we'll need to figure out which element the DOF assigned to the current device thread
        elemIdx = floor((double)i / (double)nDOFsElem);
        if (chooseFV)
        {
            if (FV_Elems[elemIdx] == FVE)
            {
                FVidx = FV_Elems[elemIdx];
            }
            else
            {
                return; // We are choosing which elements to work on and this is not one of them, exit
            }               
        }
        else // FV is on, but we aren't choosing, so use the status of this element no matter what
        {
            FVidx = FV_Elems[elemIdx];
        }
#endif

        // Now that we know the element status for FV, we can offset the Jacobian array to acccount for it
        // All this does is chooses the starting point in the FV_SIZE dim and then set this thread to proper
        // DOF from that starting point.
        // sJ(0:Nloc,0:Nloc,0:Nloc,nElems,0:FV_SIZE)
        J_offset = FVidx*nDOFs + threadID;

        if (toPhysical)
        {
            for (int n = 0; n < TP_nVar; n++)
            {
                U[U_offset+n] *= sJ[J_offset];
            }
        }
        else // ! toPhysical
        {
            for (int n = 0; n < TP_nVar; n++)
            {
                U[U_offset+n] /= sJ[J_offset];
            }
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------
// Interface - Call kernels from HOST
//----------------------------------------------------------------------------------------------------------------------------------
void ApplyJacobian_Device(int TP_nVar, int Nloc, int nElems, int d_U, bool toPhysical, int d_sJ, int FVE
#if FV_ENABLED
                           ,int d_FV_Elems, bool chooseFV
#endif                                  
                           ,int streamID
                         )
{
#ifdef USE_NVTX /* If profiling, add range */
    nvtxRangePushA("ApplyJacobian");
#endif

    int nDOFs = (Nloc+1)*(Nloc+1)*(Nloc+1)*nElems;
    INVOKE_KERNEL(ApplyJacobian_Kernel, nDOFs/256+1, 256, 0, streams[streamID],
                                     TP_nVar, Nloc, nDOFs, (double*)DeviceVars[d_U], 
                                     toPhysical, (double*)DeviceVars[d_sJ], FVE
#if FV_ENABLED
                                     ,(int*)DeviceVars[d_FV_Elems], chooseFV
#endif
                 );

#ifdef USE_NVTX
    nvtxRangePop();
#endif
}