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
    void FillFlux_Device(int Nloc, int nSides, double t, int firstInnerSide, int lastInnerSide, int firstMPISide_MINE, int lastMPISide_MINE,
                            int firstBCSide, int nBCSides, bool doMPISides,
                            int d_Flux_master, int d_Flux_slave, int d_U_master, int d_U_slave, int d_UPrim_master, int d_UPrim_slave,
                            int d_NormVec, int d_TangVec1, int d_TangVec2, int d_SurfElem, int d_Face_xGP
                            ,int streamID
                        );
}


//----------------------------------------------------------------------------------------------------------------------------------
// Kernels
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Step 1: Compute flux for non-BC sides
 **********************************************************************************************************************************/
__global__ void __launch_bounds__(256, 1) FillFlux_Fluxes_Kernel(int nDOFs, int Nloc, int firstSideID_wo_BC, int woBCFlux_offset, int woBCFluxPrim_offset, int woBCVec_offset,
                                    double* Flux_master, double* U_master, double* U_slave, double* UPrim_master, double* UPrim_slave, 
                                    double* NormVec, double* TangVec1, double* TangVec2
                                )
{

    int locFlux_offset = 0;
    int locFluxPrim_offset = 0;
    int locVec_offset = 0;

    // This device thread takes a single DOF somewhere among the sides we are working on
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nDOFs) return;

    //               start elem offset + offset to the DOF for this thread
    locFlux_offset = woBCFlux_offset + i*PP_nVar;
    locFluxPrim_offset = woBCFluxPrim_offset + i*PP_nVarPrim;
    locVec_offset = woBCVec_offset + i*3;

    // 1.1) Advective flux from Riemann solver
    Riemann_Point_Device(&Flux_master[locFlux_offset], &U_master[locFlux_offset], &U_slave[locFlux_offset], 
                            &UPrim_master[locFluxPrim_offset], &UPrim_slave[locFluxPrim_offset], 
                            &NormVec[locVec_offset], &TangVec1[locVec_offset], &TangVec2[locVec_offset], false);

#if PARABOLIC
    // 1.2) Viscous flux
#endif /*PARABOLIC*/

    // 1.3) add up viscous flux
    // THIS WILL BE RE-ADDED AFTER VISCOUS FLUX IS PORTED AND CALLED HERE
    // Flux_master[locFlux_offset] += FluxV_loc

}

/***********************************************************************************************************************************
 * @brief Step 2: Compute the fluxes at the boundaries
 **********************************************************************************************************************************/
__global__ void FillFlux_BCs_Kernel()
{
    // This device thread takes a single DOF somewhere among the sides we are working on
    // int i = blockIdx.x * blockDim.x + threadIdx.x;

    // if (i < nDOFs)
    // {

    // }
}

/***********************************************************************************************************************************
 * @brief Step 3/4: Final flux computation and copy to slave sides
 **********************************************************************************************************************************/
__global__ void __launch_bounds__(256, 1) FillFlux_PopulateFluxes_Kernel(int nDOFs, int Nloc, int firstSideID,
                                    int firstSideFlux_offset, int firstSideSurfElem_offset, 
                                    double* Flux_master, double* Flux_slave, double* SurfElem)
{

    int locFlux_offset = 0;
    int locSurfElem_offset = 0;

    // This device thread takes a single DOF somewhere among the sides we are working on
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nDOFs) return;

    //               start elem offset + offset to the DOF for this thread
    locFlux_offset = firstSideFlux_offset + i*PP_nVar;
    locSurfElem_offset = firstSideSurfElem_offset + i;

    for(int n = 0; n < PP_nVar; n++)
    {
        Flux_master[locFlux_offset+n] *= SurfElem[locSurfElem_offset];
        Flux_slave[locFlux_offset+n] = Flux_master[locFlux_offset+n]; 
    }
}

#if FV_ENABLED == 1
/***********************************************************************************************************************************
 * @brief Step 5: For FV switching, convert flux on FV points to DG points for all DG faces at mixed interfaces
 **********************************************************************************************************************************/
__global__ void FillFlux_ConvertFVFluxes_Kernel()
{
    // This device thread takes a single DOF somewhere among the sides we are working on
    // int i = blockIdx.x * blockDim.x + threadIdx.x;

    // if (i < nDOFs)
    // {

    // }
}
#endif

//----------------------------------------------------------------------------------------------------------------------------------
// Interface to Host (Fortran)
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Entry point to launch kernels for device side computation of FillFlux.
 * @remark THIS METHOD RUNS ON THE HOST!!!
 **********************************************************************************************************************************/
void FillFlux_Device(int Nloc, int nSides, double t, int firstInnerSide, int lastInnerSide, int firstMPISide_MINE, int lastMPISide_MINE,
                        int firstBCSide, int nBCSides, bool doMPISides,
                        int d_Flux_master, int d_Flux_slave, int d_U_master, int d_U_slave, int d_UPrim_master, int d_UPrim_slave,
                        int d_NormVec, int d_TangVec1, int d_TangVec2, int d_SurfElem, int d_Face_xGP
                        ,int streamID
                    )
{
    int firstSideID_wo_BC;
    int firstSideID, lastSideID;
    int nDOFs;
    int woBCFlux_offset = 0;
    int woBCVec_offset = 0;
    int woBCFluxPrim_offset = 0;
    int firstSideFlux_offset = 0;
    int firstSideSurfElem_offset = 0;

    // Fill flux for sides ranging between firstSideID and lastSideID using Riemann solver for advection and viscous terms
    // Set the side range according to MPI or no MPI
    // WE ADJUST LOCAL SIDEIDs TO BEING INDEXED FROM 0 HERE FOR EASE CALCULATING FLATTENED MEM OFFSETS LATER
    if (doMPISides)
    {
        // fill only flux for MINE MPISides (where the local proc is master)
        firstSideID_wo_BC = firstMPISide_MINE - 1;
        firstSideID = firstMPISide_MINE - 1;
        lastSideID =  lastMPISide_MINE - 1;
    }
    else
    {
        // fill only InnerSides that do not need communication
        firstSideID_wo_BC = firstInnerSide - 1; // for fluxes
        firstSideID = firstBCSide - 1;    // include BCs for master sides
        lastSideID = lastInnerSide - 1;
    }

    // We CANNOT do this (allocate a host array and then just pass it to a kernel)
    // For now will just use zero in the kernels
    // FV_Elems_Max = (int*)malloc(nSides*sizeof(int));
    // for (int n = 0; n < nSides; n++)
    // {
    //     FV_Elems_Max[n] = 0;
    // }
#if FV_ENABLED
#endif

    /* =============================
    * Workflow:
    *
    *  1.  compute flux for non-BC sides
    *  1.1) advective flux
    *  1.2) viscous flux
    *  1.3) add up viscous flux to Flux_master
    *  2.  compute flux for BC sides
    *  3.  multiply by SurfElem
    *  4.  copy flux from Flux_master to Flux_slave
    *  5.  convert FV flux to DG flux at mixed interfaces
    *============================== */

    // 2.) Compute the fluxes at the boundary conditions: 1..nBCSides
    // ATTENTION: This kernel is typically very small.
    //            Compute first when there is still volume work to fill the GPU
    // Here we use the BC side offsets, which just start from side 0 in the arrays
    // THIS WILL BE ADDED AFTER BoundaryFlux IS PORTED
    // if (firstSideID < firstSideID_wo_BC)
    // {
    //     nDOFs = (Nloc+1)*(Nloc+1)*(nBCSides);
// #ifdef USE_NVTX /* If profiling, add range */
//          nvtxRangePushA("FillFlux_BCs");
// #endif
    //     INVOKE_KERNEL( FillFlux_BCs_Kernel, nDOFs/256+1, 256, streams[streamID], ...);
// #ifdef USE_NVTX
//          nvtxRangePop();
// #endif
    // }


    if (firstSideID_wo_BC <= lastSideID)
    {
        // 1.) Compute flux for non-BC sides
        // ------------------------------------------------
        woBCFlux_offset = FindSideMasterSlaveOffset(PP_nVar, Nloc, firstSideID_wo_BC); // NOTE: All offset methods assume 0-indexed sideID
        woBCFluxPrim_offset = FindSideMasterSlaveOffset(PP_nVarPrim, Nloc, firstSideID_wo_BC);
        woBCVec_offset = FindSideSurfaceDataOffset(3, Nloc, nSides, 0, firstSideID_wo_BC); // NOTE: Need to sub in FV_Elem_Max here
        nDOFs = (Nloc+1)*(Nloc+1)*(lastSideID-firstSideID_wo_BC+1);

#ifdef USE_NVTX /* If profiling, add range */
        nvtxRangePushA("FillFlux_Fluxes");
#endif
        INVOKE_KERNEL( FillFlux_Fluxes_Kernel, nDOFs/256+1, 256, 0, streams[streamID], nDOFs, Nloc, firstSideID_wo_BC, woBCFlux_offset, woBCFluxPrim_offset, woBCVec_offset,
                                                                                (double*)DeviceVars[d_Flux_master], (double*)DeviceVars[d_U_master], 
                                                                                (double*)DeviceVars[d_U_slave], (double*)DeviceVars[d_UPrim_master], (double*)DeviceVars[d_UPrim_slave], 
                                                                                (double*)DeviceVars[d_NormVec], (double*)DeviceVars[d_TangVec1], (double*)DeviceVars[d_TangVec2]
                    );
#ifdef USE_NVTX
        nvtxRangePop();
#endif
    }

    if (firstSideID <= lastSideID)
    {
        // 3.) & 4.) Final flux computation and copy to slave sides
        firstSideFlux_offset = FindSideMasterSlaveOffset(PP_nVar, Nloc, firstSideID);
        firstSideSurfElem_offset = FindSideSurfaceDataOffset(1, Nloc, nSides, 0, firstSideID); // NOTE: Need to sub in FV_Elem_Max here
        nDOFs = (Nloc+1)*(Nloc+1)*(lastSideID-firstSideID+1);
#ifdef USE_NVTX /* If profiling, add range */
        nvtxRangePushA("FillFlux_PopulateFluxes");
#endif
        INVOKE_KERNEL( FillFlux_PopulateFluxes_Kernel, nDOFs/256+1, 256, 0 , streams[streamID], nDOFs, Nloc, firstSideID,
                                        firstSideFlux_offset, firstSideSurfElem_offset, 
                                        (double*)DeviceVars[d_Flux_master], (double*)DeviceVars[d_Flux_slave], 
                                        (double*)DeviceVars[d_SurfElem] 
                 );
#ifdef USE_NVTX
        nvtxRangePop();
#endif
    }

#if FV_ENABLED == 1
    // 5. convert flux on FV points to DG points for all DG faces at mixed interfaces
    // only inner sides can be mixed (BC do not require a change basis)
    // nDOFs = (Nloc+1)*(Nloc+1)*(lastSideID-firstSideID_wo_BC+1);
// #ifdef USE_NVTX /* If profiling, add range */
//      nvtxRangePushA("FillFlux_ConvertFVFluxes");
// #endif
    // INVOKE_KERNEL(FillFlux_ConvertFVFluxes_Kernel, nDOFs/256+1, 256, 0, streams[streamID], ... );
// #ifdef USE_NVTX
//      nvtxRangePop();
// #endif
#endif

    // Make sure to free the allocated mem
    // free(FV_Elems_Max);
}

