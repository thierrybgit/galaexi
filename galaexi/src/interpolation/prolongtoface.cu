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
    void ProlongToFace_Device(int TP_nVar, int Nloc, int nElems, int nSides, int firstMPISide_YOUR, int lastMPISide_MINE,
                                int d_Uvol, int d_Uface_master, int d_Uface_slave, int d_L_Minus, int d_L_Plus, 
                                int d_SideToElem, int d_S2V2
#if FV_ENABLED
                                ,int d_isFV
#endif
                                , int streamID
                             );

}



//----------------------------------------------------------------------------------------------------------------------------------
// Kernels
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Interpolates the element volume data stored at Gauss points
 * @returns Value for Uface at U(n,:,:,:,ElemID)
 * @param n Index of flow variable being worked on
 * @param i First S2V2 calc'd index for Uvol offset
 * @param j Second S2V2 cald'd index for Uvol offset
 * @param TP_nVar Number of variables in 1st dim of Uvol. Used for offsets in Uvol.
 * @param Nloc Local polynomial order
 * @param ElemID Index for the element that the side being worked on. Used for offsets in Uvol.
 * @param Uvol Device pointer to volume data array
 * @param L_Minus L for boundary flux computation at plus side
 * @param L_Plus L for boundary flux computation at minus side
 * @param locSide Index for the side being operated on. WARNING! This is always passed as 1-indexed.
 * @remark Due to inability for C to slice arrays, this method works on only a single variable at once.
 **********************************************************************************************************************************/
__device__ double EvalElemFaceG_Device(int n, int i, int j, int TP_nVar, int Nloc, int ElemID, double* Uvol, double* L_Minus, double* L_Plus, int locSide)
{
    double Uface = 0.0;
    int Uvol_offset = 0;

    // locSide is passed as 1-indexed, so no modification is needed to match it to macro defined directions
    switch(locSide)
    {
        case XI_MINUS:
            for (int l = 0; l <= Nloc; l++)
            {
                Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, l, i, j);
                Uface += Uvol[Uvol_offset+n]*L_Minus[l];
            }
            break;
        case ETA_MINUS:
            for (int l = 0; l <= Nloc; l++)
            {
                Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, i, l, j);
                Uface += Uvol[Uvol_offset+n]*L_Minus[l];
            }
            break;
        case ZETA_MINUS:
            for (int l = 0; l <= Nloc; l++)
            {
                Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, i, j, l);
                Uface += Uvol[Uvol_offset+n]*L_Minus[l];
            }
            break;
        case XI_PLUS:
            for (int l = 0; l <= Nloc; l++)
            {
                Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, l, i, j);
                Uface += Uvol[Uvol_offset+n]*L_Plus[l];
            }
            break;
        case ETA_PLUS:
            for (int l = 0; l <= Nloc; l++)
            {
                Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, i, l, j);
                Uface += Uvol[Uvol_offset+n]*L_Plus[l];
            }
            break;
        case ZETA_PLUS:
            for (int l = 0; l <= Nloc; l++)
            {
                Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, i, j, l);
                Uface += Uvol[Uvol_offset+n]*L_Plus[l];
            }
            break;
    }

    return Uface;
}

/***********************************************************************************************************************************
 * @brief Interpolates the element volume data stored at Gauss-Lobatto points
 * @returns Value for Uface at U(n,:,:,:,ElemID)
 * @param n Index of flow variable being worked on
 * @param i First S2V2 calc'd index for Uvol offset
 * @param j Second S2V2 cald'd index for Uvol offset
 * @param TP_nVar Number of variables in 1st dim of Uvol. Only used to find offsets.
 * @param Nloc Local polynomial order
 * @param Uvol Device pointer to volume data array
 * @param locSide Index for the side being operated on. WARNING! This is always passed as 1-indexed.
 * @remark Due to inability for C to slice arrays, this method works on only a single variable at once.
 **********************************************************************************************************************************/
__device__ double EvalElemFaceGL_Device(int n, int i, int j, int TP_nVar, int Nloc, int ElemID, double* Uvol, int locSide)
{
    double Uface = 0.0;
    int Uvol_offset = 0;

    // locSide is passed as 1-indexed, so no modification is needed to match it to macro defined directions
    switch(locSide)
    {
        case XI_MINUS:
            Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, 0, i, j);
            Uface = Uvol[Uvol_offset+n];
            break;
        case ETA_MINUS:
            Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, i, 0, j);
            Uface = Uvol[Uvol_offset+n];
            break;
        case ZETA_MINUS:
            Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, i, j, 0);
            Uface = Uvol[Uvol_offset+n];
            break;
        case XI_PLUS:
            Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, Nloc, i, j);
            Uface = Uvol[Uvol_offset+n];
            break;
        case ETA_PLUS:
            Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, i, Nloc, j);
            Uface = Uvol[Uvol_offset+n];
            break;
        case ZETA_PLUS:
            Uvol_offset = FindPointOffset(TP_nVar, Nloc, ElemID, i, j, Nloc);
            Uface = Uvol[Uvol_offset+n];
            break;
    }

    return Uface;
}


/***********************************************************************************************************************************
 * @brief Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface integration points, 
 *        using fast 1D Interpolation and store in global side structure.
 **********************************************************************************************************************************/
__global__ void __launch_bounds__(256, 1) ProlongToFace_Kernel(int TP_nVar, int nDOFs, int Nloc, int firstSideID, double* Uvol, double* Uface_master, double* Uface_slave,
                                     double* L_Minus, double* L_Plus, int* SideToElem, int* S2V2
#if FV_ENABLED
                                     , int* isFV
#endif
                                    )
{
    double Uface;
    struct ThreadIndicesSide t_info;
    int thisDOF_offset = 0;
    int sideToElem_offset = 0;
    int SideToVol_offset = 0;
    int locSide, ElemID, nbElemID;
    int nbLocSide, flip;
    int i,j;
    
    t_info = FindThreadIndicesSide(Nloc, firstSideID);
    // FindThreadIndicesSide returns t_info.SideID as a 1-indexed, so need to account for that throughout
    thisDOF_offset = FindPointMasterSlaveOffset(TP_nVar, Nloc, t_info.SideID-1, t_info.p, t_info.q);

    // Find element that owns the side that owns this threads DOF then use that info to find the
    // offset for the volume array
    // IndexFlatFortranArr assumes 1-indexing, so don't modify t_info.SideID value
    sideToElem_offset = IndexFlatFortranArr(5, S2E_ELEM_ID, t_info.SideID);
    ElemID = SideToElem[sideToElem_offset];

    if (t_info.threadID > nDOFs) return;

    // Master sides
    if (ElemID > 0)
    {
        sideToElem_offset = IndexFlatFortranArr(5, S2E_LOC_SIDE_ID, t_info.SideID);
        locSide = SideToElem[sideToElem_offset];
        flip = 0;
        // We we make the final assignment, we will need locSide's offset in S2V2
        // locSide was assigned from SideToElem, so it is still the original 1-indexed Fortran value
        SideToVol_offset = FindPointMappingOffset(locSide-1, Nloc, 0, t_info.p, t_info.q, flip, true);
        i = S2V2[SideToVol_offset];
        j = S2V2[SideToVol_offset+1];

        for (int n = 0; n < TP_nVar; n++)
        {
#if FV_ENABLED
            if (PP_NodeType == 1 && isFV(ElemID) == 0)
            {
#else
            if (PP_NodeType == 1)
            {
#endif
                Uface = EvalElemFaceG_Device(n, i, j, TP_nVar, Nloc, ElemID, Uvol, L_Minus, L_Plus, locSide);
            }
            else
            {
                Uface = EvalElemFaceGL_Device(n, i, j, TP_nVar, Nloc, ElemID, Uvol, locSide);
            }

            Uface_master[thisDOF_offset+n] = Uface;
        }
    }

    // Repeat finding the volume element information, this time for the neighbor element to this thread's DOF's side
    sideToElem_offset = IndexFlatFortranArr(5, S2E_NB_ELEM_ID, t_info.SideID);
    nbElemID  = SideToElem[sideToElem_offset];

    // Slave side (ElemID,locSide and flip =-1 if not existing)
    if (nbElemID > 0)
    {
        sideToElem_offset = IndexFlatFortranArr(5, S2E_NB_LOC_SIDE_ID, t_info.SideID);
        nbLocSide = SideToElem[sideToElem_offset];
        sideToElem_offset = IndexFlatFortranArr(5, S2E_FLIP, t_info.SideID);
        flip = SideToElem[sideToElem_offset];
        SideToVol_offset = FindPointMappingOffset(nbLocSide-1, Nloc, 0, t_info.p, t_info.q, flip, true);
        i = S2V2[SideToVol_offset];
        j = S2V2[SideToVol_offset+1];

        for (int n = 0; n < TP_nVar; n++)
        {
#if FV_ENABLED
            if (PP_NodeType == 1 && isFV[nbElemID-1] == 0)
            {
#else
            if (PP_NodeType == 1)
            {
#endif
                Uface = EvalElemFaceG_Device(n, i, j, TP_nVar, Nloc, nbElemID, Uvol, L_Minus, L_Plus, nbLocSide);
            }
            else
            {
                Uface = EvalElemFaceGL_Device(n, i, j, TP_nVar, Nloc, nbElemID, Uvol, nbLocSide);
            }
            Uface_slave[thisDOF_offset+n] = Uface;
        }
        
    }
}


//----------------------------------------------------------------------------------------------------------------------------------
// Interface - Call kernels from HOST
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Device entry point for ProlongToFace kernel
 * @remark nSides and nElems are not used here, but are still passed for API consistency with the host backend.
 **********************************************************************************************************************************/
void ProlongToFace_Device(int TP_nVar, int Nloc, int nElems, int nSides, int firstSideID, int lastSideID,
                            int d_Uvol, int d_Uface_master, int d_Uface_slave, int d_L_Minus, int d_L_Plus, 
                            int d_SideToElem, int d_S2V2
#if FV_ENABLED
                            ,int d_isFV
#endif
                            , int streamID
                         )
{

#ifdef USE_NVTX /* If profiling, add range */
    nvtxRangePushA("ProlongToFace");
#endif

    int nDOFs = (Nloc+1)*(Nloc+1)*(lastSideID-firstSideID+1);
    INVOKE_KERNEL( ProlongToFace_Kernel, nDOFs/256+1, 256, 0, streams[streamID], TP_nVar, nDOFs, Nloc, firstSideID,
                                     (double*)DeviceVars[d_Uvol], (double*)DeviceVars[d_Uface_master], (double*)DeviceVars[d_Uface_slave],
                                     (double*)DeviceVars[d_L_Minus], (double*)DeviceVars[d_L_Plus],
                                     (int*)DeviceVars[d_SideToElem], (int*)DeviceVars[d_S2V2]
#if FV_ENABLED
                                     , (int*)DeviceVars[d_isFV]
#endif
                 );

#ifdef USE_NVTX
    nvtxRangePop();
#endif
}
