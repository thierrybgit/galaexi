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
  void SurfInt_Device(
#if FV_ENABLED
                    int d_FV_Elems_master, int d_FV_Elems_slave,
#endif
                    int Nloc, int nSides, int nElems, int firstMPISide_YOUR, int lastMPISide_MINE, bool doMPISides, int streamID,
                    int d_Flux_master, int d_Flux_slave, int d_Ut,
                    int d_L_HatMinus, int d_L_HatPlus,
                    int d_ElemToSide, int d_SideToElem, int d_S2V2
#if (PP_NodeType == 1 && defined(SPLIT_DG))
                    ,int d_U, int d_UPrim,
                    int d_U_master, int d_UPrim_master,
                    int d_U_slave, int d_UPrim_slave,
                    int d_L_minus, int d_L_plus,
                    int d_S2V, int d_Metrics_fTilde, 
                    int d_Metrics_gTilde, int d_Metrics_hTilde, int d_Ja_Face, int d_Ja_slave
#endif
                    );
}

/***********************************************************************************************************************************
 * @brief Compute the surface integral for either unsplit Gauss nodes or Gauss-Lobatto nodes
 **********************************************************************************************************************************/
__device__ void DoSurfInt_Device(int Nloc, int ElemID, int SideID, int locSideID, int flip, bool isMaster, int i, int j, int k, int* S2V2,
                                 double* Flux_master, double* Flux_slave, double* Ut, double* L_HatMinus, double* L_HatPlus
                                )
{
  // NOTE: SINCE, RIGHT NOW, THE ONLY WAY THIS IS ACTUALLY CALLED IS VIA SurfIntCons WHERE TP_nVars = PP_nVar
  double Ut_loc[PP_nVar] = {0.0,0.0,0.0,0.0,0.0};
  int p,q;
  int mappingOffset = 0;
  int fluxOffset = 0;
  int UtOffset = 0;

  switch (locSideID+1)
  {
    case XI_MINUS:
#if (PP_NodeType == 2 && !defined(EXACT_MM))
      if (i != 0) return;
#endif
      mappingOffset = FindPointMappingOffset(locSideID, Nloc, 0, j, k, flip, true);
      p = S2V2[mappingOffset];
      q = S2V2[mappingOffset + 1];
      fluxOffset = FindPointMasterSlaveOffset(PP_nVar, Nloc, SideID, p, q);
      for (int n = 0; n < PP_nVar; n++)
      {
        if (isMaster)
        {
          Ut_loc[n] += Flux_master[fluxOffset+n]*L_HatMinus[i];
        }
        else // slave side
        {
          Ut_loc[n] -= Flux_slave[fluxOffset+n]*L_HatMinus[i];
        }
      }
      break;
    case ETA_MINUS:
#if (PP_NodeType == 2 && !defined(EXACT_MM))
      if (j != 0) return;
#endif
      mappingOffset = FindPointMappingOffset(locSideID, Nloc, 0, i, k, flip, true);
      p = S2V2[mappingOffset];
      q = S2V2[mappingOffset + 1];
      fluxOffset = FindPointMasterSlaveOffset(PP_nVar, Nloc, SideID, p, q);
      for (int n = 0; n < PP_nVar; n++)
      {
        if (isMaster)
        {
          Ut_loc[n] += Flux_master[fluxOffset+n]*L_HatMinus[j];
        }
        else // slave side
        {
          Ut_loc[n] -= Flux_slave[fluxOffset+n]*L_HatMinus[j];
        }
      }
      break;
    case ZETA_MINUS:
#if (PP_NodeType == 2 && !defined(EXACT_MM))
      if (k != 0) return;
#endif
      mappingOffset = FindPointMappingOffset(locSideID, Nloc, 0, i, j, flip, true);
      p = S2V2[mappingOffset];
      q = S2V2[mappingOffset + 1];
      fluxOffset = FindPointMasterSlaveOffset(PP_nVar, Nloc, SideID, p, q);
      for (int n = 0; n < PP_nVar; n++)
      {
        if (isMaster)
        {
          Ut_loc[n] += Flux_master[fluxOffset+n]*L_HatMinus[k];
        }
        else // slave side
        {
          Ut_loc[n] -= Flux_slave[fluxOffset+n]*L_HatMinus[k];
        }
      }
      break;
    case XI_PLUS:
#if (PP_NodeType == 2 && !defined(EXACT_MM))
      if (i != Nloc) return;
#endif
      mappingOffset = FindPointMappingOffset(locSideID, Nloc, 0, j, k, flip, true);
      p = S2V2[mappingOffset];
      q = S2V2[mappingOffset + 1];
      fluxOffset = FindPointMasterSlaveOffset(PP_nVar, Nloc, SideID, p, q);
      for (int n = 0; n < PP_nVar; n++)
      {
        if (isMaster)
        {
          Ut_loc[n] += Flux_master[fluxOffset+n]*L_HatPlus[i];
        }
        else // slave side
        {
          Ut_loc[n] -= Flux_slave[fluxOffset+n]*L_HatPlus[i];
        }
      }
      break;
    case ETA_PLUS:
#if (PP_NodeType == 2 && !defined(EXACT_MM))
      if (j != Nloc) return;
#endif
      mappingOffset = FindPointMappingOffset(locSideID, Nloc, 0, i, k, flip, true);
      p = S2V2[mappingOffset];
      q = S2V2[mappingOffset + 1];
      fluxOffset = FindPointMasterSlaveOffset(PP_nVar, Nloc, SideID, p, q);
      for (int n = 0; n < PP_nVar; n++)
      {
        if (isMaster)
        {
          Ut_loc[n] += Flux_master[fluxOffset+n]*L_HatPlus[j];
        }
        else // slave side
        {
          Ut_loc[n] -= Flux_slave[fluxOffset+n]*L_HatPlus[j];
        }
      }
      break;
    case ZETA_PLUS:
#if (PP_NodeType == 2 && !defined(EXACT_MM))
      if (k != Nloc) return;
#endif
      mappingOffset = FindPointMappingOffset(locSideID, Nloc, 0, i, j, flip, true);
      p = S2V2[mappingOffset];
      q = S2V2[mappingOffset + 1];
      fluxOffset = FindPointMasterSlaveOffset(PP_nVar, Nloc, SideID, p, q);
      for (int n = 0; n < PP_nVar; n++)
      {
        if (isMaster)
        {
          Ut_loc[n] += Flux_master[fluxOffset+n]*L_HatPlus[k];
        }
        else // slave side
        {
          Ut_loc[n] -= Flux_slave[fluxOffset+n]*L_HatPlus[k];
        }
      }
      break;
  }

  UtOffset = FindPointOffset(PP_nVar, Nloc, ElemID, i, j, k);
  for (int n = 0; n < PP_nVar; n++)
  {
    Ut[UtOffset+n] += Ut_loc[n];
  }
}

//----------------------------------------------------------------------------------------------------------------------------------
// SplitDG_Info data structure to simplify calculation of offsets for arrays during prolog step of SurfInt with Gauss nodes and SplitDG
//----------------------------------------------------------------------------------------------------------------------------------
#if (PP_NodeType == 1 && defined(SPLIT_DG))
/***********************************************************************************************************************************
 * Struct to hold all data needed to pass information to SplitFlux method
 * when using SplitDG with Gauss nodes. This allows all of this information
 * to be set in the function SetSplitDGInfo and returned easily
 **********************************************************************************************************************************/
struct SplitDG_Info
{
  int Metrics_offset = 0;
  int Mapping_offset = 0;
  int U_offset = 0;
  int Uside_offset = 0;
  int UPrim_offset = 0;
  int UPrimSide_offset = 0;
  int Ja_idx = 0;
  int ijk[3];

  double Ja_Face_loc[3];

  double* Metrics;
  double* Leg;
  double* LegHat;
};

/***********************************************************************************************************************************
 * @brief Set up all information needed to pass to SplitFlux formulation when using Gauss nodes
 **********************************************************************************************************************************/
__device__ SplitDG_Info SetSplitDGInfo(int locSideID, int SideID, int Nloc, int ElemID, int nElems, int l, int p, int q, int flip, 
                                       double* Metrics_fTilde, double* Metrics_gTilde, double* Metrics_hTilde,
                                       double* L_minus, double* L_plus, double* L_HatMinus, double* L_HatPlus,
                                       double* Ja_Face, int* S2V
                                      )
{
  struct SplitDG_Info newInfo;
  int Ja_offset = 0;

  // Find offsets for flattened memory
  // locSideID and SideID passed as 0-indexed
  // ElemID passed as 1-indexed
  newInfo.Mapping_offset = FindPointMappingOffset(locSideID, Nloc, l, p, q, flip, false);
  newInfo.ijk[0] = S2V[newInfo.Mapping_offset];
  newInfo.ijk[1] = S2V[newInfo.Mapping_offset+1];
  newInfo.ijk[2] = S2V[newInfo.Mapping_offset+2];
  newInfo.Metrics_offset = FindPointMetricsOffset(Nloc, ElemID, nElems, 0, newInfo.ijk[0], newInfo.ijk[1], newInfo.ijk[2]);
  newInfo.U_offset = FindPointOffset(PP_nVar, Nloc, ElemID, newInfo.ijk[0], newInfo.ijk[1], newInfo.ijk[2]);
  newInfo.Uside_offset = FindPointMasterSlaveOffset(PP_nVar, Nloc, SideID, p, q);
  newInfo.UPrim_offset = FindPointOffset(PP_nVarPrim, Nloc, ElemID, newInfo.ijk[0], newInfo.ijk[1], newInfo.ijk[2]);
  newInfo.UPrimSide_offset = FindPointMasterSlaveOffset(PP_nVarPrim, Nloc, SideID, p, q);

  switch(locSideID+1) // Add one to match directives which are defined assuming Fortran indexing
  {
    case XI_MINUS:
      // Pull the surface metrics we want out of flat memory
      newInfo.Ja_idx = 0;
      for (int m = 0; m < 3; m++)
      {
        Ja_offset = SideID*(Nloc+1)*(Nloc+1)*3*3 + q*(Nloc+1)*3*3 + p*3*3 + m*3 + newInfo.Ja_idx; // This is the offset to the start of 3x3 first 2 dimensions
        newInfo.Ja_Face_loc[m] = Ja_Face[Ja_offset];
      }
      newInfo.Metrics = Metrics_fTilde;
      newInfo.Leg = L_minus;
      newInfo.LegHat = L_HatMinus;
      break;
    case ETA_MINUS:
      newInfo.Ja_idx = 1;
      for (int m = 0; m < 3; m++)
      {
        Ja_offset = SideID*(Nloc+1)*(Nloc+1)*3*3 + q*(Nloc+1)*3*3 + p*3*3 + m*3 + newInfo.Ja_idx; // This is the offset to the start of 3x3 first 2 dimensions
        newInfo.Ja_Face_loc[m] = Ja_Face[Ja_offset];
      }
      newInfo.Metrics = Metrics_gTilde;
      newInfo.Leg = L_minus;
      newInfo.LegHat = L_HatMinus;
      break;
    case ZETA_MINUS:
      newInfo.Ja_idx = 2;
      for (int m = 0; m < 3; m++)
      {
        Ja_offset = SideID*(Nloc+1)*(Nloc+1)*3*3 + q*(Nloc+1)*3*3 + p*3*3 + m*3 + newInfo.Ja_idx; // This is the offset to the start of 3x3 first 2 dimensions
        newInfo.Ja_Face_loc[m] = Ja_Face[Ja_offset];
      }
      newInfo.Metrics = Metrics_hTilde;
      newInfo.Leg = L_minus;
      newInfo.LegHat = L_HatMinus;
      break;
    case XI_PLUS:
      newInfo.Ja_idx = 0;
      for (int m = 0; m < 3; m++)
      {
        Ja_offset = SideID*(Nloc+1)*(Nloc+1)*3*3 + q*(Nloc+1)*3*3 + p*3*3 + m*3 + newInfo.Ja_idx; // This is the offset to the start of 3x3 first 2 dimensions
        newInfo.Ja_Face_loc[m] = Ja_Face[Ja_offset];
      }
      newInfo.Metrics = Metrics_fTilde;
      newInfo.Leg = L_plus;
      newInfo.LegHat = L_HatPlus;
      break;
    case ETA_PLUS:
      newInfo.Ja_idx = 1;
      for (int m = 0; m < 3; m++)
      {
        Ja_offset = SideID*(Nloc+1)*(Nloc+1)*3*3 + q*(Nloc+1)*3*3 + p*3*3 + m*3 + newInfo.Ja_idx; // This is the offset to the start of 3x3 first 2 dimensions
        newInfo.Ja_Face_loc[m] = Ja_Face[Ja_offset];
      }
      newInfo.Metrics = Metrics_gTilde;
      newInfo.Leg = L_plus;
      newInfo.LegHat = L_HatPlus;
      break;
    case ZETA_PLUS:
      newInfo.Ja_idx = 2;
      for (int m = 0; m < 3; m++)
      {
        Ja_offset = SideID*(Nloc+1)*(Nloc+1)*3*3 + q*(Nloc+1)*3*3 + p*3*3 + m*3 + newInfo.Ja_idx; // This is the offset to the start of 3x3 first 2 dimensions
        newInfo.Ja_Face_loc[m] = Ja_Face[Ja_offset];
      }
      newInfo.Metrics = Metrics_hTilde;
      newInfo.Leg = L_plus;
      newInfo.LegHat = L_HatPlus;
      break;
  }

  return newInfo;
}

/***********************************************************************************************************************************
 * @brief Compute the surface integral on the device for Split Gauss
 * @remark This kernel is SIDE parallel, where as the other 
 **********************************************************************************************************************************/
__global__ void __launch_bounds__(128, 1) SurfInt_Kernel(
#if FV_ENABLED
                                int* FV_Elems_master, int* FV_Elems_slave,
#endif
                                int Nloc, int nSides, int nElems, int firstSideID, int lastSideID,
                                double* Flux_master, double* Flux_slave, double* Ut,
                                double* L_HatMinus, double* L_HatPlus,
                                int* ElemToSide,  int* SideToElem, int* S2V2
                                ,double* U, double* UPrim, 
                                double* U_master, double* UPrim_master,
                                double* U_slave, double* UPrim_slave,
                                double* L_minus, double* L_plus,
                                int* S2V,
                                double* Metrics_fTilde, double* Metrics_gTilde,
                                double* Metrics_hTilde, double* Ja_Face, double* Ja_slave
                                )
{
  // Declare local variables C style
  struct ThreadIndicesSide t_info;
  int ElemID, flip;
  int ElemToSide_offset = 0;
  struct SplitDG_Info splitInfo;
  double FluxB_sum[PP_nVar] = {0.0,0.0,0.0,0.0,0.0};
  double FluxB;
  int SideToElem_offset = 0;
  int nbElemID = 0;
  int locSideID = 0;

  // Select the DOF indices for the local thread
  // Because it is easier for the other SurfInt method, firstSideID and lastSideID are both adjust to 0-indexed
  // So account for that here.
  t_info = FindThreadIndicesSide(Nloc, firstSideID+1);
  if (t_info.SideID > lastSideID+1) return;

  SideToElem_offset = IndexFlatFortranArr(5, S2E_ELEM_ID, t_info.SideID);
  ElemID = SideToElem[SideToElem_offset];
  SideToElem_offset = IndexFlatFortranArr(5,S2E_NB_ELEM_ID,t_info.SideID);
  nbElemID = SideToElem[SideToElem_offset];

#if FV_ENABLED
  if (FV_Elems_master[t_info.SideID-1] == 0) // Only do surface integral if this is a DG element
  {
#endif /* FV_ENABLED */

    if (ElemID > 0)
    {
      SideToElem_offset = IndexFlatFortranArr(5,S2E_LOC_SIDE_ID,t_info.SideID);
      locSideID = SideToElem[SideToElem_offset];
      flip = 0;

      // Reset sum vector for the new side
      for (int n = 0; n < PP_nVar; n++)
      {
        FluxB_sum[n] = 0.0;
      }

      // Find sum of fluxes for all DOFs along the i-line
      // FluxB_sum is thread local
      for (int l = 0; l <= Nloc; l++)
      {
        // Do logic to set metrics pointer for SPLIT_DG
        splitInfo = SetSplitDGInfo(locSideID-1, t_info.SideID-1, Nloc, ElemID, nElems, l, t_info.p, t_info.q, flip, 
                                    Metrics_fTilde, Metrics_gTilde, Metrics_hTilde,
                                    L_minus, L_plus, L_HatMinus, L_HatPlus,
                                    Ja_Face, S2V
                                  );

        // Sum the split flux -- "Similar to prolong to face of the two-point flux"         
        for (int n = 0; n < PP_nVar; n++)
        {
          FluxB = SplitDGVolume_Pointer_Device(n, &U[splitInfo.U_offset], &UPrim[splitInfo.UPrim_offset],
                                                                &U_master[splitInfo.Uside_offset], &UPrim_master[splitInfo.UPrimSide_offset],
                                                                &splitInfo.Metrics[splitInfo.Metrics_offset], &splitInfo.Ja_Face_loc[0]
                                              );
          FluxB_sum[n] +=  FluxB * splitInfo.Leg[splitInfo.ijk[splitInfo.Ja_idx]];
        }
      }

      for (int l = 0; l <= Nloc; l++)
      {
        // Do logic to set metrics pointers to global indices
        splitInfo = SetSplitDGInfo(locSideID-1, t_info.SideID-1, Nloc, ElemID, nElems, l, t_info.p, t_info.q, flip, 
                                    Metrics_fTilde, Metrics_gTilde, Metrics_hTilde,
                                    L_minus, L_plus, L_HatMinus, L_HatPlus,
                                    Ja_Face, S2V
                                  );

        for (int n = 0; n < PP_nVar; n++)
        {
          FluxB = SplitDGVolume_Pointer_Device(n, &U[splitInfo.U_offset], &UPrim[splitInfo.UPrim_offset],
                                                                &U_master[splitInfo.Uside_offset], &UPrim_master[splitInfo.UPrimSide_offset],
                                                                &splitInfo.Metrics[splitInfo.Metrics_offset], &splitInfo.Ja_Face_loc[0]
                                              );

          atomicAdd(&Ut[splitInfo.U_offset+n], (Flux_master[splitInfo.Uside_offset+n] - 
                                                       0.5*(FluxB_sum[n] - FluxB) * d_NormalSigns[locSideID-1]) * splitInfo.LegHat[splitInfo.ijk[splitInfo.Ja_idx]]);
        }
      }
    }

#if FV_ENABLED
  } // if (Fv_Elems_master[locSideID] == 0)

  if (FV_Elems_slave[t_info.SideID-1] == 0) // Only account for neighbor if it is a DG element
  {
#endif /* FV_ENABLED */

    if (nbElemID > 0)
    {

      SideToElem_offset = IndexFlatFortranArr(5,S2E_NB_LOC_SIDE_ID,t_info.SideID);
      locSideID = SideToElem[SideToElem_offset];
      SideToElem_offset = IndexFlatFortranArr(5,S2E_FLIP,t_info.SideID);
      flip = SideToElem[SideToElem_offset];

      // Reset sum vector for the new side
      for (int n = 0; n < PP_nVar; n++)
      {
        FluxB_sum[n] = 0.0;
      }

      // Find sum of fluxs for all DOFs along the i-line
      // FluxB_sum is thread local
      for (int l = 0; l <= Nloc; l++)
      {
        // Do logic to set metrics pointer for SPLIT_DG
        splitInfo = SetSplitDGInfo(locSideID-1, t_info.SideID-1, Nloc, nbElemID, nElems, l, t_info.p, t_info.q, flip, 
                                    Metrics_fTilde, Metrics_gTilde, Metrics_hTilde,
                                    L_minus, L_plus, L_HatMinus, L_HatPlus,
                                    Ja_slave, S2V
                                  );

        // Sum the split flux -- "Similar to prolong to face of the two-point flux"         
        for (int n = 0; n < PP_nVar; n++)
        {

          FluxB = SplitDGVolume_Pointer_Device(n, &U[splitInfo.U_offset], &UPrim[splitInfo.UPrim_offset],
                                                &U_slave[splitInfo.Uside_offset], &UPrim_slave[splitInfo.UPrimSide_offset],
                                                &splitInfo.Metrics[splitInfo.Metrics_offset], &splitInfo.Ja_Face_loc[0]
                                              );
          FluxB_sum[n] += FluxB * splitInfo.Leg[splitInfo.ijk[splitInfo.Ja_idx]];
        }
      }

      for (int l = 0; l <= Nloc; l++)
      {
        // Do logic to set metrics pointers to global indices
        splitInfo = SetSplitDGInfo(locSideID-1, t_info.SideID-1, Nloc, nbElemID, nElems, l, t_info.p, t_info.q, flip, 
                                    Metrics_fTilde, Metrics_gTilde, Metrics_hTilde,
                                    L_minus, L_plus, L_HatMinus, L_HatPlus,
                                    Ja_slave, S2V
                                  );

        // With full FluxB_sum value available for the local thread calculate d_Ut atomically using full element context
        // This is functionally the same as summing over i,j,k on whole element
        for (int n = 0; n < PP_nVar; n++)
        {
          FluxB = SplitDGVolume_Pointer_Device(n, &U[splitInfo.U_offset], &UPrim[splitInfo.UPrim_offset],
                                                &U_slave[splitInfo.Uside_offset], &UPrim_slave[splitInfo.UPrimSide_offset],
                                                &splitInfo.Metrics[splitInfo.Metrics_offset], &splitInfo.Ja_Face_loc[0]
                                              );

          atomicAdd(&Ut[splitInfo.U_offset+n], (-Flux_slave[splitInfo.Uside_offset+n] - 
                                                         0.5*(FluxB_sum[n] - FluxB) * d_NormalSigns[locSideID-1]) * splitInfo.LegHat[splitInfo.ijk[splitInfo.Ja_idx]]);
        }
      }
    }

#if FV_ENABLED
  // LEAVE OUT FV stuff for now, will add that back in later when doing FV work
  } // if (Fv_Elems_slave[locSideID] == 0)
#endif /* FV_ENABLED */
}

#else /* ! (PP_NodeType == 1 && defined(SPLIT_DG)) */

/***********************************************************************************************************************************
 * @brief Compute the surface integral on the device for Unsplit Gauss or Split Gauss-Lobatto
 **********************************************************************************************************************************/
__global__ void __launch_bounds__(128, 1) SurfInt_Kernel(
#if FV_ENABLED
                                int* FV_Elems_master, int* FV_Elems_slave,
#endif
                                int Nloc, int nSides, int nElems, int firstSideID, int lastSideID,
                                double* Flux_master, double* Flux_slave, double* Ut,
                                double* L_HatMinus, double* L_HatPlus,
                                int* ElemToSide,  int* SideToElem, int* S2V2
                                )
{
  // Declare local variables C style
  struct ThreadIndicesVolume t_info;
  int SideID, flip;
  int ElemToSide_offset = 0;
  bool isMaster;

  // Select the DOF indices for the local thread
  // Index threads from 1 as having a thread 0 causing i = -1 for that thread
  t_info = FindThreadIndicesVolume(Nloc);

  if (t_info.ElemID <= nElems)
  {
    // Loop through sides of thread local element
    // NOTE:: firstSideID and lastSideID were shifted to 0-base index in wrapper method
    for (int locSideID = 0; locSideID < 6; locSideID++)
    {
      // MASTER SIDES
      ElemToSide_offset = IndexFlatFortranArr(2, 6, E2S_SIDE_ID, locSideID+1, t_info.ElemID);
      // Subtract 1 to make sure that both side indices are indexed from 0 for consistency
      // It also makes it less verbose in calculating offsets
      SideID = ElemToSide[ElemToSide_offset] - 1; 
      if (SideID < firstSideID || SideID > lastSideID) continue;
      ElemToSide_offset = IndexFlatFortranArr(2, 6, E2S_FLIP, locSideID+1, t_info.ElemID);
      flip = ElemToSide[ElemToSide_offset];

      isMaster = flip == 0;

#if FV_ENABLED
    if ((isMaster && FV_Elems_master[SideID] == 0) || (! isMaster && FV_Elems_slave[SideID] == 0)) // Only do surface integral if this is a DG element
    {
#endif /* FV_ENABLED */

      DoSurfInt_Device(Nloc, t_info.ElemID, SideID, locSideID, flip, isMaster, t_info.i, t_info.j, t_info.k, S2V2, Flux_master, Flux_slave, Ut, L_HatMinus, L_HatPlus);

#if FV_ENABLED
    // LEAVE OUT FV stuff for now, will add that back in later when doing FV work
    }
    else
    {
    }      
#endif /* FV_ENABLED */

    }

  } // if (t_info.ElemID <= nElems)

}
#endif /* PP_NodeType == 1 && defined(SPLIT_DG) */

#if (FV_ENABLED == 2) /* Blending FV */
  // LEAVE OUT FV stuff for now, will add that back in later when doing FV work
  // THIS WILL BE THE SURFINT_BLEND METHOD -- INVESTIGATE IF THIS CAN BE MERGED INTO THE ABOVE METHOD CLEANLY
#endif /* Blending FV */

/***********************************************************************************************************************************
 * @brief Entry-point method from Fortran for the surface integral calculation that calls C kernel
 **********************************************************************************************************************************/
void SurfInt_Device(
#if FV_ENABLED
                    int d_FV_Elems_master, int d_FV_Elems_slave,
#endif
                    int Nloc, int nSides, int nElems, int firstMPISide_YOUR, int lastMPISide_MINE, bool doMPISides, int streamID,
                    int d_Flux_master, int d_Flux_slave, int d_Ut,
                    int d_L_HatMinus, int d_L_HatPlus,
                    int d_ElemToSide,  int d_SideToElem, int d_S2V2
#if (PP_NodeType == 1 && defined(SPLIT_DG))
                    ,int d_U, int d_UPrim,
                    int d_U_master, int d_UPrim_master,
                    int d_U_slave, int d_UPrim_slave,
                    int d_L_minus, int d_L_plus,
                    int d_S2V, int d_Metrics_fTilde, 
                    int d_Metrics_gTilde, int d_Metrics_hTilde, int d_Ja_Face, int d_Ja_slave
#endif
                    )
{

  int firstSideID = 0;
  int lastSideID = 0;
  int nDOF = 0;

  // Do logic for first AND lastSideID here so consistent across threads on device
  if (doMPISides)
  {
    firstSideID = firstMPISide_YOUR - 1;
    lastSideID = nSides - 1;
  }
  else
  {
    firstSideID = 0;
    lastSideID = lastMPISide_MINE - 1;
  }

#ifdef USE_NVTX /* If profiling, add range */
  nvtxRangePushA("SurInt");
#endif

// When using split Gauss nodes, the kernel is side parallel
// When using split GL or unsplit Gauss, the kernel is element parallel
#if (PP_NodeType == 1 && defined(SPLIT_DG))
  nDOF = nSides*(Nloc+1)*(Nloc+1);
#else
  nDOF = nElems*(Nloc+1)*(Nloc+1)*(Nloc+1);
#endif

  INVOKE_KERNEL( SurfInt_Kernel, nDOF/128+1, 128, 0, streams[streamID],
#if FV_ENABLED
                                (int*)DeviceVars[d_FV_Elems_master], (int*)DeviceVars[d_FV_Elems_slave],
#endif
                                Nloc, nSides, nElems, firstSideID, lastSideID,
                                (double*)DeviceVars[d_Flux_master], (double*)DeviceVars[d_Flux_slave], (double*)DeviceVars[d_Ut],
                                (double*)DeviceVars[d_L_HatMinus], (double*)DeviceVars[d_L_HatPlus],
                                (int*)DeviceVars[d_ElemToSide], (int*)DeviceVars[d_SideToElem], (int*)DeviceVars[d_S2V2] 
#if (PP_NodeType == 1 && defined(SPLIT_DG))
                                ,(double*)DeviceVars[d_U], (double*)DeviceVars[d_UPrim], 
                                (double*)DeviceVars[d_U_master], (double*)DeviceVars[d_UPrim_master],
                                (double*)DeviceVars[d_U_slave], (double*)DeviceVars[d_UPrim_slave],
                                (double*)DeviceVars[d_L_minus], (double*)DeviceVars[d_L_plus],
                                (int*)DeviceVars[d_S2V],
                                (double*)DeviceVars[d_Metrics_fTilde], (double*)DeviceVars[d_Metrics_gTilde],
                                (double*)DeviceVars[d_Metrics_hTilde], (double*)DeviceVars[d_Ja_Face], (double*)DeviceVars[d_Ja_slave]
#endif
                                );

#ifdef USE_NVTX /* If profiling, end range */
  nvtxRangePop(); 
#endif

}
