/*================================================================================================================================
  Copyright (c) 2022-2024 Prof. Andrea Beck
  Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
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
    void ApplyDMatrix_Device(int Nloc, int nElems, int d_Ut, int d_F, int d_G, 
                            int d_H, int d_D_Hat_T, int streamID, bool doOverwrite);
}

//----------------------------------------------------------------------------------------------------------------------------------
// Kernels
//----------------------------------------------------------------------------------------------------------------------------------
__global__ void ApplyDMatrix_Kernel(int Nloc, int nElems, double* Ut, double* F, double* G, 
                            double* H, double* D_Hat_T, bool doOverwrite)
{
    struct ThreadIndicesVolume t_info;
    int D_offsets[3] = {0,0,0};
    int Flux_offsets[3] = {0,0,0};
    int Ut_offset = 0;
    double Ut_loc[PP_nVar] = {0.0,0.0,0.0,0.0,0.0};

    // Get ElemID and i,j,k indices for the DOF assigned to this device thread.
    t_info = FindThreadIndicesVolume(Nloc);

    if (t_info.ElemID <= nElems)
    {
        for (int l = 0; l <= Nloc; l++)
        {
            // Get flattened memory offsets (in bytes) for D_Hat_T
            // Must +1, as IndexFlatFortranArr assumes array is indexed from 1 on Fortran side but
            // D_Hat_T is allocated as 0:PP_N
            D_offsets[0] = IndexFlatFortranArr(Nloc+1, l+1, t_info.i+1);
            D_offsets[1] = IndexFlatFortranArr(Nloc+1, l+1, t_info.j+1);
            D_offsets[2] = IndexFlatFortranArr(Nloc+1, l+1, t_info.k+1);
            for (int n = 0; n < PP_nVar; n++)
            {
                // Get flattened memory offsets (in bytes) for fluxes
                Flux_offsets[0] = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, l, t_info.j, t_info.k) + n;
                Flux_offsets[1] = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, t_info.i, l, t_info.k) + n;
                Flux_offsets[2] = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, t_info.i, t_info.j, l) + n;
                
                Ut_loc[n] += D_Hat_T[D_offsets[0]]*F[Flux_offsets[0]]
                            + D_Hat_T[D_offsets[1]]*G[Flux_offsets[1]]
                            + D_Hat_T[D_offsets[2]]*H[Flux_offsets[2]];
            }
        }

        // Get flattened memory for the DOF of Ut assigned to this thread for the final Ut calc
        Ut_offset = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, t_info.i, t_info.j, t_info.k);
        for (int n = 0; n < PP_nVar; n++)
        {
            if (doOverwrite)
            {
                Ut[Ut_offset+n] = Ut_loc[n];
            }
            else
            {
                Ut[Ut_offset+n] += Ut_loc[n];
            }
        }
    }
}


//----------------------------------------------------------------------------------------------------------------------------------
// Interface - Call kernels from HOST
//----------------------------------------------------------------------------------------------------------------------------------
void ApplyDMatrix_Device(int Nloc, int nElems, int d_Ut, int d_F, int d_G, 
                            int d_H, int d_D_Hat_T, int streamID, bool doOverwrite)
{
#ifdef USE_NVTX /* If profiling, add range */
    nvtxRangePushA("ApplyDMatrix");
#endif
    
    int nDOF = (Nloc+1)*(Nloc+1)*(Nloc+1)*nElems;
    INVOKE_KERNEL( ApplyDMatrix_Kernel,  nDOF/256+1, 256, 0, streams[streamID], Nloc, nElems,
                                    (double*)DeviceVars[d_Ut], (double*)DeviceVars[d_F], (double*)DeviceVars[d_G],
                                    (double*)DeviceVars[d_H], (double*)DeviceVars[d_D_Hat_T], doOverwrite
                  );

#ifdef USE_NVTX
    nvtxRangePop();
#endif
}