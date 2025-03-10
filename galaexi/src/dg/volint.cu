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
    void VolInt_splitForm_Device(int Nloc, int nElems, int d_Ut, int d_U, int d_UPrim, int d_Mf, int d_Mg,
                                 int d_Mh, int d_DVolSurf
                               , int streamID);
}

__global__ void __launch_bounds__(256, 1) VolInt_splitForm_Kernel(int Nloc, int nElems, double* Ut, double* U, double* UPrim,
                                        double* Mf, double* Mg, double* Mh, double* DVolSurf) {
    
struct ThreadIndicesVolume t_info;
t_info = FindThreadIndicesVolume(Nloc);

int U_offset = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, t_info.i, t_info.j, t_info.k);
int UPrim_offset = FindPointOffset(PP_nVarPrim, Nloc, t_info.ElemID, t_info.i, t_info.j, t_info.k);
int Metric_offset = FindPointMetricsOffset(Nloc, t_info.ElemID, nElems, 0, t_info.i, t_info.j, t_info.k);

if (t_info.ElemID <= nElems) {

  for(int n = 0; n < PP_nVar; n++) {
#if FV_ENABLED || !PARABOLIC
    Ut[U_offset + n] = 0.0;
#endif
    double Ut_tmp = 0.0;
    for(int l = 0; l <= Nloc; l++) {
      int VolSurf_offsetx = IndexFlatFortranArr(Nloc+1, l+1, t_info.i+1);
      int U_offsetx = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, l, t_info.j, t_info.k);
      int UPrim_offsetx = FindPointOffset(PP_nVarPrim, Nloc, t_info.ElemID, l, t_info.j, t_info.k);
      int Metric_offsetx = FindPointMetricsOffset(Nloc, t_info.ElemID, nElems, 0, l, t_info.j, t_info.k);
      int VolSurf_offsety = IndexFlatFortranArr(Nloc+1, l+1, t_info.j+1);
      int U_offsety = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, t_info.i, l, t_info.k);
      int UPrim_offsety = FindPointOffset(PP_nVarPrim, Nloc, t_info.ElemID, t_info.i, l, t_info.k);
      int Metric_offsety = FindPointMetricsOffset(Nloc, t_info.ElemID, nElems, 0, t_info.i, l, t_info.k);
      int VolSurf_offsetz = IndexFlatFortranArr(Nloc+1, l+1, t_info.k+1);
      int U_offsetz = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, t_info.i, t_info.j, l);
      int UPrim_offsetz = FindPointOffset(PP_nVarPrim, Nloc, t_info.ElemID, t_info.i, t_info.j, l);
      int Metric_offsetz = FindPointMetricsOffset(Nloc, t_info.ElemID, nElems, 0, t_info.i, t_info.j, l);

      Ut_tmp += DVolSurf[VolSurf_offsetx]*SplitDGVolume_Pointer_Device(n,&U[U_offset], &UPrim[UPrim_offset],
                &U[U_offsetx], &UPrim[UPrim_offsetx], &Mf[Metric_offset], &Mf[Metric_offsetx]);
      Ut_tmp += DVolSurf[VolSurf_offsety]*SplitDGVolume_Pointer_Device(n,&U[U_offset], &UPrim[UPrim_offset],
                &U[U_offsety], &UPrim[UPrim_offsety], &Mg[Metric_offset], &Mg[Metric_offsety]);
      Ut_tmp += DVolSurf[VolSurf_offsetz]*SplitDGVolume_Pointer_Device(n,&U[U_offset], &UPrim[UPrim_offset],
                &U[U_offsetz], &UPrim[UPrim_offsetz], &Mh[Metric_offset], &Mh[Metric_offsetz]);
    }
    Ut[U_offset + n] += Ut_tmp;
  }
}
}

//----------------------------------------------------------------------------------------------------------------------------------
// Kernel
//----------------------------------------------------------------------------------------------------------------------------------
// THIS METHOD CURRENTLY DOES NOT INCLUDE HANDLING FOR ANY OF THE FV STUFF, WHICH WILL BE INCORPORATED WHEN FV WORK IS DONE
__global__ void __launch_bounds__(256, 1) VolInt_splitForm_KernelA(int Nloc, int nElems, double* Ut, double* U, double* UPrim,
                                        double* Mf, double* Mg, double* Mh, double* DVolSurf
                                        )
{
    struct ThreadIndicesVolume t_info;
    double Ut_tmp[PP_nVar] = {0.0,0.0,0.0,0.0,0.0};
    int VolSurf_offset = 0;
    int U_offset1 = 0;
    int U_offset2 = 0;
    int UPrim_offset1 = 0;
    int UPrim_offset2 = 0;
    int Metric_offset1 = 0;
    int Metric_offset2 = 0;

    // Get ElemID and i,j,k indices for the DOF assigned to this device thread.
    t_info = FindThreadIndicesVolume(Nloc);

    // There are two offsets for the variable arrays and metrics
    // offset1 --> corresponds to the values passed to "ref" args of the SplitFlux method
    // offset2 --> corresponds to the vlaues passed to "non-ref" args
    // offset1's are calc'd here as they are dependent on DOF indices for this device thread
    // offset2's are calc'd inside the for loops below
    U_offset1 = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, t_info.i, t_info.j, t_info.k);
    UPrim_offset1 = FindPointOffset(PP_nVarPrim, Nloc, t_info.ElemID, t_info.i, t_info.j, t_info.k);
    Metric_offset1 = FindPointMetricsOffset(Nloc, t_info.ElemID, nElems, 0, t_info.i, t_info.j, t_info.k);

    if (t_info.ElemID <= nElems)
    {
        // Compute the split flux in the x-direction
        for(int l = 0; l <= Nloc; l++)
        {
            // Index flattened DVolSurf array
            // Must +1 as IndexFlatFortranArr assumes array is indexed from 1 on Fortran side but
            // DVolSurf is allocated as 0:PP_N
            VolSurf_offset = IndexFlatFortranArr(Nloc+1, l+1, t_info.i+1);
            // Find offsets that are local to this loop for the other arrays
            U_offset2 = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, l, t_info.j, t_info.k);
            UPrim_offset2 = FindPointOffset(PP_nVarPrim, Nloc, t_info.ElemID, l, t_info.j, t_info.k);
            Metric_offset2 = FindPointMetricsOffset(Nloc, t_info.ElemID, nElems, 0, l, t_info.j, t_info.k);
            
            for(int n = 0; n < PP_nVar; n++)
            {
                Ut_tmp[n] += DVolSurf[VolSurf_offset]*SplitDGVolume_Pointer_Device(n,&U[U_offset1], &UPrim[UPrim_offset1],
                                                                                   &U[U_offset2], &UPrim[UPrim_offset2],
                                                                                   &Mf[Metric_offset1], &Mf[Metric_offset2]
                                                                                  );
            }
        }

        // Repeat for the other two dimensions
        for(int l = 0; l <= Nloc; l++)
        {
            VolSurf_offset = IndexFlatFortranArr(Nloc+1, l+1, t_info.j+1);
            U_offset2 = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, t_info.i, l, t_info.k);
            UPrim_offset2 = FindPointOffset(PP_nVarPrim, Nloc, t_info.ElemID, t_info.i, l, t_info.k);
            Metric_offset2 = FindPointMetricsOffset(Nloc, t_info.ElemID, nElems, 0, t_info.i, l, t_info.k);
            
            for(int n = 0; n < PP_nVar; n++)
            {
                Ut_tmp[n] += DVolSurf[VolSurf_offset]*SplitDGVolume_Pointer_Device(n,&U[U_offset1], &UPrim[UPrim_offset1],
                                                                                   &U[U_offset2], &UPrim[UPrim_offset2],
                                                                                   &Mg[Metric_offset1], &Mg[Metric_offset2]
                                                                                  );
            }
        }

        for(int l = 0; l <= Nloc; l++)
        {
            VolSurf_offset = IndexFlatFortranArr(Nloc+1, l+1, t_info.k+1);
            U_offset2 = FindPointOffset(PP_nVar, Nloc, t_info.ElemID, t_info.i, t_info.j, l);
            UPrim_offset2 = FindPointOffset(PP_nVarPrim, Nloc, t_info.ElemID, t_info.i, t_info.j, l);
            Metric_offset2 = FindPointMetricsOffset(Nloc, t_info.ElemID, nElems, 0, t_info.i, t_info.j, l);
            
            for(int n = 0; n < PP_nVar; n++)
            {
                Ut_tmp[n] += DVolSurf[VolSurf_offset]*SplitDGVolume_Pointer_Device(n,&U[U_offset1], &UPrim[UPrim_offset1],
                                                                                   &U[U_offset2], &UPrim[UPrim_offset2],
                                                                                   &Mh[Metric_offset1], &Mh[Metric_offset2]
                                                                                  );
            }
        }

        // Perform final assignment to Ut array
        // U_offset1 can be reused to index Ut to current DOF for this thread
        for(int n = 0; n < PP_nVar; n++)
        {
#if PARABOLIC && !FV_ENABLED
            Ut[U_offset1 + n] += Ut_tmp[n];
#else
            Ut[U_offset1 + n] = Ut_tmp[n];
#endif
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------
// DEVICE Interface - Call kernel
//----------------------------------------------------------------------------------------------------------------------------------
void VolInt_splitForm_Device(int Nloc, int nElems, int d_Ut, int d_U, int d_UPrim, int d_Mf, int d_Mg,
                                 int d_Mh, int d_DVolSurf
                               , int streamID)
{
#ifdef USE_NVTX /* If profiling, add range */
    nvtxRangePushA("VolInt_splitForm");
#endif
    int nDOF = (Nloc+1)*(Nloc+1)*(Nloc+1)*nElems;

    INVOKE_KERNEL( VolInt_splitForm_Kernel, nDOF/256+1, 256, 0, streams[streamID], Nloc, nElems, 
                                                (double*)DeviceVars[d_Ut], (double*)DeviceVars[d_U],(double*)DeviceVars[d_UPrim],
                                                (double*)DeviceVars[d_Mf], (double*)DeviceVars[d_Mg], (double*)DeviceVars[d_Mh],
                                                (double*)DeviceVars[d_DVolSurf]
                 );
#ifdef USE_NVTX
    nvtxRangePop();
#endif
}
