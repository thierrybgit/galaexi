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
#include "device.h"
#include "eos.h"
#include "../src/device/offsets.cu"
#include "../src/equations/navierstokes/idealgas/eos.cu"
#include "../src/equations/navierstokes/splitflux.cu"
#include "../equations/navierstokes/flux.cu"
#include "../src/equations/navierstokes/riemann.cu"

//----------------------------------------------------------------------------------------------------------------------------------
// External C prototypes for procedures that bind to Fortran
//----------------------------------------------------------------------------------------------------------------------------------
extern "C"
{
    void TestMethods_Device(int riemannID, double kappa, int d_U_L, int d_U_R, int d_UPrim_L, int d_UPrim_R, int d_nv, 
                                                                                int d_t1, int d_t2, int d_Flux);
}

__global__ void TestMethods_Kernel(int riemannID, double kappa, double* U_L, double* U_R, double* UPrim_L, double* UPrim_R, 
                                                                    double* nv, double* t1, double* t2, double* Flux)
{   
    d_Kappa = kappa;
    d_KappaM1 = d_Kappa - 1.0;

    // For some reason that I (Spencer Starr -- 08.2024) cannot figure out, if we rely on the InitSplitDG_Device
    // and InitRiemann_Device methods to set the function pointers for this test, they will sometimes be NULL. 
    // Whether they is NULL or not seems to be completely random. I believe it has something to do with scope of the
    // pointers themselves (library vs. this test). Seeing as this is a unit test, bending how the function pointers work
    // (which we know to be working correctly in "normal" operation as that is implicitly tested in the SurfInt test)
    // to pass here is silly. Instead this workaround just resets the device pointer here in place, which removes
    // any question whether the pointer is NULL or not.
    // There is a similar workaround in the unit test for SplitFlux.
#ifdef SPLIT_DG
    // We aren't testing the splitflux methods here, so just use the default
    SplitDGSurface_Pointer_Device = SplitSurfaceFluxSD_Device;
#endif

    switch (riemannID)
    {
        case PRM_RIEMANN_LF:
            Riemann_Pointer_Device = Riemann_LF_Device;
            break;
        case PRM_RIEMANN_ROE:
            Riemann_Pointer_Device = Riemann_Roe_Device;
            break;
        case PRM_RIEMANN_ROEENTROPYFIX:
            Riemann_Pointer_Device = Riemann_RoeEntropyFix_Device;
            break;
        case PRM_RIEMANN_ROEL2:
            Riemann_Pointer_Device = Riemann_RoeL2_Device;
            break;
#ifndef SPLIT_DG
        case PRM_RIEMANN_HLL:
            Riemann_Pointer_Device = Riemann_HLL_Device;
            break;
        case PRM_RIEMANN_HLLC:
            Riemann_Pointer_Device = Riemann_HLLC_Device;
            break;
        case PRM_RIEMANN_HLLE:
            Riemann_Pointer_Device = Riemann_HLLE_Device;
            break;
        case PRM_RIEMANN_HLLEM:
            Riemann_Pointer_Device = Riemann_HLLEM_Device;
            break;
#else /* SPLIT_DG */
        case PRM_RIEMANN_CH:
            Riemann_Pointer_Device = Riemann_CH_Device;
            break;
        case PRM_RIEMANN_Average:
            Riemann_Pointer_Device = Riemann_FluxAverage_Device;
            break;
#endif /* SPLIT_DG */
    }

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Riemann_Point_Device(Flux, U_L, U_R, UPrim_L, UPrim_R, nv, t1, t2, false);
}


void TestMethods_Device(int riemannID, double kappa, int d_U_L, int d_U_R, int d_UPrim_L, int d_UPrim_R, int d_nv, 
                                                                                int d_t1, int d_t2, int d_Flux)
{
    INVOKE_KERNEL( TestMethods_Kernel, 1, 1, 0, 0, riemannID, kappa, (double*)DeviceVars[d_U_L], (double*)DeviceVars[d_U_R], 
                                    (double*)DeviceVars[d_UPrim_L], (double*)DeviceVars[d_UPrim_R], (double*)DeviceVars[d_nv],
                                    (double*)DeviceVars[d_t1], (double*)DeviceVars[d_t2], (double*)DeviceVars[d_Flux]  );
}