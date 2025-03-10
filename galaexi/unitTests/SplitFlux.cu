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

//----------------------------------------------------------------------------------------------------------------------------------
// External C prototypes for procedures that bind to Fortran
//----------------------------------------------------------------------------------------------------------------------------------
extern "C"
{
    void TestVolumeMethods_Device(int splitDGOpt, int d_U, int d_UPrim, int d_U_master, int d_UPrim_master, int d_Metrics_fTilde, int d_Ja_Face, int d_Flux);
    void TestSurfaceMethods_Device(int splitDGOpt, int d_U_LL, int d_U_RR, int d_Flux);
}

__global__ void TestVolumeMethods_Kernel(int splitDGOpt, double* Uref, double* UPrimRef, double* U, double* UPrim, double* Mref, double* M, double* Flux)
{
    double Ja_face[3];

    Ja_face[0] = M[0];
    Ja_face[1] = M[3];
    Ja_face[2] = M[6];

    // For some reason that I (Spencer Starr -- 08.2024) cannot figure out, if we rely on the InitSplitDG_Device
    // method from splitflux.cu to set the SplitDG function pointer for this test, it will sometimes be NULL.
    // Whether it is NULL or not seems to be completely random. I believe it has something to do with scope of the
    // pointers themselves (library vs. this test). Seeing as this is a unit test, bending how the function pointers work
    // (which we know to be working correctly in "normal" operation as that is implicitly tested in the SurfInt test)
    // to pass here is silly. Instead this workaround just resets the device pointer here in place, which removes
    // any question whether the pointer is NULL or not.
    // There is a repeat of this issue in the Riemann solver test.
    switch (splitDGOpt)
    {
        case PRM_SPLITDG_SD:
            SplitDGVolume_Pointer_Device  = SplitVolumeFluxSD_Device;
            break;
        case PRM_SPLITDG_MO:
            SplitDGVolume_Pointer_Device  = SplitVolumeFluxMO_Device;
            break;
        case PRM_SPLITDG_DU:
            SplitDGVolume_Pointer_Device  = SplitVolumeFluxDU_Device;
            break;
        case PRM_SPLITDG_KG:
            SplitDGVolume_Pointer_Device  = SplitVolumeFluxKG_Device;
            break;
        case PRM_SPLITDG_PI:
            SplitDGVolume_Pointer_Device  = SplitVolumeFluxPI_Device;
            break;
        case PRM_SPLITDG_CH:
            SplitDGVolume_Pointer_Device  = SplitVolumeFluxCH_Device;
            break;
    }
    
    for (int n = 0; n < PP_nVar; n ++)
    {
        Flux[n] = SplitDGVolume_Pointer_Device( n, Uref, UPrimRef, U, UPrim, Mref, &Ja_face[0] );
    }    
}

__global__ void TestSurfaceMethods_Kernel(int splitDGOpt, double* U_LL, double* U_RR, double* Flux)
{

    // See comment in above in volume method test for why this is here.
    switch (splitDGOpt)
    {
        case PRM_SPLITDG_SD:
            SplitDGSurface_Pointer_Device = SplitSurfaceFluxSD_Device;
            break;
        case PRM_SPLITDG_MO:
            SplitDGSurface_Pointer_Device = SplitSurfaceFluxMO_Device;
            break;
        case PRM_SPLITDG_DU:
            SplitDGSurface_Pointer_Device = SplitSurfaceFluxDU_Device;
            break;
        case PRM_SPLITDG_KG:
            SplitDGSurface_Pointer_Device = SplitSurfaceFluxKG_Device;
            break;
        case PRM_SPLITDG_PI:
            SplitDGSurface_Pointer_Device = SplitSurfaceFluxPI_Device;
            break;
        case PRM_SPLITDG_CH:
            SplitDGSurface_Pointer_Device = SplitSurfaceFluxCH_Device;
            break;
    }

    for (int n = 0; n < PP_nVar; n++)
    {
        Flux[n] = SplitDGSurface_Pointer_Device( n, U_LL, U_RR );
    }    
}


void TestVolumeMethods_Device(int splitDGOpt, int d_U, int d_UPrim, int d_U_master, int d_UPrim_master, int d_Metrics_fTilde, int d_Ja_Face, int d_Flux)
{
    INVOKE_KERNEL( TestVolumeMethods_Kernel, 1, 1, 0, 0, splitDGOpt, (double*)DeviceVars[d_U], (double*)DeviceVars[d_UPrim], (double*)DeviceVars[d_U_master], 
                                    (double*)DeviceVars[d_UPrim_master], (double*)DeviceVars[d_Metrics_fTilde],
                                    (double*)DeviceVars[d_Ja_Face], (double*)DeviceVars[d_Flux] );
}

void TestSurfaceMethods_Device(int splitDGOpt, int d_U_LL, int d_U_RR, int d_Flux)
{
    INVOKE_KERNEL( TestSurfaceMethods_Kernel, 1, 1, 0, 0, splitDGOpt, (double*)DeviceVars[d_U_LL], (double*)DeviceVars[d_U_RR], (double*)DeviceVars[d_Flux] );
}
