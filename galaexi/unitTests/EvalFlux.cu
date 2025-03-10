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
#include "device.h"
#include "eos.h"
#include "../src/device/offsets.cu"
#include "../src/equations/navierstokes/flux.cu"

//----------------------------------------------------------------------------------------------------------------------------------
// External C prototypes for procedures that bind to Fortran
//----------------------------------------------------------------------------------------------------------------------------------
extern "C"
{
    void TestEvalEulerFlux1D_Device(int d_U_Long, int d_Flux);
}

__global__ void TestEvalEulerFlux1D_Kernel(double* U_Long, double* Flux)
{
    EvalEulerFlux1D_Device(U_Long, Flux);
}


void TestEvalEulerFlux1D_Device(int d_U_Long, int d_Flux)
{
    INVOKE_KERNEL( TestEvalEulerFlux1D_Kernel, 1, 1, 0, 0, (double*)DeviceVars[d_U_Long], (double*)DeviceVars[d_Flux] );
}