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

/**
 * @brief This file includes all CUDA files from the entire project in one place.
 * We do this because of how the CUDA compiler works. It optimizes code the best when 
 * it all appears in the same compute unit.
 * 
 * Early CUDA did not have a linker so this type of file used to be required for CUDA codes.
 * Now (08.2024) CUDA does have a linker (not sure about HIP), but using it (-rdc=true)
 * can, and usually does, cause significant losses in performance.
 */

#include "flexi.h"
#include "device.h"
#include "eos.h"

// WARNING!!!: The following list of includes is ordered as such for dependency reasons
//             DO NOT CHANGE THE ORDER unless you know what you are doing.
#include "device_manage.cu"
#include "offsets.cu"
#include "../mesh/mesh.cu"
#include "../globals/vector.cu"
#include "../equations/navierstokes/idealgas/eos.cu"
#include "../equations/navierstokes/calctimestep.cu"
#include "../interpolation/applyjacobian.cu"
#include "../interpolation/prolongtoface.cu"
#include "../equations/navierstokes/splitflux.cu"
#include "../dg/surfint.cu"
#include "../equations/navierstokes/flux.cu"
#include "../dg/applydmatrix.cu"
#include "../dg/volint.cu"
#include "../equations/navierstokes/riemann.cu"
#include "../dg/fillflux.cu"