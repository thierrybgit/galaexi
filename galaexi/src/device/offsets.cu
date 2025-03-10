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

/***********************************************************************************************************************
 * @file offsets.cu
 * @brief Functions that wrap the logic of indexing flattened Fortran arrays
***********************************************************************************************************************/

/***********************************************************************************************************************************
 * The ThreadIndices structure holds index information about the specific DOF assigned to a device thread.
 * It hides the logic for calculating these IDs and indices with the method FindThreadIndices.
 * 
 * Nominially this is a class with FindThreadIndices as a member method. We do it this way because classes CANNOT be used in kernels.
 **********************************************************************************************************************************/
struct ThreadIndicesVolume
{
    int threadID;    /// ID of the current device thread. Indexed from 1.
    int ElemID;      /// ID of the element that has the DOF assigned to this thread. Indexed from 1.
    int i, j, k;     /// Indices within element ElemID for the DOF assigned to thread. Indexed from 0.
};

struct ThreadIndicesSide
{
    int threadID;    /// ID of the current device thread. Indexed from 1.
    int SideID;      /// ID of the side that has the DOF assigned to this thread. Indexed from 1.
    int p,q;         /// Indices within side SideID for the DOF assigned to thread. Indexed from 0.
};

//----------------------------------------------------------------------------------------------------------------------------------
// Device thread indexing methods
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Wraps the often used logic of finding the element and DOF indices for a device thread.
 * @returns An instance of the ThreadIndices class initialized with the element ID and DOF indices
 *          of the calling device thread.
 * @param Nloc Local polynomial order
 * @remark WARNING! Indexes threadIDs from 1, not 0. This is due to the logic for finding i,j,k.
 **********************************************************************************************************************************/
__device__ ThreadIndicesVolume FindThreadIndicesVolume(int Nloc)
{

    struct ThreadIndicesVolume t_info;
    int rest = 0;

    t_info.threadID = blockIdx.x * blockDim.x + threadIdx.x + 1;
    t_info.ElemID = ((t_info.threadID-1) / ((Nloc+1)*(Nloc+1)*(Nloc+1))) + 1;
    rest = t_info.threadID - (t_info.ElemID-1)*(Nloc+1)*(Nloc+1)*(Nloc+1);
    t_info.k = (rest-1)/((Nloc+1)*(Nloc+1));
    rest = rest - t_info.k*(Nloc+1)*(Nloc+1);
    t_info.j = (rest-1)/(Nloc+1);
    rest = rest - t_info.j*(Nloc+1);
    t_info.i = rest - 1;

    return t_info;
}

/***********************************************************************************************************************************
 * @brief Wraps the often used logic of finding the side and DOF indices for a device thread.
 * @returns An instance of the ThreadIndices class initialized with the side ID and DOF indices
 *          of the calling device thread.
 * @param Nloc Local polynomial order
 * @param startSide For kernels that don't work on a full sides array, this allows offseting SideID.
 * @remark WARNING! Indexes threadIDs from 1, not 0. This is due to the logic for finding p and q.
 **********************************************************************************************************************************/
__device__ ThreadIndicesSide FindThreadIndicesSide(int Nloc, int startSide)
{

    struct ThreadIndicesSide t_info;
    int rest = 0;

    t_info.threadID = blockIdx.x * blockDim.x + threadIdx.x + 1;
    t_info.SideID = ((t_info.threadID-1) / ((Nloc+1)*(Nloc+1))) + 1;
    rest = t_info.threadID - (t_info.SideID-1)*(Nloc+1)*(Nloc+1);
    t_info.q = (rest-1)/(Nloc+1);
    rest = rest - t_info.q*(Nloc+1);
    t_info.p = rest - 1;

    t_info.SideID += (startSide-1);

    return t_info;
}


//----------------------------------------------------------------------------------------------------------------------------------
// Offset flattened memory methods -- ELEMENT
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Find the offset for a specific element in 1D flattened device mem
 * @param nVar Number of flow variables
 * @param Nloc Local polynomial order
 * @param elemIdx Index of the desired element (assumed to be indexed from 1)
 * @returns Computed offset in bytes to the starting index of the desired element
 **********************************************************************************************************************************/
__device__ __host__ int FindElemOffset(int nVar, int Nloc, int elemIdx)
{
    return (Nloc+1)*(Nloc+1)*(Nloc+1)*nVar*(elemIdx-1);
}

/***********************************************************************************************************************************
 * @brief Find the offset of an element in the Metrics arrays
 * @param Nloc Local polynomial order
 * @param elemIdx Index of the current element (assumed to be indexed from 1)
 * @param FV_idx Idx for FV_SIZE  (assumed to be indexed from 0)
 **********************************************************************************************************************************/
__device__ __host__ int FindElemMetricsOffset(int Nloc, int elemIdx, int FV_idx)
{
    return (Nloc+1)*(Nloc+1)*(Nloc+1)*3*(elemIdx-1)*FV_idx + (Nloc+1)*(Nloc+1)*(Nloc+1)*3*(elemIdx-1);
}

/***********************************************************************************************************************************
 * @brief Find the offset of an element in the Jacobian array
 * @param Nloc Local polynomial order
 * @param elemIdx Index of the current element (assumed to be indexed from 1)
 * @param FV_idx Idx for FV_SIZE  (assumed to be indexed from 0)
 **********************************************************************************************************************************/
__device__ __host__ int FindElemJacobianOffset(int Nloc, int elemIdx, int FV_idx)
{
    return (Nloc+1)*(Nloc+1)*(Nloc+1)*(elemIdx-1)*FV_idx + (Nloc+1)*(Nloc+1)*(Nloc+1)*(elemIdx-1);
}


//----------------------------------------------------------------------------------------------------------------------------------
// Offset flattened memory methods -- SIDE
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Find the offset for a specific side of an element in 1D flattened device mem
 * @param nVar Number of flow variables
 * @param Nloc Local polynomial order
 * @param elemIdx Index of the desired element  (assumed to be indexed from 1)
 * @param sideIdx Index of the side. (assumed to be indexed from 0)
 * @returns Computed offset in bytes to the starting index of the side within the desired element
 * @remark This method assumes that the sideID is indexed from zero.
 **********************************************************************************************************************************/
__device__ __host__ int FindSideOffset(int nVar, int Nloc, int elemIdx, int sideIdx)
{
    return (Nloc+1)*(Nloc+1)*(Nloc+1)*nVar*(elemIdx-1) + (Nloc+1)*(Nloc+1)*nVar*sideIdx;
}

/***********************************************************************************************************************************
 * @brief Find the offset to the start of a side in a master/slave sides array
 * @param nVar Number of variables at each DOF
 * @param Nloc Local polynomial order
 * @param sideIdx Index of the desired side (assumed to be indexed from 0)
 * @returns Computed offset in bytes to the starting index of the desired side
 **********************************************************************************************************************************/
__device__ __host__ int FindSideMasterSlaveOffset(int nVar, int Nloc, int sideIdx)
{
    return (Nloc+1)*(Nloc+1)*nVar*sideIdx;
}

/***********************************************************************************************************************************
 * @brief Find the offset to the start of a side in the surface vector arrays
 * @param nVars Number of vars in the 1st dim in the array. Used to allow this method to offset arrays like NormVec (nDim=3) and 
 *             SurfElem (nDim = 1)
 * @param Nloc Local polynomial order
 * @param sideIdx Index of the side to offset to (assumed to be indexed from 0)
 * @returns Computed offset in bytes to the starting index of the desired side
 **********************************************************************************************************************************/
__device__ __host__ int FindSideSurfaceDataOffset(int nVars, int Nloc, int nSides, int FV_idx, int sideIdx)
{
    // (3,0:PP_N,0:PP_N,1:nSides,0:FV_SIZE)
    return FV_idx*nSides*(Nloc+1)*(Nloc+1)*nVars + sideIdx*(Nloc+1)*(Nloc+1)*nVars;
}


//----------------------------------------------------------------------------------------------------------------------------------
// Offset flattened memory methods -- POINT/DOF
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Find the offset for a specific DOF in 1D flattened device mem
 * @param nVar Number of flow variables
 * @param Nloc Local polynomial order
 * @param elemIdx Index of the desired element (Indexed from 1)
 * @param i Index of the first dimension DOF of the element. (Indexed from 0)
 * @param j Index of the second dimension DOF of the element. (Indexed from 0)
 * @param k Index of the third dimension DOF of the element. (Indexed from 0)
 * @returns Computed offset in bytes to the starting index of desired single DOF
 **********************************************************************************************************************************/
__device__ __host__ int FindPointOffset(int nVar, int Nloc, int elemIdx, int i, int j, int k)
{
    return (Nloc+1)*(Nloc+1)*(Nloc+1)*nVar*(elemIdx-1) + (Nloc+1)*(Nloc+1)*nVar*k + (Nloc+1)*nVar*j + nVar*i;
}

/***********************************************************************************************************************************
 * @brief Find the offset of a point in a flattened master/slave sides array
 * @param nVar Number of variables at each DOF
 * @param Nloc Local polynomial order
 * @param sideIdx Index of the desired side (assumed to be indexed from 0)
 * @param i Index of the first dimension DOF of the element. (assumed to be indexed from 0)
 * @param j Index of the second dimension DOF of the element. (assumed to be indexed from 0)
 * @returns Computed offset in bytes to the starting index of the desired DOF
 **********************************************************************************************************************************/
__device__ __host__ int FindPointMasterSlaveOffset(int nVar, int Nloc, int sideIdx, int i, int j)
{
    return (Nloc+1)*(Nloc+1)*nVar*sideIdx + (Nloc+1)*nVar*j + nVar*i;
}

/***********************************************************************************************************************************
 * @brief Find the offset of a point in a flattened metrics array
 * @param Nloc Local polynomial order
 * @param elemIdx Index of the desired element (assumed to be indexed from 1)
 * @param FV_idx Idx for FV_SIZE  (assumed to be indexed from 0)
 * @param i Index of the first dimension DOF of the element. (assumed to be indexed from 0)
 * @param j Index of the second dimension DOF of the element. (assumed to be indexed from 0)
 * @param k Index of the third dimension DOF of the element. (assumed to be indexed from 0)
 * @returns Computed offset in bytes to the starting index of the desired DOF
 **********************************************************************************************************************************/
__device__ __host__ int FindPointMetricsOffset(int Nloc, int elemIdx, int nElems, int FV_idx, int i, int j, int k)
{
    return (Nloc+1)*(Nloc+1)*(Nloc+1)*3*nElems*FV_idx + (Nloc+1)*(Nloc+1)*(Nloc+1)*3*(elemIdx-1) + (Nloc+1)*(Nloc+1)*3*k + (Nloc+1)*3*j + 3*i;
}

/***********************************************************************************************************************************
 * @brief Find the offset of a point in a flattened Jacobian array
 * @param Nloc Local polynomial order
 * @param elemIdx Index of the desired element (assumed to be indexed from 1)
 * @param FV_idx Idx for FV_SIZE (assumed to be indexed from 0)
 * @param i Index of the first dimension DOF of the element. (assumed to be indexed from 0)
 * @param j Index of the second dimension DOF of the element. (assumed to be indexed from 0)
 * @param k Index of the third dimension DOF of the element. (assumed to be indexed from 0)
 * @returns Computed offset in bytes to the starting index of the desired DOF
 **********************************************************************************************************************************/
__device__ __host__ int FindPointJacobianOffset(int Nloc, int elemIdx, int nElems, int FV_idx, int i, int j, int k)
{
    return (Nloc+1)*(Nloc+1)*(Nloc+1)*nElems*FV_idx + (Nloc+1)*(Nloc+1)*(Nloc+1)*(elemIdx-1) + (Nloc+1)*(Nloc+1)*k + (Nloc+1)*j + i;
}

/***********************************************************************************************************************************
 * @brief Find the offset of a point in one flattened mapping arrays (S2V, V2S, S2V2, etc.)
 * @param sideID Index of the current side of the desired element. (assumed to be indexed from 0)
 * @param Nloc Local polynomial order 
 * @param flip Index of flip side (assumed to be indexed from 0)
 * @param i Index of the first dimension DOF of the element. (assumed to be indexed from 0)
 * @param j Index of the second dimension DOF of the element. (assumed to be indexed from 0)
 * @param k Index of the third dimension DOF of the element. (assumed to be indexed from 0)
 * @param isMapping2 Is the mapping array being indexed S2V2 (true) or S2V (false)?
 **********************************************************************************************************************************/
__device__ __host__ int FindPointMappingOffset(int sideID, int Nloc, int i, int j, int k, int flip, bool isMapping2)
{
    if (! isMapping2)
    {
        // S2V(3,0:Nloc,0:Nloc,0:Nloc,0:4,1:6)
        return sideID*5*(Nloc+1)*(Nloc+1)*(Nloc+1)*3 + flip*(Nloc+1)*(Nloc+1)*(Nloc+1)*3 + k*(Nloc+1)*(Nloc+1)*3 + j*(Nloc+1)*3 + i*3;
    }
    else
    {
        // S2V2(2,0:Nloc,0:Nloc,0:4,1:6)
        return sideID*5*(Nloc+1)*(Nloc+1)*2 + flip*(Nloc+1)*(Nloc+1)*2 + k*(Nloc+1)*2 + j*2;
    }
    
}


//----------------------------------------------------------------------------------------------------------------------------------
// Offset flattened memory methods -- STANDARD FORTRAN ARRAYS
//----------------------------------------------------------------------------------------------------------------------------------
/***********************************************************************************************************************************
 * @brief Overload for 2D memory that is allocated from 1 in both dimensions
 * @param xSize Size of the 1st dimension of the array
 * @param xIdx Fortran index of the 1st dimension of the array you would like the offset for (assumed to be indexed from 1)
 * @param yIdx Fortran index of the 2nd dimension of the array you would like the offset for (assumed to be indexed from 1)
 * @returns Calculated offset of the desired element in flattened memory
 * @remark ASSUMES THAT xIdx AND yIdx ARE INDEXED STARTING AT 1
 **********************************************************************************************************************************/
__device__ __host__ int IndexFlatFortranArr(int xSize, int xIdx, int yIdx)
{
    return (yIdx-1)*xSize + (xIdx-1);
}

/***********************************************************************************************************************************
 * @brief Overload for 3D memory that is allocated from 1 in all dimensions
 * @param xSize Size of the 1st dimension of the array
 * @param ySize Size of the 2nd dimension of the array
 * @param xIdx Fortran index of the 1st dimension of the array you would like the offset for (assumed to be indexed from 1)
 * @param yIdx Fortran index of the 2nd dimension of the array you would like the offset for (assumed to be indexed from 1)
 * @param zIdx Fortran index of the 3rd dimension of the array you would like the offset for (assumed to be indexed from 1)
 * @returns Calculated offset of the desired element in flattened memory
 * @remark ASSUMES THAT xIdx, yIdx AND zIdx ARE INDEXED STARTING AT 1
 **********************************************************************************************************************************/
__device__ __host__ int IndexFlatFortranArr(int xSize, int ySize, int xIdx, int yIdx, int zIdx)
{
    return ((zIdx-1)*ySize*xSize) + ((yIdx-1)*xSize) + (xIdx-1);
}


// THIS METHOD WAS USED FOR NORMVEC, SURFELEM, ETC. WHEN THEY HAD THE OLD ORDER (FV_SIZE BEFORE NSIDES).
// THE ALLOCATED ORDER OF THOSE DIMENSIONS WERE SWITCHED TO MAKE OFFSETTING THOSE ARRAYS EASIER.
// THAT MADE THIS METHOD UNNECESSARY, BUT IT IS KEPT HERE JUST IN CASE A SIMILAR SITUATION IS NEEDED LATER OR 
// THOSE ARRAYS HAVE TO BE CHANGED BACK TO THEIR ORIGINAL ALLOCATION ORDER.
// /***********************************************************************************************************************************
//  * @brief Find offset for a device thread to a DOF accounting for FV
//  * 
//  * Many of the arrays with surface information a dimension that accounts for neighboring FV elements. When offsetting to a DOF in
//  * those arrays, we have to use the correct "FV side". This method finds the offset in bytes within those array accounting for
//  * the "FV side" choice.
//  * 
//  * Main use is for offsetting NormVec, TangVec1 and TangVec2, but can be used for other similarly structure arrays.
//  * Long term, it would be better to restructure these arrays in the Fortran code to have the FV_SIZE dim come last.
//  * 
//  * @param threadID ID of the thread we are finding the offset for
//  * @param Nloc Local polynomial order
//  * @param startSideID Index-0 side ID for the first side of the array slice currently being worked on
//  * @param initOffset Offset in bytes to the start of side startSideID
//  * @param FV_Elems_Max Flag array storing which "FV side" to use for each element
//  * @param nVars Number of vars in the 1st dim in the array. Used to allow this method to offset arrays like NormVec (nDim=3) and 
//  *             SurfElem (nDim = 1)
//  * @returns Computed offset in bytes for the DOF within the correct "FV side"
//  * @remark When FV is ported, will need to re-add passing of FV_Elems_Max in second term of return
//  **********************************************************************************************************************************/
// __device__ int FindSurfaceDataOffsetWithFV(int threadID, int Nloc, int startSideID, int initOffset, int nVars)
// {
//     int nDOFsSide = (Nloc+1)*(Nloc+1);
//     int thisDOF_SideOffset = floor((double)threadID / (double)(nDOFsSide));
//     int thisDOF_SideID = startSideID + thisDOF_SideOffset;
//                  start elem + offset to this DOF's side + offset to FV side + offset to DOF within FV side
//     return initOffset + thisDOF_SideOffset*(FV_SIZE+1)*nDOFsSide*nVars + 
//                                     FV_Elems_Max[thisDOF_SideID]*nDOFsSide*nVars + thisDOF_SideOffset*nVars;
// }