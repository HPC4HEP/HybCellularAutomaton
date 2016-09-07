/** \file ca.h*/


#include <vector>
#include <array>

#include <stdio.h>


#include "GPUCACell.h"

/** This kernel initializes the cells and the hit indexes.
*/
template<int numberOfLayers>
__global__
void kernel_create(const GPULayerDoublets* gpuDoublets,
		GPUCACell<numberOfLayers>** cells,
                GPUSimpleVector<64, GPUCACell<numberOfLayers>* > ** isOuterHitOfCell)
{

	unsigned int layerPairIndex = blockIdx.y;
	unsigned int cellIndexInLayerPair = threadIdx.x + blockIdx.x * blockDim.x;
	if(layerPairIndex < numberOfLayers-1)
	{

		for(int i = cellIndexInLayerPair; i < gpuDoublets[layerPairIndex].size; i+=gridDim.x * blockDim.x)
		{

			cells[layerPairIndex][i].init(&gpuDoublets[layerPairIndex],layerPairIndex,i,gpuDoublets[layerPairIndex].indices[2*i], gpuDoublets[layerPairIndex].indices[2*i+1]);

			if(layerPairIndex < 2){
				isOuterHitOfCell[layerPairIndex][cells[layerPairIndex][i].get_outer_hit_id()].push_back_ts(& (cells[layerPairIndex][i]));

			}
		}
	}

}

/** This kernel connects compatible cells of adjacent layer pairs.
*/
template<int numberOfLayers>
__global__
void kernel_connect(const GPULayerDoublets* gpuDoublets,
		GPUCACell<numberOfLayers>** cells,
                GPUSimpleVector<64, GPUCACell<numberOfLayers>* > ** isOuterHitOfCell,
		float ptmin,
		float region_origin_x,
		float region_origin_y,
		float region_origin_radius,
		float thetaCut,
		float phiCut)
{
	unsigned int layerPairIndex = blockIdx.y+1;
	unsigned int cellIndexInLayerPair = threadIdx.x + blockIdx.x * blockDim.x;
	if(layerPairIndex < numberOfLayers-1)
	{
		for (int i = cellIndexInLayerPair; i < gpuDoublets[layerPairIndex].size; i += gridDim.x * blockDim.x)
		{

           for (int j = 0; j < isOuterHitOfCell[layerPairIndex-1][cells[layerPairIndex][i].get_inner_hit_id()].size(); ++j)
			{
        	   GPUCACell<numberOfLayers>* otherCell= isOuterHitOfCell[layerPairIndex-1][cells[layerPairIndex][i].get_inner_hit_id()].m_data[j];

        	   if (cells[layerPairIndex][i].check_alignment_and_tag(otherCell,
								ptmin, region_origin_x, region_origin_y,
								region_origin_radius, thetaCut, phiCut))
        	   {
        		   	   cells[layerPairIndex][i].theInnerNeighbors.push_back_ts(otherCell);
        	   }
			}
		}
	}
}

/** This kernel retrieves quadruplets by DFS over the connections
*/
template<int numberOfLayers, int maxNumberOfQuadruplets>
__global__
void kernel_find_ntuplets(const GPULayerDoublets* gpuDoublets,
		GPUCACell<numberOfLayers>** cells,
		GPUSimpleVector<maxNumberOfQuadruplets, int4>* foundNtuplets,
		unsigned int minHitsPerNtuplet)
{

	unsigned int cellIndexInLastLayerPair = threadIdx.x + blockIdx.x * blockDim.x;
	constexpr unsigned int lastLayerPairIndex = numberOfLayers - 2;


	GPUSimpleVector<4, GPUCACell<4>*> stack;
	for (int i = cellIndexInLastLayerPair; i < gpuDoublets[lastLayerPairIndex].size;
			i += gridDim.x * blockDim.x)
	{


		stack.reset();
		stack.push_back(&cells[lastLayerPairIndex][i]);
		cells[lastLayerPairIndex][i].find_ntuplets(foundNtuplets, stack, minHitsPerNtuplet);


	}


}

