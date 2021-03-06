/** \file GPUCACell.h*/

#ifndef GPU_CACELL_H_
#define GPU_CACELL_H_

#include "GPUHitsAndDoublets.h"
#include "GPUSimpleVector.h"
#include <cuda.h>
#include <cmath>
#include <array>

//! The class that is used for representing the Cells-doublets
template<int numberOfLayers>
class GPUCACell
{
public:

	using CAntuplet = GPUSimpleVector<numberOfLayers, GPUCACell<numberOfLayers>>;
	__device__ GPUCACell()
	{

	}
	/** Initialize the cell to match the doublet info in the layer pair structures.
	\param doublets The layer pair structure
	\param layerId The layer of the doublet in the pairs
	\param doubletId The identifier of the doublet within the layer pair
	\param innerHitId The identifier of the inner hit in its layer
	\param outerHitId The identifier of the outer hit in its layer
	*/
	__device__
	void init(const GPULayerDoublets* doublets, int layerId, int doubletId,
			int innerHitId, int outerHitId)
	{

		theInnerHitId = innerHitId;
		theOuterHitId = outerHitId;

		theDoublets = doublets;

		theDoubletId = doubletId;
		theLayerIdInFourLayers = layerId;
		theInnerX = doublets->layers[0].x[doublets->indices[2 * doubletId]];
		theOuterX = doublets->layers[1].x[doublets->indices[2 * doubletId + 1]];

		theInnerY = doublets->layers[0].y[doublets->indices[2 * doubletId]];
		theOuterY = doublets->layers[1].y[doublets->indices[2 * doubletId + 1]];

		theInnerZ = doublets->layers[0].z[doublets->indices[2 * doubletId]];
		theOuterZ = doublets->layers[1].z[doublets->indices[2 * doubletId + 1]];
		theInnerR = hypot(theInnerX, theInnerY);
		theOuterR = hypot(theOuterX, theOuterY);
		theInnerNeighbors.reset();
	}

	__device__
	float get_inner_x() const
	{
		return theInnerX;
	}
	__device__
	float get_outer_x() const
	{
		return theOuterX;
	}
	__device__
	float get_inner_y() const
	{
		return theInnerY;
	}
	__device__
	float get_outer_y() const
	{
		return theOuterY;
	}
	__device__
	float get_inner_z() const
	{
		return theInnerZ;
	}
	__device__
	float get_outer_z() const
	{
		return theOuterZ;
	}
	__device__
	float get_inner_r() const
	{
		return theInnerR;
	}
	__device__
	float get_outer_r() const
	{
		return theOuterR;
	}
	__device__
	unsigned int get_inner_hit_id() const
	{
		return theInnerHitId;
	}
	__device__
	unsigned int get_outer_hit_id() const
	{
		return theOuterHitId;
	}

	__device__
	void print_cell() const
	{

		printf(
				"\n printing cell: %d, on layer: %d, innerHitId: %d, outerHitId: %d, innerradius %f, outerRadius %f ",
				theDoubletId, theLayerIdInFourLayers, theInnerHitId,
				theOuterHitId, theInnerR, theOuterR);

	}
	/**Check the compatibility with another cell
	\param innerCell The Cell checked for compatibility with "this" cell. 
	*/
	__device__
	bool check_alignment_and_tag(const GPUCACell<numberOfLayers>* innerCell,
			const float ptmin, const float region_origin_x,
			const float region_origin_y, const float region_origin_radius,
			const float thetaCut, const float phiCut)
	{

		return (are_aligned_RZ(innerCell, ptmin, thetaCut)
				&& have_similar_curvature(innerCell, region_origin_x,
						region_origin_y, region_origin_radius, phiCut));

	}

	__device__
	bool are_aligned_RZ(const GPUCACell<numberOfLayers>* otherCell,
			const float ptmin, const float thetaCut) const
	{

		float r1 = otherCell->get_inner_r();
		float z1 = otherCell->get_inner_z();
		float distance_13_squared = (r1 - theOuterR) * (r1 - theOuterR)
				+ (z1 - theOuterZ) * (z1 - theOuterZ);
		float tan_12_13_half = fabs(
				z1 * (theInnerR - theOuterR) + theInnerZ * (theOuterR - r1)
						+ theOuterZ * (r1 - theInnerR)) / distance_13_squared;
		return tan_12_13_half * ptmin <= thetaCut;
	}

	__device__
	bool have_similar_curvature(const GPUCACell<numberOfLayers>* otherCell,
			const float region_origin_x, const float region_origin_y,
			const float region_origin_radius, const float phiCut) const
	{
		auto x1 = otherCell->get_inner_x();
		auto y1 = otherCell->get_inner_y();

		auto x2 = get_inner_x();
		auto y2 = get_inner_y();

		auto x3 = get_outer_x();
		auto y3 = get_outer_y();

		auto precision = 0.5f;
		auto offset = x2 * x2 + y2 * y2;

		auto bc = (x1 * x1 + y1 * y1 - offset) / 2.f;

		auto cd = (offset - x3 * x3 - y3 * y3) / 2.f;

		auto det = (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);

		//points are aligned
		if (fabs(det) < precision)
			return true;

		auto idet = 1.f / det;

		auto x_center = (bc * (y2 - y3) - cd * (y1 - y2)) * idet;
		auto y_center = (cd * (x1 - x2) - bc * (x2 - x3)) * idet;

		auto radius = std::sqrt(
				(x2 - x_center) * (x2 - x_center)
						+ (y2 - y_center) * (y2 - y_center));
		auto centers_distance_squared = (x_center - region_origin_x)
				* (x_center - region_origin_x)
				+ (y_center - region_origin_y) * (y_center - region_origin_y);

		auto minimumOfIntesectionRange = (radius - region_origin_radius)
				* (radius - region_origin_radius) - phiCut;

		if (centers_distance_squared >= minimumOfIntesectionRange)
		{
			auto maximumOfIntesectionRange = (radius + region_origin_radius)
					* (radius + region_origin_radius) + phiCut;
			return centers_distance_squared <= maximumOfIntesectionRange;
		}
		else
		{

			return false;
		}

	}

	//! Retrieves quadruplets out of connected doublets
	/*! This function gets the connected doublets and performs a DFS to produce the quadruplets.
	\param foundNtuplets The result returned for our problem
	\param tmpNtuplet The current path of the search
	\param minHitsPerNtuplet The number of layers involved
	*/
	template<int maxNumberOfQuadruplets>
	__device__
	void find_ntuplets(
			GPUSimpleVector<maxNumberOfQuadruplets,int4>* foundNtuplets,
			GPUSimpleVector<4, GPUCACell<4>*>& tmpNtuplet,
			const unsigned int minHitsPerNtuplet
	) const
	{
		GPUCACell<numberOfLayers>* otherCell;
		int4 found;

		if (theInnerNeighbors.size() == 0)
		{
			if (tmpNtuplet.size() >= minHitsPerNtuplet - 1)
			{
				found.x=tmpNtuplet.m_data[2]->get_inner_hit_id();
				found.y=tmpNtuplet.m_data[2]->get_outer_hit_id();
				found.z=tmpNtuplet.m_data[1]->get_outer_hit_id();
				found.w=tmpNtuplet.m_data[0]->get_outer_hit_id();
				foundNtuplets->push_back_ts(found);
			}
			else
			return;
		}
		else
		{
			for (int j = 0; j < theInnerNeighbors.size(); ++j)
			{
				otherCell = theInnerNeighbors.m_data[j];
				tmpNtuplet.push_back(otherCell);
				otherCell->find_ntuplets(foundNtuplets, tmpNtuplet, minHitsPerNtuplet);
				tmpNtuplet.pop_back();

			}

		}
	}
	GPUSimpleVector<64, GPUCACell<numberOfLayers>*> theInnerNeighbors;

private:

	unsigned int theInnerHitId;
	unsigned int theOuterHitId;
	const GPULayerDoublets* theDoublets;
	int theDoubletId;
	int theLayerIdInFourLayers;
	float theInnerX;
	float theOuterX;
	float theInnerY;
	float theOuterY;
	float theInnerZ;
	float theOuterZ;
	float theInnerR;
	float theOuterR;

};

#endif /*CACELL_H_ */
