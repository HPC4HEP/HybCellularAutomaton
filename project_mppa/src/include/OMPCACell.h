/** \file OMPCACell.h*/


#ifndef OMP_CACELL_H_
#define OMP_CACELL_H_
 
#include "OMPHitsAndDoublets.h"
#include "OMPSimpleVector.h"
#include "OMPResultVector.h"
#include "OMPMixedList.h"

#include <math.h>

const int numberOfLayers = 4;

//! The Cellular Automaton cell
/** Not used here because it creates duplicates of values. Wastes precious memory
*/
typedef struct OMPCACell
{	
	int theInnerHitId; /**< integer index of inner hit in its layer*/
	int theOuterHitId; /**< integer index of outer hit in its layer*/
	int theDoubletId; /**< integer index of doublets in the layer pair structure*/
	int theLayerIdInFourLayers; /**< the id of the layer*/
	float theInnerX; /**< x value of inner hit*/
	float theOuterX; /**< x value of outer hit*/
	float theInnerY; /**< y value of inner hit*/
	float theOuterY; /**< y value of outer hit*/
	float theInnerZ; /**< z value of inner hit*/
	float theOuterZ; /**< z value of outer hit*/
	float theInnerR; /**< r value of inner hit*/
	float theOuterR; /**< r value of outer hit*/

} OMPCACell;

void init(OMPCACell* this, const OMPLayerDoublets* doublets, int layerId, int doubletId, int innerHitId, int outerHitId) {
	this->theInnerHitId = innerHitId;
	this->theOuterHitId = outerHitId;

	this->theDoubletId = doubletId;
	this->theLayerIdInFourLayers = layerId;
	this->theInnerX = doublets->layers[0].p[3*(doublets->indices[2 * doubletId])];
	this->theOuterX = doublets->layers[1].p[3*(doublets->indices[2 * doubletId + 1])];

	this->theInnerY = doublets->layers[0].p[1+3*(doublets->indices[2 * doubletId])];
	this->theOuterY = doublets->layers[1].p[1+3*(doublets->indices[2 * doubletId + 1])];

	this->theInnerZ = doublets->layers[0].p[2+3*(doublets->indices[2 * doubletId])];
	this->theOuterZ = doublets->layers[1].p[2+3*(doublets->indices[2 * doubletId + 1])];
	this->theInnerR = hypot(this->theInnerX, this->theInnerY);
	this->theOuterR = hypot(this->theOuterX, this->theOuterY);
}

/** Get x value of inner hit. Provides an abstaction
* for getting the value.
\param doublets The structure that holds all pairs of layers
\param thisInnerIdx The integer index of the inner hit in the layer pair
\param thisLayer The integer index of the layer pair
*/
float get_inner_x(const OMPLayerDoublets* doublets, const int thisInnerIdx, const int thisLayer) {
	return doublets[thisLayer].layers[0].p[3*thisInnerIdx];
}
/** Get x value of outer hit. Provides an abstaction
* for getting the value.
\param doublets The structure that holds all pairs of layers
\param thisOuterIdx The integer index of the outer hit in the layer pair
\param thisLayer The integer index of the layer pair
*/
float get_outer_x(const OMPLayerDoublets* doublets, const int thisOuterIdx, const int thisLayer) {
	return doublets[thisLayer].layers[1].p[3*thisOuterIdx];
}
/** Get y value of inner hit. Provides an abstaction
* for getting the value.
\param doublets The structure that holds all pairs of layers
\param thisInnerIdx The integer index of the inner hit in the layer pair
\param thisLayer The integer index of the layer pair
*/
float get_inner_y(const OMPLayerDoublets* doublets, const int thisInnerIdx, const int thisLayer) {
	return doublets[thisLayer].layers[0].p[3*thisInnerIdx+1];
}
/** Get y value of outer hit. Provides an abstaction
* for getting the value.
\param doublets The structure that holds all pairs of layers
\param thisOuterIdx The integer index of the outer hit in the layer pair
\param thisLayer The integer index of the layer pair
*/
float get_outer_y(const OMPLayerDoublets* doublets, const int thisOuterIdx, const int thisLayer) {
	return doublets[thisLayer].layers[1].p[3*thisOuterIdx+1];
}
/** Get z value of inner hit. Provides an abstaction
* for getting the value.
\param doublets The structure that holds all pairs of layers
\param thisInnerIdx The integer index of the inner hit in the layer pair
\param thisLayer The integer index of the layer pair
*/
float get_inner_z(const OMPLayerDoublets* doublets, const int thisInnerIdx, const int thisLayer) {
	return doublets[thisLayer].layers[0].p[3*thisInnerIdx+2];
}
/** Get z value of outer hit. Provides an abstaction
* for getting the value.
\param doublets The structure that holds all pairs of layers
\param thisOuterIdx The integer index of the outer hit in the layer pair
\param thisLayer The integer index of the layer pair
*/
float get_outer_z(const OMPLayerDoublets* doublets, const int thisOuterIdx, const int thisLayer) {
	return doublets[thisLayer].layers[1].p[3*thisOuterIdx+2];
}
/** Get r value of inner hit. Provides an abstaction
* for getting the value.
\param doublets The structure that holds all pairs of layers
\param thisInnerIdx The integer index of the inner hit in the layer pair
\param thisLayer The integer index of the layer pair
*/
float get_inner_r(const OMPLayerDoublets* doublets, const int thisIdx, const int thisLayer) {
	return doublets[thisLayer].r[2*thisIdx];
}
/** Get r value of outer hit. Provides an abstaction
* for getting the value.
\param doublets The structure that holds all pairs of layers
\param thisOuterIdx The integer index of the outer hit in the layer pair
\param thisLayer The integer index of the layer pair
*/
float get_outer_r(const OMPLayerDoublets* doublets, const int thisIdx, const int thisLayer) {
	return doublets[thisLayer].r[2*thisIdx+1];
}
/** Get integer index of inner hit. Provides an abstaction
* for getting the index of inner hit.
\param doublets The structure that holds all pairs of layers
\param thisIdx The integer index of the doublet in the layer pair
\param thisLayer The integer index of the layer pair
*/
unsigned int get_inner_hit_id(const OMPLayerDoublets* doublets, const int thisIdx, const int thisLayer) {
	return doublets[thisLayer].indices[2*thisIdx];
}
/** Get integer index of outer hit. Provides an abstaction
* for getting the index of outer hit.
\param doublets The structure that holds all pairs of layers
\param thisIdx The integer index of the doublet in the layer pair
\param thisLayer The integer index of the layer pair
*/
unsigned int get_outer_hit_id(const OMPLayerDoublets* doublets, const int thisIdx, const int thisLayer) {
	return doublets[thisLayer].indices[2*thisIdx+1];
}


int are_aligned_RZ(const OMPLayerDoublets* doublets, const int* this, const int thisLayer, const int* otherCell, const float ptmin, const float thetaCut) {
        float r1 = get_inner_r(doublets, otherCell[0], thisLayer-1);
        float z1 = get_inner_z(doublets, otherCell[1], thisLayer-1);
        float thisOuterR = get_outer_r(doublets, this[0], thisLayer);
        float thisOuterZ = get_outer_z(doublets, this[2], thisLayer);
        float thisInnerR = get_inner_r(doublets, this[0], thisLayer);
        float thisInnerZ = get_inner_z(doublets, this[1], thisLayer);

        float distance_13_squared = (r1 - thisOuterR) * (r1 - thisOuterR) 
        	+ (z1 - thisOuterZ) * (z1 - thisOuterZ);
        float tan_12_13_half = fabs(z1 * (thisInnerR - thisOuterR) + thisInnerZ * (thisOuterR - r1) + thisOuterZ * (r1 - thisInnerR)) / distance_13_squared;
        return tan_12_13_half * ptmin <= thetaCut;
}

int have_similar_curvature(const OMPLayerDoublets* doublets, const int* this, const int thisLayer, const int* otherCell, const float region_origin_x, const float region_origin_y, const float region_origin_radius, const float phiCut) {
	float x1 = get_inner_x(doublets, otherCell[1], thisLayer-1);
	float y1 = get_inner_y(doublets, otherCell[1], thisLayer-1);

	float x2 = get_inner_x(doublets, this[1], thisLayer);
	float y2 = get_inner_y(doublets, this[1], thisLayer);

	float x3 = get_outer_x(doublets, this[2], thisLayer);
	float y3 = get_outer_y(doublets, this[2], thisLayer);

	float precision = 0.5f;
	float offset = x2 * x2 + y2 * y2;

	float bc = (x1 * x1 + y1 * y1 - offset) / 2.f;

	float cd = (offset - x3 * x3 - y3 * y3) / 2.f;

	float det = (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);

                //points are aligned
	if (fabs(det) < precision)
		return 1;

	float idet = 1.f / det;

	float x_center = (bc * (y2 - y3) - cd * (y1 - y2)) * idet;
	float y_center = (cd * (x1 - x2) - bc * (x2 - x3)) * idet;

	float radius = sqrt((x2 - x_center) * (x2 - x_center)+ (y2 - y_center) * (y2 - y_center));
	float centers_distance_squared = (x_center - region_origin_x) * (x_center - region_origin_x) + 
		(y_center - region_origin_y) * (y_center - region_origin_y);

	float minimumOfIntesectionRange = (radius - region_origin_radius) * (radius - region_origin_radius) - phiCut;

	if (centers_distance_squared >= minimumOfIntesectionRange)
	{
		float maximumOfIntesectionRange = (radius + region_origin_radius)
			* (radius + region_origin_radius) + phiCut;
		return centers_distance_squared <= maximumOfIntesectionRange;
	}
	else
	{
		return 0;
	}

}


/** Get boolean value (as integer) of compatibility predicate. Receives 2 cells and checks if they are compatible in respect with the parameters provides.
\param doublets The structure that holds all pairs of layers
\param this The integer index of the outer cell doublet in the layer pair
\param thisLayer The integer index of the layer pair of the outer doublet
\param innerCell The integer index of the inner cell doublet in the layer pair
*/
int check_alignment_and_tag(const OMPLayerDoublets* doublets, const int* this, const int thisLayer, const int* innerCell, const float ptmin, const float region_origin_x, const float region_origin_y, const float region_origin_radius, const float thetaCut, const float phiCut) {

	return (are_aligned_RZ(doublets, this, thisLayer, innerCell, ptmin, thetaCut) && 
		have_similar_curvature(doublets, this, thisLayer, innerCell, region_origin_x, region_origin_y, region_origin_radius, phiCut));

}

#endif /*CACELL_H_ */
