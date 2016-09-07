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
/** Keeps all the relevant info for compatibility checks.
*/
typedef struct OMPCACell
{	
	unsigned int theInnerHitId; /**< integer index of inner hit in its layer*/
	unsigned int theOuterHitId; /**< integer index of outer hit in its layer*/
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


//! Initialize the values of the cell.
/** The cells correspond 1-1 to a given doublet.
\param this The cell in question
\param doublets The structure of the layer pair
\param layerId The layer in which the doublet resides
\param innerHitId The identifier of the inner hit in its layer
\param outerHitId The identifier of the outer hit in its layer
*/
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
\param this Pointer to the cell in question
*/
float get_inner_x(const OMPCACell* this) {
	return this->theInnerX;
}
/** Get x value of outer hit. Provides an abstaction
* for getting the value.
\param this Pointer to the cell in question
*/
float get_outer_x(const OMPCACell* this) {
	return this->theOuterX;
}
/** Get y value of inner hit. Provides an abstaction
* for getting the value.
\param this Pointer to the cell in question
*/
float get_inner_y(const OMPCACell* this) {
	return this->theInnerY;
}
/** Get y value of outer hit. Provides an abstaction
* for getting the value.
\param this Pointer to the cell in question
*/
float get_outer_y(const OMPCACell* this) {
	return this->theOuterY;
}
/** Get z value of inner hit. Provides an abstaction
* for getting the value.
\param this Pointer to the cell in question
*/
float get_inner_z(const OMPCACell* this) {
	return this->theInnerZ;
}
/** Get z value of outer hit. Provides an abstaction
* for getting the value.
\param this Pointer to the cell in question
*/
float get_outer_z(const OMPCACell* this) {
	return this->theOuterZ;
}
/** Get r value of inner hit. Provides an abstaction
* for getting the value.
\param this Pointer to the cell in question
*/
float get_inner_r(const OMPCACell* this) {
	return this->theInnerR;
}
/** Get r value of outer hit. Provides an abstaction
* for getting the value.
\param this Pointer to the cell in question
*/
float get_outer_r(const OMPCACell* this) {
	return this->theOuterR;
}
/** Get id of inner hit. Provides an abstaction
* for getting the value.
\param this Pointer to the cell in question
*/
unsigned int get_inner_hit_id(const OMPCACell* this) {
	return this->theInnerHitId;
}
/** Get id of outer hit. Provides an abstaction
* for getting the value.
\param this Pointer to the cell in question
*/
unsigned int get_outer_hit_id(const OMPCACell* this) {
	return this->theOuterHitId;
}

void print_cell(const OMPCACell* this) {

		printf(
				"\n printing cell: %d, on layer: %d, innerHitId: %d, outerHitId: %d, innerradius %f, outerRadius %f ",
				this->theDoubletId, this->theLayerIdInFourLayers, this->theInnerHitId,
				this->theOuterHitId, this->theInnerR, this->theOuterR);

}

int are_aligned_RZ(const OMPCACell* this, const OMPCACell* otherCell, const float ptmin, const float thetaCut) {
        float r1 = get_inner_r(otherCell);
        float z1 = get_inner_z(otherCell);
        float distance_13_squared = (r1 - this->theOuterR) * (r1 - this->theOuterR) + (z1 - this->theOuterZ) * (z1 - this->theOuterZ);
        float tan_12_13_half = fabs(z1 * (this->theInnerR - this->theOuterR) + this->theInnerZ * (this->theOuterR - r1) + this->theOuterZ * (r1 - this->theInnerR)) / distance_13_squared;
        return tan_12_13_half * ptmin <= thetaCut;
}

int have_similar_curvature(const OMPCACell* this, const OMPCACell* otherCell, const float region_origin_x, const float region_origin_y, const float region_origin_radius, const float phiCut) {
                float x1 = get_inner_x(otherCell);
                float y1 = get_inner_y(otherCell);

                float x2 = get_inner_x(this);
                float y2 = get_inner_y(this);

                float x3 = get_outer_x(this);
                float y3 = get_outer_y(this);

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
                float centers_distance_squared = (x_center - region_origin_x) * (x_center - region_origin_x) + (y_center - region_origin_y) * (y_center - region_origin_y);

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
\param this The integer index of the outer cell doublet in the layer pair
\param innerCell The integer index of the inner cell doublet in the layer pair
*/
int check_alignment_and_tag(const OMPCACell* this, const OMPCACell* innerCell, const float ptmin, const float region_origin_x, const float region_origin_y, const float region_origin_radius, const float thetaCut, const float phiCut) {

	return (are_aligned_RZ(this, innerCell, ptmin, thetaCut) && have_similar_curvature(this, innerCell, region_origin_x, region_origin_y, region_origin_radius, phiCut));

}


#endif /*CACELL_H_ */
