/** \file GPUHitsAndDoublets.h*/


#ifndef RecoPixelVertexing_PixelTriplets_GPUHitsAndDoublets_h
#define RecoPixelVertexing_PixelTriplets_GPUHitsAndDoublets_h

//! The structure for a layer
/** This structure is used for representing each layer of the detector
* Basically it is a vector of the hits in the layer
*/
struct GPULayerHits
{
	size_t size; /**< The number of hits in this layer*/
	float * x;
	float * y;
	float * z;
};

//! The structure that holds a pair of layers
/** This structure holds a pair of layer structures as previously described,
* and a vector of the doublets formed by the hits.
*/
struct GPULayerDoublets
{
	size_t size; /**< The number of doublets in the pair of layers*/
	int * indices; /**< The indices of hits in doublets. For the i-th doublet, indices[2i] is inner hit and indices[2i+1] is outer hit*/
	GPULayerHits layers[2]; /**< the pair of layers in question*/
};



#endif // not defined RecoPixelVertexing_PixelTriplets_GPUHitsAndDoublets_h
