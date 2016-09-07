/** \file test.cu The main functionality for the CUDA program*/

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#include <time.h>

#include "cuda.h"

#include "ca.h"

#include "GPUSimpleVector.h"
#include "GPUCACell.h"
#include "OMPHitsAndDoublets.h"


//! A wrapper struct for holding CA parameters
struct GPUCellularAutomaton {
	float thePtMin;
	float theRegionOriginX;
	float theRegionOriginY;
	float theRegionOriginRadius;
	float theThetaCut;
	float thePhiCut;

	GPUCellularAutomaton (float arg1, float arg2, float arg3, float arg4, float arg5, float arg6) 
		: thePtMin (arg1), theRegionOriginX (arg2), theRegionOriginY (arg3), theRegionOriginRadius (arg4), theThetaCut (arg5), thePhiCut (arg6)
		{

		}
};

//! The function that launches the CUDA kernels
/** This function allocates the memory and runs the CUDA code.
\param grid_size Number of blocks in grid
\param block_size Number of threads in block
\param ca Parameters of the CA
\param doublets The pairs of layers
\param quadruplets The results of the CA
*/
template<unsigned int theNumberOfLayers, unsigned int maxNumberOfQuadruplets>
void run (int grid_size, int block_size,GPUCellularAutomaton& ca, std::array<const GPULayerDoublets *, theNumberOfLayers - 1> const & doublets, std::vector<std::array<int, 4>> & quadruplets) {
	GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* > * hostIsOuterHitOfCell[theNumberOfLayers-2];
    GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* > ** isOuterHitOfCell;

    int numberOfChunksIn1stArena = 0;
    std::array<int, theNumberOfLayers - 1> numberOfKeysIn1stArena;

    for (size_t i = 0; i < theNumberOfLayers - 1; ++i) {
    	numberOfKeysIn1stArena[i] = doublets[i]->layers[1].size;
    	numberOfChunksIn1stArena += doublets[i]->size;
    }
    GPULayerDoublets* gpu_doublets;

    cudaMalloc(&gpu_doublets, 3 * sizeof(GPULayerDoublets));

    for (int i = 0; i < 3; ++i) {
		cudaMemcpy(&gpu_doublets[i], doublets[i], sizeof(GPULayerDoublets), cudaMemcpyHostToDevice);
	}

	cudaMalloc(& hostIsOuterHitOfCell[0], numberOfKeysIn1stArena[0] * sizeof(GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* >) );
	cudaMemset(hostIsOuterHitOfCell[0], 0, numberOfKeysIn1stArena[0] * sizeof(GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* >) );
	cudaMalloc(& hostIsOuterHitOfCell[1], numberOfKeysIn1stArena[1] * sizeof(GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* >) );
	cudaMemset(hostIsOuterHitOfCell[1], 0, numberOfKeysIn1stArena[1] * sizeof(GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* >) );
	cudaMalloc(& isOuterHitOfCell, 2 * sizeof(void*));
	cudaMemcpy(isOuterHitOfCell, hostIsOuterHitOfCell, 2 * sizeof(void*), cudaMemcpyHostToDevice);

	int numberOfChunksIn2ndArena = 0;
	std::array<int, theNumberOfLayers - 2> numberOfKeysIn2ndArena;
	for (size_t i = 1; i < theNumberOfLayers - 1; ++i) {
		numberOfKeysIn2ndArena[i - 1] = doublets[i]->size;
		numberOfChunksIn2ndArena += doublets[i - 1]->size;
	}

	GPUCACell < theNumberOfLayers > **theCells;
	cudaMalloc(&theCells, (theNumberOfLayers - 1) * sizeof(GPUCACell<theNumberOfLayers> *));
	GPUCACell < theNumberOfLayers > *hostCells[theNumberOfLayers - 1];


	for (unsigned int i = 0; i < theNumberOfLayers - 1; ++i)
		cudaMalloc(&hostCells[i],doublets[i]->size* sizeof(GPUCACell<theNumberOfLayers> ));

	cudaMemcpy(theCells, hostCells,
			(theNumberOfLayers - 1) * sizeof(GPUCACell<theNumberOfLayers> *),
			cudaMemcpyHostToDevice);

	GPUSimpleVector<maxNumberOfQuadruplets, int4>* foundNtuplets;
	cudaMallocManaged(&foundNtuplets,sizeof(GPUSimpleVector<maxNumberOfQuadruplets,	int4>));
	cudaMemset(foundNtuplets, 0, sizeof(int));
	
	
	/*create the cells*/
	kernel_create<4><<<dim3(grid_size,3),block_size>>>(gpu_doublets, theCells, isOuterHitOfCell);
	/*connect compatible cells*/
	kernel_connect<4><<<dim3(grid_size,2),block_size>>>(gpu_doublets, theCells, isOuterHitOfCell, ca.thePtMin, ca.theRegionOriginX, ca.theRegionOriginY, ca.theRegionOriginRadius, ca.theThetaCut, ca.thePhiCut);
	/*find the quadrupltes*/
	kernel_find_ntuplets<4,1000><<<16,1024>>>(gpu_doublets, theCells, foundNtuplets, 4);
	cudaDeviceSynchronize();
	
	quadruplets.resize(foundNtuplets->size());
	memcpy(quadruplets.data(), foundNtuplets->m_data, foundNtuplets->size() * sizeof(std::array<int, 4>));

	cudaFree(foundNtuplets);
	for (unsigned int i = 0; i< theNumberOfLayers-1; ++i)
		cudaFree(hostCells[i]);
	cudaFree(theCells);
}

//! \brief A function for loading a track seeding doublets structure from a file
/*! This function retrieves a track seeding doublets structure from a file. It also retrieves some floats that are used for compatibility 
* testing which are also included in the file.File name format is specified in settings.h (e.g. input_[var]_[var].txt)
\param x Identifier of the event (in)
\param y Identifier of the barrel (in)
\param fptr An array of 6 float used in the compatibility check afterwards (out)
\param iptr An array of sizes for each layer of the doublets structure that is sent to the compute cluster processing (out)
\param doublets The data structure for storing the doublets and hits of the problem (out)
*/
void dataLoad (int x, int y, const char* inpath, float* fptr, int* iptr, OMPLayerDoublets* doublets) {
	char filename[100];

	int theNumberOfLayers = 4;

	float data;

	sprintf (filename, inpath, x, y);
	printf ("File: %s\n", filename);
	FILE* fin = fopen(filename, "r");

	int i;

	for (i = 0; i < 6; i++) {
		fscanf (fin, "%f", &data);
		fptr[i] = data;		
	}

	int j,k;

	int numi = 0;
	int nump = 0;

	for (i = 0; i < theNumberOfLayers-1; i++) {
		fscanf (fin, "%d", &(doublets[i].size));
		numi += doublets[i].size;
		iptr[3*i] = doublets[i].size;
		doublets[i].indices = (int*) malloc(2*doublets[i].size*sizeof(int));
		for (j = 0; j < 2*doublets[i].size; j++)
			fscanf (fin, "%d", &(doublets[i].indices[j]));
		
		for (k = 0; k < 2; k++) {
			fscanf (fin, "%d", &(doublets[i].layers[k].size));
			nump += doublets[i].layers[k].size;
			iptr[3*i+k+1] = doublets[i].layers[k].size;
			doublets[i].layers[k].x = (float*) malloc(doublets[i].layers[k].size*sizeof(float));
			doublets[i].layers[k].y = (float*) malloc(doublets[i].layers[k].size*sizeof(float));
			doublets[i].layers[k].z = (float*) malloc(doublets[i].layers[k].size*sizeof(float));			

			for (j = 0; j < doublets[i].layers[k].size; j++) {
				fscanf (fin, "%f", &(doublets[i].layers[k].x[j]));
				fscanf (fin, "%f", &(doublets[i].layers[k].y[j]));
				fscanf (fin, "%f", &(doublets[i].layers[k].z[j]));
			}
		}
	}

	fclose (fin);

}

//! Reads input files and calls the CA algorithm for every one of the problems.
int main (int argc, char** argv) {
	const char* base = "input/input_%d_%d.txt";
	float* fargs[700];
	int* iargs[700];
	OMPLayerDoublets* doubletz[700];

	int theNumberOfLayers = 4;

	cudaSetDevice(1);

	for (int i = 0; i < 700; i++) {
		fargs[i] = (float*) malloc (6*sizeof(float));
		iargs[i] = (int*) malloc (3*(theNumberOfLayers-1)*sizeof(int));
		doubletz[i] = (OMPLayerDoublets*) malloc ((theNumberOfLayers-1)*sizeof(OMPLayerDoublets));
		dataLoad (i/7,i%7,base,fargs[i],iargs[i], doubletz[i]);
	}
	
	clock_t start = clock ();
	for (int i = 0; i < 700; i++) {
		GPULayerDoublets doublets[theNumberOfLayers-1][1];

		for (int j = 0; j < theNumberOfLayers-1; j++) {
			doublets[j]->size = doubletz[i][j].size;
			cudaMalloc (&(doublets[j]->indices), 2*doubletz[i][j].size*sizeof(int));
			cudaMemcpy (doublets[j]->indices, doubletz[i][j].indices, 2*doubletz[i][j].size*sizeof(int), cudaMemcpyHostToDevice);
		

			for (int k = 0; k < 2; k++) {
				doublets[j]->layers[k].size = doubletz[i][j].layers[k].size;
				cudaMalloc (&(doublets[j]->layers[k].x), doubletz[i][j].layers[k].size*sizeof(float));
				cudaMemcpy (doublets[j]->layers[k].x, doubletz[i][j].layers[k].x, doubletz[i][j].layers[k].size*sizeof(float), cudaMemcpyHostToDevice);
				cudaMalloc (&(doublets[j]->layers[k].y), doubletz[i][j].layers[k].size*sizeof(float));
				cudaMemcpy (doublets[j]->layers[k].y, doubletz[i][j].layers[k].y, doubletz[i][j].layers[k].size*sizeof(float), cudaMemcpyHostToDevice);
				cudaMalloc (&(doublets[j]->layers[k].z), doubletz[i][j].layers[k].size*sizeof(float));
				cudaMemcpy (doublets[j]->layers[k].z, doubletz[i][j].layers[k].z, doubletz[i][j].layers[k].size*sizeof(float), cudaMemcpyHostToDevice);
			}
		}

		std::array<const GPULayerDoublets*, 3> doublets_c;
		for (int j = 0; j < theNumberOfLayers-1; j++)
			doublets_c[j] = doublets[j];
		std::vector<std::array<int, 4>> quadruplets;

		GPUCellularAutomaton ca (fargs[i][0], fargs[i][1], fargs[i][2], fargs[i][3], fargs[i][4], fargs[i][5]);
		run <4,1000> (atoi(argv[1]), atoi(argv[2]), ca, doublets_c, quadruplets);

		printf ("%d: Size of results is %d\n", i, (int) quadruplets.size());
	}
	

	clock_t end = clock();
	
	printf ("Total execution time: %f\n", (((double)(end-start))/CLOCKS_PER_SEC));
	
	return 0;
}
