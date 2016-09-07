/** \file test.cu*/

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

/** Structure holding essential parameters*/
struct GPUCellularAutomaton {
	float thePtMin;
	float theRegionOriginX;
	float theRegionOriginY;
	float theRegionOriginRadius;
	float theThetaCut;
	float thePhiCut;

	GPUCellularAutomaton () {}

	GPUCellularAutomaton (float arg1, float arg2, float arg3, float arg4, float arg5, float arg6) 
		: thePtMin (arg1), theRegionOriginX (arg2), theRegionOriginY (arg3), theRegionOriginRadius (arg4), theThetaCut (arg5), thePhiCut (arg6)
		{

		}
};

/** Main work is done here. The function gets 2 problems at once and alternates between the two devices on loading the data,
* executing the CA and getting the data back. Thus, it is executed in parallel in two GPUs.
*/
template<unsigned int theNumberOfLayers, unsigned int maxNumberOfQuadruplets>
void run (int ndevice, int grid_size, int block_size, GPUCellularAutomaton ca[2], std::array<const GPULayerDoublets *, theNumberOfLayers - 1> const doublets[2], GPUSimpleVector<maxNumberOfQuadruplets, int4>* foundNtuplets[2]) {

	GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* > * hostIsOuterHitOfCell[2][theNumberOfLayers-2];
	GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* > ** isOuterHitOfCell[2];

	int numberOfChunksIn1stArena[2];
	std::array<int, theNumberOfLayers - 1> numberOfKeysIn1stArena[2];
	GPULayerDoublets* gpu_doublets[2];
	int numberOfChunksIn2ndArena[2];
	std::array<int, theNumberOfLayers - 2> numberOfKeysIn2ndArena[2];
	GPUCACell < theNumberOfLayers > **theCells[2];
	GPUCACell < theNumberOfLayers > *hostCells[2][theNumberOfLayers - 1];

	for (int device = 0; device < 2; device++) {
		cudaSetDevice(device);
		numberOfChunksIn1stArena[device] = 0;
		for (size_t i = 0; i < theNumberOfLayers - 1; ++i) {
			numberOfKeysIn1stArena[device][i] = doublets[device][i]->layers[1].size;
    			numberOfChunksIn1stArena[device] += doublets[device][i]->size;
    		}
    
    		cudaMalloc(&gpu_doublets[device], 3 * sizeof(GPULayerDoublets));
		
    		for (int i = 0; i < 3; ++i) {
			cudaMemcpy(&gpu_doublets[device][i], doublets[device][i], sizeof(GPULayerDoublets), cudaMemcpyHostToDevice);
		}

		cudaMalloc(& hostIsOuterHitOfCell[device][0], numberOfKeysIn1stArena[device][0] 
					* sizeof(GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* >) );
		cudaMemset(hostIsOuterHitOfCell[device][0], 0, numberOfKeysIn1stArena[device][0] 
					* sizeof(GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* >) );
		cudaMalloc(& hostIsOuterHitOfCell[device][1], numberOfKeysIn1stArena[device][1] 
					* sizeof(GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* >) );
		cudaMemset(hostIsOuterHitOfCell[device][1], 0, numberOfKeysIn1stArena[device][1] 
					* sizeof(GPUSimpleVector<64, GPUCACell<theNumberOfLayers>* >) );
		cudaMalloc(& isOuterHitOfCell[device], 2 * sizeof(void*));
		cudaMemcpy(isOuterHitOfCell[device], hostIsOuterHitOfCell[device], 2 * sizeof(void*), cudaMemcpyHostToDevice);

		numberOfChunksIn2ndArena[device] = 0;	
		for (size_t i = 1; i < theNumberOfLayers - 1; ++i) {
			numberOfKeysIn2ndArena[device][i - 1] = doublets[device][i]->size;
			numberOfChunksIn2ndArena[device] += doublets[device][i - 1]->size;
		}


		cudaMalloc(&theCells[device], (theNumberOfLayers - 1) * sizeof(GPUCACell<theNumberOfLayers> *));

		for (unsigned int i = 0; i < theNumberOfLayers - 1; ++i) {
			cudaMalloc(&hostCells[device][i],doublets[device][i]->size* sizeof(GPUCACell<theNumberOfLayers> ));
		}

		cudaMemcpy(theCells[device], hostCells[device], (theNumberOfLayers - 1) 
			* sizeof(GPUCACell<theNumberOfLayers> *), cudaMemcpyHostToDevice);

		cudaMemset(foundNtuplets[device], 0, sizeof(int));
	}
	cudaSetDevice(0);

	kernel_create<4><<<dim3(grid_size,3),block_size>>>(gpu_doublets[0], theCells[0], isOuterHitOfCell[0]);

	cudaSetDevice(1);

	kernel_create<4><<<dim3(grid_size,3),block_size>>>(gpu_doublets[1], theCells[1], isOuterHitOfCell[1]);

	cudaSetDevice(0);

	kernel_connect<4><<<dim3(grid_size,2),block_size>>>(gpu_doublets[0], theCells[0], isOuterHitOfCell[0], ca[0].thePtMin, ca[0].theRegionOriginX, ca[0].theRegionOriginY, ca[0].theRegionOriginRadius, ca[0].theThetaCut, ca[0].thePhiCut);
	
	cudaSetDevice(1);

        kernel_connect<4><<<dim3(grid_size,2),block_size>>>(gpu_doublets[1], theCells[1], isOuterHitOfCell[1], ca[1].thePtMin, ca[1].theRegionOriginX, ca[1].theRegionOriginY, ca[1].theRegionOriginRadius, ca[1].theThetaCut, ca[1].thePhiCut);

	cudaSetDevice(0);

	kernel_find_ntuplets<4,1000><<<16,1024>>>(gpu_doublets[0], theCells[0], foundNtuplets[0], 4);

	cudaSetDevice(1);

        kernel_find_ntuplets<4,1000><<<16,1024>>>(gpu_doublets[1], theCells[1], foundNtuplets[1], 4);


	cudaSetDevice(0);

	cudaDeviceSynchronize ();

	cudaSetDevice(1);

	cudaDeviceSynchronize ();

	for (unsigned int i = 0; i< theNumberOfLayers-1; ++i) {
		cudaFree(hostCells[0][i]);
		cudaFree(hostCells[1][i]);
	}
	cudaFree(theCells[0]);
	cudaFree(theCells[1]);
}

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

int main (int argc, char** argv) {
	const char* base = "input/input_%d_%d.txt";
	float* fargs[700];
	int* iargs[700];
	OMPLayerDoublets* doubletz[700];

	int nDevices;
	cudaGetDeviceCount(&nDevices);
	printf ("%d CUDA devices found\n", nDevices);
	for (int i = 0; i < nDevices; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties (&prop, i);
		printf ("%d->%s\n", i, prop.name);
	}
	GPUSimpleVector<1000, int4> found;
	GPUSimpleVector<1000, int4>* foundNtuplets[2];
	for (int i = 0; i < nDevices; i++) {
		cudaSetDevice(i);
		cudaMallocManaged(&foundNtuplets[i],sizeof(GPUSimpleVector<1000,int4>));
	}
	int theNumberOfLayers = 4;

	for (int i = 0; i < 700; i++) {
		fargs[i] = (float*) malloc (6*sizeof(float));
		iargs[i] = (int*) malloc (3*(theNumberOfLayers-1)*sizeof(int));
		doubletz[i] = (OMPLayerDoublets*) malloc ((theNumberOfLayers-1)*sizeof(OMPLayerDoublets));
		dataLoad (i/7,i%7,base,fargs[i],iargs[i], doubletz[i]);
	}
	
	clock_t start = clock ();
	for (int i = 0; i < 350; i++) {
		/*load data on the two devices*/
		GPULayerDoublets doublets[2][theNumberOfLayers-1][1];

		cudaSetDevice(0);

		for (int j = 0; j < theNumberOfLayers-1; j++) {
			doublets[0][j]->size = doubletz[2*i][j].size;
			cudaMalloc (&(doublets[0][j]->indices), 2*doubletz[2*i][j].size*sizeof(int));
			cudaMemcpy (doublets[0][j]->indices, doubletz[2*i][j].indices, 2*doubletz[2*i][j].size*sizeof(int), cudaMemcpyHostToDevice);
		

			for (int k = 0; k < 2; k++) {
				doublets[0][j]->layers[k].size = doubletz[2*i][j].layers[k].size;
				cudaMalloc (&(doublets[0][j]->layers[k].x), doubletz[2*i][j].layers[k].size*sizeof(float));
				cudaMemcpy (doublets[0][j]->layers[k].x, doubletz[2*i][j].layers[k].x, doubletz[2*i][j].layers[k].size*sizeof(float), cudaMemcpyHostToDevice);
				cudaMalloc (&(doublets[0][j]->layers[k].y), doubletz[2*i][j].layers[k].size*sizeof(float));
				cudaMemcpy (doublets[0][j]->layers[k].y, doubletz[2*i][j].layers[k].y, doubletz[2*i][j].layers[k].size*sizeof(float), cudaMemcpyHostToDevice);
				cudaMalloc (&(doublets[0][j]->layers[k].z), doubletz[2*i][j].layers[k].size*sizeof(float));
				cudaMemcpy (doublets[0][j]->layers[k].z, doubletz[2*i][j].layers[k].z, doubletz[2*i][j].layers[k].size*sizeof(float), cudaMemcpyHostToDevice);
			}
		}

		cudaSetDevice(1);

		for (int j = 0; j < theNumberOfLayers-1; j++) {
			doublets[1][j]->size = doubletz[2*i+1][j].size;
			cudaMalloc (&(doublets[1][j]->indices), 2*doubletz[2*i+1][j].size*sizeof(int));
			cudaMemcpy (doublets[1][j]->indices, doubletz[2*i+1][j].indices, 2*doubletz[2*i+1][j].size*sizeof(int), cudaMemcpyHostToDevice);
		

			for (int k = 0; k < 2; k++) {
				doublets[1][j]->layers[k].size = doubletz[2*i+1][j].layers[k].size;
				cudaMalloc (&(doublets[1][j]->layers[k].x), doubletz[2*i+1][j].layers[k].size*sizeof(float));
				cudaMemcpy (doublets[1][j]->layers[k].x, doubletz[2*i+1][j].layers[k].x, doubletz[2*i+1][j].layers[k].size*sizeof(float), cudaMemcpyHostToDevice);
				cudaMalloc (&(doublets[1][j]->layers[k].y), doubletz[2*i+1][j].layers[k].size*sizeof(float));
				cudaMemcpy (doublets[1][j]->layers[k].y, doubletz[2*i+1][j].layers[k].y, doubletz[2*i+1][j].layers[k].size*sizeof(float), cudaMemcpyHostToDevice);
				cudaMalloc (&(doublets[1][j]->layers[k].z), doubletz[2*i+1][j].layers[k].size*sizeof(float));
				cudaMemcpy (doublets[1][j]->layers[k].z, doubletz[2*i+1][j].layers[k].z, doubletz[2*i+1][j].layers[k].size*sizeof(float), cudaMemcpyHostToDevice);
			}
		}

		std::array<const GPULayerDoublets*, 3> doublets_c[2];
		for (int j = 0; j < theNumberOfLayers-1; j++)
			doublets_c[0][j] = doublets[0][j];
		for (int j = 0; j < theNumberOfLayers-1; j++)
			doublets_c[1][j] = doublets[1][j];
		std::vector<std::array<int, 4>> quadruplets;

		GPUCellularAutomaton ca1 (fargs[2*i][0], fargs[2*i][1], fargs[2*i][2], fargs[2*i][3], fargs[2*i][4], fargs[2*i][5]);
		GPUCellularAutomaton ca2 (fargs[2*i+1][0], fargs[2*i+1][1], fargs[2*i+1][2], fargs[2*i+1][3], fargs[2*i+1][4], fargs[2*i+1][5]);
		GPUCellularAutomaton ca[2];
		ca[0] = ca1;
		ca[1] = ca2;

		run <4,1000> (nDevices, atoi(argv[1]), atoi(argv[2]), ca, doublets_c, foundNtuplets);

		cudaSetDevice(0);
		cudaMemcpy (&found, foundNtuplets[0], sizeof(GPUSimpleVector<1000,int4>), cudaMemcpyDeviceToHost);
		quadruplets.resize(found.size());
	        memcpy(quadruplets.data(), found.m_data, found.size() * sizeof(std::array<int, 4>));
		printf ("%d: Size of results is %d\n", 2*i, (int) quadruplets.size());
		cudaSetDevice(1);
                cudaMemcpy (&found, foundNtuplets[1], sizeof(GPUSimpleVector<1000,int4>), cudaMemcpyDeviceToHost);
                quadruplets.resize(found.size());
                memcpy(quadruplets.data(), found.m_data, found.size() * sizeof(std::array<int, 4>));
                printf ("%d: Size of results is %d\n", 2*i+1, (int) quadruplets.size());
	}
	
	
	clock_t end = clock();
	
	printf ("Total execution time: %f\n", (((double)(end-start))/CLOCKS_PER_SEC));
	
	return 0;
}
