/** \file main.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <omp.h>

#include "shared_defs.h"
#include "settings.h"
#include "OMPResultVector.h"
#include "OMPHitsAndDoublets.h"
#include "OMPCACell.h"
#include "OMPSimpleVector.h"
#include "OMPMixedList.h"

//! Retrieves quadruplets out of connected doublets
/*! This function gets the connected doublets and performs a DFS to produce the quadruplets.
\param doublets The data structure of the hit and doublet data
\param layer The current layer in the search
\param idx The identifier within the layer for the current doublet we are visiting
\param ml The lists of the connections between doublets
\param top An array that holds the indices of the start of every successor list
\param foundNtuplets The result returned for our problem
\param tmpNtuplet The current path of the search
\param minHitsPerNtuplet The number of layers involved
*/
void find_ntuplets(OMPCACell** cell, unsigned int layer, unsigned int idx, MixedList* ml, int** top, 
			OMPResultVector* foundNtuplets, OMPSimpleVector* tmpNtuplet, const unsigned int minHitsPerNtuplet) {
	int j;

	int4 found;
	const OMPCACell* otherCell;

	if (layer == 0) {
		if (sizesv(tmpNtuplet) >= minHitsPerNtuplet - 1) {
			found.elem[0]=get_inner_hit_id(tmpNtuplet->m_data[2]);
			found.elem[1]=get_outer_hit_id(tmpNtuplet->m_data[2]);
			found.elem[2]=get_outer_hit_id(tmpNtuplet->m_data[1]);
			found.elem[3]=get_outer_hit_id(tmpNtuplet->m_data[0]);
			push_backtsrv(foundNtuplets, &found);
		}else
			return;
	} else {
		int ptr = top[layer][idx];
		while (ptr != -1) {
			int otherIdx = fetch_ml (ml, ptr);
			push_backsv(tmpNtuplet, &cell[layer-1][otherIdx]);
			find_ntuplets(cell, layer-1, otherIdx, ml, top, foundNtuplets, tmpNtuplet, minHitsPerNtuplet);
			pop_backsv(tmpNtuplet, &otherCell);
			ptr = next_ml (ml, ptr);
		}		
	}
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
void dataLoad (int x, int y, float* fptr, int* iptr, OMPLayerDoublets* doublets) {
	char filename[100];

	float data;

	sprintf (filename, inpath, x, y);
	printf ("File: %s\n", filename);
	FILE* fin = fopen(filename, "r");
	if (fin == NULL)
		printf ("wtf\n");

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
			doublets[i].layers[k].p = (float*) malloc(3*doublets[i].layers[k].size*sizeof(float));
			
			for (j = 0; j < doublets[i].layers[k].size; j++) {
				fscanf (fin, "%f", &(doublets[i].layers[k].p[3*j]));
				fscanf (fin, "%f", &(doublets[i].layers[k].p[3*j+1]));
				fscanf (fin, "%f", &(doublets[i].layers[k].p[3*j+2]));
			}
		}
	}
	fclose (fin);

}

//! Main work of CA is done in the main function
int
main(int argc __attribute__ ((unused)), char *argv[])
{
	int loop;
	int i;
	int j;
	int thread_num = atoi(argv[1]);

	OMPLayerDoublets* doubletss[700];
	float* fptr[700];
	int* iptr[700];

	for (i = 0; i < 100; i++) 
		for (j = 0; j < 7; j++) {
			doubletss[7*i+j] = malloc((theNumberOfLayers-1)*sizeof(OMPLayerDoublets));
			fptr[7*i+j] = malloc(6*sizeof(float));
			iptr[7*i+j] = malloc(3*(theNumberOfLayers-1)*sizeof(int));
			dataLoad (i, j, fptr[7*i+j], iptr[7*i+j], doubletss[7*i+j]);
		}

	for (loop = 0; loop < 700; loop++) {
	OMPLayerDoublets* doublets = doubletss[loop];

	float* fargs = fptr[loop];

	OMPCACell* cell[3];
	OMPSimpleVector* isOuterHitOfCell[2];
	int size[3];
	int sizeOuter[2];
	int* liptr[3];

	for (i = 0; i < theNumberOfLayers-1; i++) {
		size[i] = doublets[i].size;
		cell[i] = malloc(size[i]*sizeof(OMPCACell));
		if (i < 2) {
			sizeOuter[i] = doublets[i].layers[1].size;
			isOuterHitOfCell[i] = malloc(sizeOuter[i]*sizeof(OMPSimpleVector));
			memset (isOuterHitOfCell[i], 0, sizeOuter[i]*sizeof(OMPSimpleVector));
		}
	}


	int outer;
	/*create phase*/
	for (i = 0; i < theNumberOfLayers-1; i++) {
		#pragma omp parallel for num_threads(thread_num) private(outer)
		for (j = 0; j < size[i]; j++) {
			init(&(cell[i][j]),&doublets[i],i,j,doublets[i].indices[2*j],doublets[i].indices[2*j+1]);
			if(i < 2) {
				outer = get_outer_hit_id(&cell[i][j]);
				push_back_tssv(&(isOuterHitOfCell[i][outer]), &(cell[i][j]));
			}	
		}
	}




	MixedList* ml = malloc(sizeof(MixedList));
	init_ml(ml);

	int layer;
	unsigned int s;
	int top, inner, val;
	const OMPCACell* otherCell;
	
	for (layer = 1; layer < theNumberOfLayers-1; layer++) {
		s = size[layer];
		liptr[layer] = malloc(s*sizeof(int));
	}

	/*
	connect phase of the algorithm
	find the alignments of doublets
	*/
	for (layer = 1; layer < theNumberOfLayers-1; layer++) {
		s = size[layer];
		#pragma omp parallel for num_threads(thread_num) private(j,top,inner,otherCell, val)
		for (i = 0; i < s; i++) {
			top = -1;
			inner = get_inner_hit_id(&cell[layer][i]);
			for (j = 0; j < sizesv(&isOuterHitOfCell[layer-1][inner]); ++j) {
				otherCell
					= isOuterHitOfCell[layer-1][inner].m_data[j];

        	   		if (check_alignment_and_tag(&(cell[layer][i]),otherCell,
								fargs[0], fargs[1], fargs[2],
								fargs[3], fargs[4], fargs[5])) {
        	   			//push_back_tssv(&(cell[layer][i].theInnerNeighbors), otherCell);
						val = otherCell->theDoubletId;
						top = push_back_ml (ml, top, val);
						if (top < 0) {
							break;
						}
				}
			}
			liptr[layer][i] = top;
		}
	}

	unsigned int lastLayerPairIndex = numberOfLayers - 2;
	s = size[lastLayerPairIndex];

	OMPResultVector results;
	resetrv(&results);
	/*
	DFS for getting results
	*/
	OMPSimpleVector stack;

	#pragma omp parallel for num_threads(thread_num) private(stack)
	for (i = 0; i < s; i++)
	{
		resetsv(&stack);
		push_backsv(&stack, &cell[lastLayerPairIndex][i]);
		find_ntuplets(cell, lastLayerPairIndex, i, ml, liptr, &results, &stack, 4);
	}
	
	int n = sizerv(&results);
	printf ("%d: Size of results is %d\n", loop, n);;

	}

	return 0;
}
